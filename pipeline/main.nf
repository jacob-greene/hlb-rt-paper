#!/usr/bin/env nextflow
// Replication timing of replication-dependent histone genes across ENCODE cell lines (hg19).
//
// ENCODE wavelet-smoothed RT bigWigs (108 tracks, 50 cell-line labels)
//   -> per-gene RT score (TSS +/- 5 kb)
// ENCODE RNA-seq gene quantifications (104 files, 24 labels)
//   -> top 5% expressed protein-coding genes per label
// Drug-treatment Repli-seq fractions (Zenodo 10.5281/zenodo.19093559; 4 conditions x 6 fractions)
//   -> weighted-average J7 wavelet track per condition -> per-gene RT score (exact gene body)
// All -> notebook inputs -> analysis.ipynb (sections 1-5) -> figures.
nextflow.enable.dsl = 2

include { FETCH as FETCH_REF; FETCH as FETCH_BW; FETCH as FETCH_RNA; FETCH as FETCH_UW; FETCH as FETCH_DRUG } from './modules/fetch'
include { GENESETS; GENE_RT; RNASEQ_TOP5; EXTRACT_DRUG; WAVELET_DRUG; DRUG_GENE_RT; ASSEMBLE; NOTEBOOK; REPROCESS_UW; QC_TRACKS } from './modules/analysis'

def subset(rows, key, wanted) {
    wanted == 'all' ? rows : rows.findAll { it[key] in wanted.tokenize(',') }
}

workflow {
    def gtf_name = params.gencode_url.tokenize('/').last()
    FETCH_REF(Channel.of([gtf_name.tokenize('.')[1], 'ref', params.gencode_url, params.gencode_md5, gtf_name]))
    GENESETS(FETCH_REF.out)

    def hg19_sheet = file(params.samples_hg19)
    def bw_rows = subset(hg19_sheet.splitCsv(header: true), 'cell_line', params.cell_lines)
    FETCH_BW(Channel.fromList(bw_rows).map { r ->
        [r.file_accession, 'bigwig_hg19', r.url, r.md5sum, "${r.file_accession}.bigWig"] })
    GENE_RT(FETCH_BW.out, GENESETS.out.coding.first())   // value channel: reused for every track

    def rna_sheet = file(params.samples_rnaseq)
    def rna_rows = subset(rna_sheet.splitCsv(header: true), 'cell_line', params.cell_lines)
    FETCH_RNA(Channel.fromList(rna_rows).map { r ->
        [r.file_accession, 'rnaseq', r.url, r.md5sum, "${r.file_accession}.tsv"] })
    RNASEQ_TOP5(rna_sheet, FETCH_RNA.out.map { it[1] }.collect(), GENESETS.out.coding)

    // Drug-treatment Repli-seq: the six fraction bigWigs of each condition, from the authors' Zenodo archive
    def fraction_order = ['G1b', 'S1', 'S2', 'S3', 'S4', 'G2']
    def condition_order = ['NT', 'ATRi', 'Aph06', 'ATRi+Aph06']
    def drug_sheet = file(params.samples_drug)
    def drug_rows = drug_sheet.splitCsv(header: true)
    def archive_name = params.drug_archive_url.tokenize('/')[-2]
    FETCH_DRUG(Channel.of(['zenodo', 'drug_archive', params.drug_archive_url, params.drug_archive_md5, archive_name]))
    EXTRACT_DRUG(FETCH_DRUG.out.map { it[1] }, drug_sheet)
    def drug_groups = EXTRACT_DRUG.out.flatMap { files ->
        def by_name = files.collectEntries { [it.name, it] }
        drug_rows.groupBy { it.condition }.collect { cond, rows ->
            [cond, "${cond.replace('+', '-')}_rep${rows[0].replicate}",
             rows.sort { fraction_order.indexOf(it.fraction) }.collect { by_name[it.file_name] }] }
    }
    WAVELET_DRUG(drug_groups)
    DRUG_GENE_RT(WAVELET_DRUG.out.toSortedList { a, b -> condition_order.indexOf(a[0]) <=> condition_order.indexOf(b[0]) },
                 WAVELET_DRUG.out.map { it[1] }.collect(), GENESETS.out.coding, GENESETS.out.histones)

    ASSEMBLE(hg19_sheet, GENESETS.out.histones, RNASEQ_TOP5.out, DRUG_GENE_RT.out, GENE_RT.out.collect())
    NOTEBOOK(ASSEMBLE.out.inputs)   // analysis.ipynb, sections 1-5

    // QC (not used by the notebook): rebuild each UW Repli-seq track from its six fraction signals
    // (UCSC ENCODE archive) and compare it with the released ENCODE track.
    if (params.qc) {
        def order = ['G1b', 'S1', 'S2', 'S3', 'S4', 'G2']
        def uw_rows = file(params.samples_uw_fractions).splitCsv(header: true)
                          .findAll { params.cell_lines == 'all' || it.cell_line in params.cell_lines.tokenize(',') }
        def uw_label = uw_rows.collectEntries { [it.wavelet_accession, it.cell_line] }
        FETCH_UW(Channel.fromList(uw_rows).map { r ->
            ["${r.wavelet_accession}__${order.indexOf(r.fraction)}", 'uw_fractions', r.url, r.md5sum, r.file_name] })
        def grouped = FETCH_UW.out
            .map { id, f -> def (acc, i) = id.tokenize('__'); [acc, i as int, f] }
            .groupTuple(size: 6)
            .map { acc, idx, files -> [acc, uw_label[acc], [idx, files].transpose().sort { it[0] }.collect { it[1] }] }
        REPROCESS_UW(grouped)
        def pairs = REPROCESS_UW.out.join(FETCH_BW.out)        // [acc, label, rebuilt, encode]
            .toList()
        def label_of = { acc, label, n -> n > 1 ? "${label} (${acc})" : label }
        def pairs_csv = pairs.map { rows ->
            def counts = rows.countBy { it[1] }
            def f = file("${workflow.workDir}/uw_pairs.csv")
            f.text = 'label,encode,reprocessed\n' + rows.collect { r ->
                "\"${label_of(r[0], r[1], counts[r[1]])}\",${r[3].name},${r[2].name}" }.join('\n') + '\n'
            f }
        QC_TRACKS(pairs_csv, pairs.map { rows -> rows.collectMany { [it[2], it[3]] } })
    }
}
