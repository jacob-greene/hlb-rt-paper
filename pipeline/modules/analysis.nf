// Python steps. Every script lives in ../scripts and runs in the uv environment
// given by params.python (default: the repository's .venv).

process GENESETS {
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(id), path(gtf)

    output:
    path 'genesets/all_coding.parquet', emit: coding
    path 'genesets/histones.parquet',   emit: histones
    path 'genesets/*.csv'

    script:
    """
    ${params.python} ${projectDir}/../scripts/build_genesets.py ${gtf} genesets
    """
}

process GENE_RT {
    tag "${id}"
    label 'small'
    publishDir "${params.outdir}/gene_rt", mode: 'copy'

    input:
    tuple val(id), path(track)
    path coding

    output:
    path "${id}.parquet"

    script:
    """
    ${params.python} ${projectDir}/../scripts/gene_rt.py ${track} ${coding} ${id}.parquet
    """
}

process RNASEQ_TOP5 {
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sheet
    path tsvs
    path coding

    output:
    path 'rnaseq_top5.parquet'

    script:
    """
    ${params.python} ${projectDir}/../scripts/rnaseq_top5.py ${sheet} . ${coding} rnaseq_top5.parquet
    """
}

process ASSEMBLE {
    label 'medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sheet
    path histones
    path top5
    path drug
    path gene_rt

    output:
    path 'notebook_inputs', emit: inputs

    script:
    """
    ${params.python} ${projectDir}/../scripts/assemble_inputs.py ${sheet} ${histones} ${top5} notebook_inputs ${gene_rt}
    mkdir -p notebook_inputs/drug && cp ${drug} notebook_inputs/drug/drug_per_gene.csv
    """
}

process NOTEBOOK {
    label 'medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path inputs

    output:
    path 'figures/*'
    path 'analysis.executed.ipynb'

    script:
    """
    # Register the pipeline's Python as a private kernel, so that no other "python3" kernel is used.
    ${params.python} -m ipykernel install --prefix kernel --name hlb-rt > /dev/null
    export JUPYTER_PATH=\$PWD/kernel/share/jupyter
    RT_DATA=\$PWD/${inputs} RT_FIGDIR=\$PWD/figures ${params.python} -m nbconvert --to notebook --execute \
        --ExecutePreprocessor.kernel_name=hlb-rt \
        ${projectDir}/../analysis.ipynb --output-dir . --output analysis.executed.ipynb
    """
}

// ---- Drug-treatment Repli-seq (Zenodo 10.5281/zenodo.19093559) ----

process EXTRACT_DRUG {
    label 'small'
    storeDir "${params.cache_dir}/drug_fractions"

    input:
    path archive
    path sheet

    output:
    path '*.bw'

    script:
    """
    tail -n +2 ${sheet} | cut -d, -f5 > members.txt
    unzip -q -j ${archive} \$(cat members.txt)
    tail -n +2 ${sheet} | awk -F, '{print \$6"  "\$4}' | md5sum -c --quiet -
    """
}

process WAVELET_DRUG {
    tag "${condition}"
    label 'small'
    publishDir "${params.outdir}/drug/tracks", mode: 'copy'

    input:
    tuple val(condition), val(name), path(fractions, stageAs: 'frac?/*')   // G1b S1 S2 S3 S4 G2, in order

    output:
    tuple val(condition), path("${name}_rt_J7.bw")

    script:
    """
    ${params.python} ${projectDir}/../scripts/wavelet_rt.py ${name}_rt_J7.bw ${fractions.join(' ')}
    """
}

process DRUG_GENE_RT {
    label 'small'
    publishDir "${params.outdir}/drug", mode: 'copy'

    input:
    val specs            // [condition, track] pairs
    path tracks
    path coding
    path histones

    output:
    path 'drug_per_gene.csv'

    script:
    def args = specs.collect { c, t -> "'${c}=${t.name}'" }.join(' ')
    """
    ${params.python} ${projectDir}/../scripts/drug_gene_rt.py ${coding} ${histones} drug_per_gene.csv ${args}
    """
}

// ---- QC: rebuild the UW Repli-seq tracks from their six fractions and compare with ENCODE ----

process REPROCESS_UW {
    tag "${id}"
    label 'small'
    publishDir "${params.outdir}/qc/reprocessed_uw", mode: 'copy'

    input:
    tuple val(id), val(label), path(fractions, stageAs: 'frac?/*')   // G1b S1 S2 S3 S4 G2, in order

    output:
    tuple val(id), val(label), path("${id}_reprocessed.bw")

    script:
    """
    ${params.python} ${projectDir}/../scripts/wavelet_rt.py ${id}_reprocessed.bw ${fractions.join(' ')}
    """
}

process QC_TRACKS {
    label 'medium'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path pairs_csv
    path bigwigs

    output:
    path 'track_qc/*'

    script:
    """
    ${params.python} ${projectDir}/../scripts/qc_tracks.py track_qc ${pairs_csv}
    """
}
