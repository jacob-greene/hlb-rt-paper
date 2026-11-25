// Download one public file and check its md5 against the sample sheet.
process FETCH {
    tag "${id}"
    label 'small'
    storeDir "${params.cache_dir}/${kind}"
    maxForks params.max_downloads
    errorStrategy { task.attempt < 3 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(id), val(kind), val(url), val(md5), val(fname)

    output:
    tuple val(id), path(fname)

    script:
    """
    curl -sSfL --retry 5 -o ${fname} '${url}'
    echo "${md5}  ${fname}" | md5sum -c --quiet -
    """
}
