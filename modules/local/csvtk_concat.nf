process CSVTK_CONCAT {
    label 'process_low'

    conda "bioconda::csvtk=0.23.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(predicted)

    output:
    tuple val(meta), path("*.tsv"), emit: concatenated
    path "versions.yml", emit: versions

    script:
    """
    csvtk concat -t $predicted > ${meta.sample}_predicted.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
