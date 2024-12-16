//and also do the allele name

process PREPARE_PREDICTION_INPUT {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcgnomes=1.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.4--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcgnomes:1.8.4--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.csv|*.tsv"), emit: prepared
    path "versions.yml", emit: versions

    script:
    //TODO handle the thresholds (parse the --tools_thresholds and --use_affinity_thresholds)

    template "prepare_prediction_input.py"

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_syfpeithi.csv
    touch ${prefix}_mhcflurry.csv
    touch ${prefix}_mhcnuggets.csv
    touch ${prefix}_netmhcpan.csv
    touch ${prefix}_netmhciipan.csv
    touch versions.yml
    """
}