process NETMHCIIPAN {
    label 'process_single'
    tag "${metadata.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(metadata), path(peptide_file), path(nonfree_tools)

    output:
    tuple val(metadata), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (metadata.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCPAN for MHC class I, or adjust the samplesheet accordingly."
    }
    def prefix = peptide_file.baseName

    """
    netmhciipan.py \
        --input ${peptide_file} \
        --output ${prefix}_predicted_netmhciipan.tsv \
        --alleles '${metadata.alleles}' \
        --sample_id ${metadata.sample}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python \$(python --version | sed 's/Python //g')

    END_VERSIONS
    """
    stub:
    """
    touch ${metadata.sample}_predicted_netmhciipan.tsv
    touch versions.yml
    """
}
