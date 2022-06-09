process EPYTOPE_PEPTIDE_PREDICTION {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::coreutils=9.1 conda-forge::tcsh=6.20.00 bioconda::epytope=3.0.0 conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3774d4f6160bf3a5a53d47875424448104ba6d22:79da0cdeada554621d205f5713988db45f65fb80-0' :
        'quay.io/biocontainers/mulled-v2-3774d4f6160bf3a5a53d47875424448104ba6d22:79da0cdeada554621d205f5713988db45f65fb80-0' }"

    input:
    tuple val(meta), path(splitted), path(software_versions)
    path netmhc_paths

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.tsv"), emit: predicted
    tuple val(meta), path("*.fasta"), emit: fasta optional true
    path "versions.yml", emit: versions

    script:
    // Additions to the argument command need to go to the beginning.
    // Argument list needs to end with --peptides or --somatic_mutation
    def argument = task.ext.args

    if (params.proteome) {
        argument = "--proteome ${params.proteome} " + argument
    }

    if (params.wild_type) {
        argument = "--wild_type " + argument
    }

    if (params.fasta_output) {
        argument = "--fasta_output " + argument
    }

    if (params.tool_thresholds) {
        argument = "--tool_thresholds ${params.tool_thresholds} " + argument
    }

    """
    # create folder for MHCflurry downloads to avoid permission problems when running pipeline with docker profile and mhcflurry selected
    mkdir -p mhcflurry-data
    export MHCFLURRY_DATA_DIR=./mhcflurry-data
    # specify MHCflurry release for which to download models, need to be updated here as well when MHCflurry will be updated
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=1.4.0
    # Add non-free software to the PATH
    shopt -s nullglob
    for p in ${netmhc_paths} ; do export PATH="\$(realpath -s "\$p"):\$PATH"; done
    shopt -u nullglob

    epaa.py --identifier ${splitted.baseName} \
        --alleles '${meta.alleles}' \
        --versions ${software_versions} \
        ${argument} ${splitted}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyvcf').version)")
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}