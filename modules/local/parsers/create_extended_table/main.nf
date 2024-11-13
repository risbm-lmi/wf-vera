process CREATE_EXTENDED_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.3.0 anaconda::seaborn=0.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' :
        'biocontainers/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' }"

    input:
    tuple val(meta), path(in_nucl), path(in_prot), path(in_cove)

    output:
    tuple val(meta), path('*_extended_table.tsv'), emit: table
    path('*_for_extension.txt')                 , emit: as_list
    path('*_two_col_coverage.tsv')                 , emit: two_col_coverage
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_extended_table.py \\
        -in_nucl $in_nucl \\
        -in_prot $in_prot \\
        -in_cove $in_cove \\
        -out_p $prefix \\
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
    END_VERSIONS    
    """

    stub:
    def out1 = "${prefix}_extended_table.tsv"
    def out2 = "${prefix}_for_extension.txt"
    def out3 = "${prefix}_two_col_coverage.tsv"
    """
    touch $out1
    touch $out2
    touch $out3
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
    END_VERSIONS    
    """
}
