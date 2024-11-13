process PARSE_MMSEQS_PROTEIN_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9 conda-forge::pandas=1.3.0 anaconda::seaborn=0.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' :
        'biocontainers/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0' }"

    input:
    tuple val(meta), path(table_lca)

    output:
    tuple val(meta), path('*_selected.tsv'), emit: table
    tuple val(meta), path('*_selected.txt')                 , emit: as_list
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_mmseqs_proteins_taxonomy.py \\
        -in_lca $table_lca \\
        -out_p $prefix \\
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
    END_VERSIONS    
    """

    stub:
    def out1 = "${prefix}_selected.tsv"
    def out2 = "${prefix}_selected.txt"
    """
    touch $out1
    touch $out2
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
    END_VERSIONS    
    """
}
