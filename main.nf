nextflow.enable.dsl=2

log.info """\
    B E T T Y Core - N F   P I P E L I N E
    ======================================
    table        : ${params.input}
    outdir       : ${params.outdir}
    """
    .stripIndent()

include { VERA } from './workflows/vera'

workflow {
    VERA ()
}

