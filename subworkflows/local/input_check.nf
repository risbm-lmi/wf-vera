//
// Check input samplesheet and get read channels
//

// #TODO: устроить проверку всех обрабатываемых исключений, пока что взято as is

// just took from random nf-core pipeline
// checks wether the input extension is correct 
// вообще можно и в скриптик обернуть питоновский, так тоже немало кто делает!
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow INPUT_CHECK {
    main:
        if(hasExtension(params.input, ".csv")){
            ch_input_rows = Channel
                .from(file(params.input))
                .splitCsv(header: true)
                .map {
                    row -> 
                        if (row.size() == 4) {
                            def id = row.sample
                            def run = row.run
                            def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false
                            def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : false 
                            // Check if given combination is valid
                            if (run != null && run == "") exit 1, "ERROR: Check input samplesheet -> Column 'run' contains an empty field."
                            if (!sr1 || !sr2) exit 1, "Invalid input PE samplesheet: short reads can not be empty."
                            return [ id, run, sr1, sr2 ]
                        } else if (row.size() == 3) {
                            def id = row.sample
                            def run = row.run
                            def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false 
                            // Check if given combination is valid
                            if (run != null && run == "") exit 1, "ERROR: Check input samplesheet -> Column 'run' contains an empty field."
                            if (!sr1 ) exit 1, "Invalid input SE samplesheet: short reads can not be empty."
                            return [ id, run, sr1 ]
                        } else {
                            exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 3 or 4."
                        }
                }
        } else {
            exit 1, "An input CSV file with sample pathes and description is expected"
        }

        ch_raw_reads = ch_input_rows.map { it -> 
                                            if (it.size() == 3) {
                                                def meta = [:]
                                                    meta.id = it[0]
                                                    meta.run = it[1]
                                                    meta.single_end = true
                                                return [meta, it[2]]
                                        } else {
                                                def meta = [:]
                                                    meta.id = it[0]
                                                    meta.run = it[1]
                                                return [ meta, [ it[2], it[3] ]]
                                        }
        }
        
        // Ensure run IDs are unique within samples (also prevents duplicated sample names)
        // пока не могу разобраться, почему тут не используется exit, но используется error :O 
        ch_input_rows
            .groupTuple(by: 0)
            .map { it -> if( it[1].size() != it[1].unique().size() ) { { error("ERROR: input samplesheet contains duplicated sample or run IDs (within a sample)! Check samplesheet for sample id: ${it[1]}") } } }
    
        // тут добавить проверку ncbi_api_token и иных ключей, а затем обернуть это в 

    emit:
        raw_short_reads = ch_raw_reads
}
