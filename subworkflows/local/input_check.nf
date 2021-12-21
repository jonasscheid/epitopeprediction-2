//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {

    def allele_string = generate_allele_string(row.alleles)
    def type = determine_input_type(row.filename)
    def meta = [:]
    meta.sample         = row.sample
    meta.alleles        = allele_string
    meta.inputtype      = type

    def array = []
    if (!file(row.filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> file does not exist!\n${row.Filename}"
    } else {
        array = [ meta, file(row.filename) ]
    }
    return array
}

def generate_allele_string(String alleles) {
    // Collect the allele information from the file
    def allele_string
    if ( alleles.endsWith(".txt") || alleles.endsWith(".alleles") )  {
        allele_string = file(alleles).readLines().join(';')
    }
    // or assign the information to a new variable
    else {
        allele_string = alleles
    }
    return allele_string
}

def determine_input_type(String filename) {
    def filetype
    def input_file = file(filename)
    def extension = input_file.extension

    if ( extension == "vcf" | extension == "vcf.gz" ) {
        filetype = "variant"
    }
    else if ( extension == "tsv" | extension == "GSvar" ) {
        // Check if it is a variant annotation file or a peptide file
        input_file.withReader {
            def first_header_col = it.readLine().split('\t')[0]
            if (first_header_col == "id") { filetype = "peptide" }
            else if (first_header_col == "#chr") {filetype = "variant"}
        }
    }
    else {
        filetype = "protein"
    }
    return filetype
}
