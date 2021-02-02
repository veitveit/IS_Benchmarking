#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/MaxQuant
========================================================================================
(NOT YET A nf-core!)
 #### Homepage / Documentation
----------------------------------------------------------------------------------------
*/

def helpMessage() {
log.info nfcoreHeader()
log.info"""
Usage: 

The typical command for running the pipeline is as follwos:
nextflow run main.nf --raws "path/to/raws" --sdrf "path/to/*.tsv" --fasta "path/to/*.fasta" -profile docker

Mandatory arguments: 
--raws      Path to input data (must be surrounded with quotes)
--sdrf      Path to sdrf.tsv file (must be surrounded with quotes, and in tsv format)
--fasta     Path tp the fasta file(must be surrounded with quotes, and in .fasta)



-profile    Configuration profile to use. Can use multiple (comma separated)
            Available: docker


Other options:
    --numthreads    each thread needs at least 2 GB of RAM,number of threads should be ≤ number of logical cores available          
    --peptidefdr    posterior error probability calculation based on target-decoy search, default=0.01 
    --proteinfdr    protein score = product of peptide PEPs (one for each sequence)', default=0.01
    --match         "True" via matching between runs to boosts number of identifications, defualt = False
    --tempfol       Path to a temp folder, if missing a default will be created in the basedirectory.
    --outdir        The output directory where the results will be saved
    --outPara       The output name for the xml for openms. (should be named './<file>.xml')
    --outExpDes     The output for the experimental design. (Should be name './<filename>.txt')
    --email         Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
    -name           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

""".stripIndent()
}

/*
* Read and setup the variables
*/ 
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.sdrf = params.sdrf ?: {log.error "No read data provided. Make sure you have used the '--sdrf' option."; exit 1}()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1}()
params.match = params.match ?: {log.warn "No match choice. Make sure to choose and write '--match'"; return "False"}()
params.tempfol = params.tempfol ?: { log.warn "No tmp directory provided. Will put the tmp into './tmp'"; return "./tmp" }()
 
/*
* Define the default paramaters. 
* The parameters are found in the xml file and is change within the file.
*/

// SDRF parameters
params.numthreads = 8
params.peptidefdr = 0.01
params.proteinfdr = 0.01
params.outPara ='./parameters.xml'
params.outExpDes = './experimentalDesign.txt'

params.proteinGroups = "./proteinGroups.txt"
params.mqResults = "MaxQuantResults"

// NormalyzerDE parameters
params.project = "PXD001819"
// Which groups to compare
params.comparisons = 'c("1-2","2-3")'
params.normalyzerDesign = "${params.outdir}/normalyzer_design.txt"
params.normalyzerMethod = "log2"

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
* Generate the channels for the raw files
*/ 
Channel 
    .fromPath (params.raws).into {input_raw; input_raw2}

/*
* Generate the channels for the sdrf files
*/ 
Channel 
    .fromPath (params.sdrf).into {input_sdrf; input_sdrf2}

/*
* Generate the channels for the fasta file
*/ 
Channel
    .fromPath (params.fasta).into {input_fasta; input_fasta2}



/* 
* STEP 1 - Generate the parameter and experimental design for maxquant through sdrf
*/
process run_sdrf {
    publishDir "${params.outdir}"
    input: 
        path sdrf_file from input_sdrf
        path fasta_file from input_fasta
    
    output:
        file "${params.outPara}" into outParameters
        file "${params.outExpDes}" into outExpDesign

    script: 
    """
    parse_sdrf convert-maxquant -s ${params.sdrf} -f ${params.fasta} -m ${params.match} -pef ${params.peptidefdr} -prf ${params.proteinfdr} -t ${params.tempfol} -r ${params.raws} -n ${params.numthreads} -o1 ${params.outPara} -o2 ${params.outExpDes}
    """
}

/*
* STEP 2 - Run maxQuant with parameters from mqpar
*/
process run_maxquant {
   
    publishDir "${params.outdir}"

    input:
        file rawfile from input_raw
        path mqparameters from outParameters

    output:
        file "${params.proteinGroups}" into input_proteinGroups	


    script:
    """
    maxquant ${mqparameters}
    mkdir -p "${params.outdir}/${params.mqResults}"
    cp -R "${params.raws}"/combined/txt/* "${params.outdir}/${params.mqResults}"
    cp "${params.raws}/combined/txt/proteinGroups.txt" "${params.proteinGroups}"
    """        
}

/*
* STEP 3 - Run NormalyzerDE
*/
process run_normalyzerde {

    publishDir "${params.outdir}"

    input:
        path sdrf_file from input_sdrf2
        path exp_file from outExpDesign
        path protein_file from input_proteinGroups
    script:
    """
#    paste ${exp_file} ${sdrf_file} | awk '{OFS="\\t"}{print \$3, \$NF}' > ${params.normalyzerDesign}
#    sed -i '1s/.*/sample\\tgroup/' ${params.normalyzerDesign}
    Rscript -e 'NormalyzerDE::normalyzer(jobName="${params.project}", designPath="${params.normalyzerDesign}", dataPath="${protein_file}", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="${params.outdir}")'
    Rscript -e 'NormalyzerDE::normalyzerDE(jobName="${params.project}", comparisons=${params.comparisons}, designPath="${params.normalyzerDesign}", dataPath="${params.outdir}/${params.project}/${params.normalyzerMethod}-normalized.txt", outputDir="${params.outdir}")'


    """   
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}
