#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/moff_workflow
========================================================================================
 nf-core/moff_workflow Analysis Pipeline.
 #### Homepage / Documentation
TODO https://github.com/nf-core/moff_workflow
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
        This tool has no proper documentation yet!
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Validate inputs
params.txts = params.txts ?: { log.error "No read data provided. Make sure you have used the '--txts' option."; exit 1 }()
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.inifile = params.inifile ?: { log.error "No ini configuration file provided. Make sure you have used the '--inifile' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


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
 * Create channels for input files
 */
Channel.fromPath( params.raws ).set { input_raw }
Channel.fromPath( params.txts ).set { input_txt }
Channel.fromPath( params.inifile ).set { input_ini; }


/*
 * STEP 0 - 100% dysfunctional for now
 */
// process get_peptideshaker_tsv {
//     publishDir "${params.outdir}"
//     input:
//         tuple file(pepshaker), file(mgffile) from peptideshaker_file

//     output:
//       file "${params.name}_${mgffile.baseName}_1_Default_PSM_Report_with_non-validated_matches.txt" into (peptideshaker_tsv_file)

//     script:
//         """
//         peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
//         peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in "./${pepshaker}" -out_reports "./" -reports "4"
//         """    
// }


/*
 * STEP 1 - Run moFF
 */
process moff_all {
    echo true
    publishDir "${params.outdir}"

    container "veitveit/moffworkflow:dev"

    input:
        file txt_filepath from input_txt
        file raw_filepath from input_raw
        file inifile from input_ini

    output:
        stdout stdout_channel

    script:
        """
        source activate moFF;

        moff_filepath=\$(which moff_all.py;);

        python3.6 \${moff_filepath} --config_file "${inifile}" --loc_in "${txt_filepath}" --raw_repo "${raw_filepath}" --peptide_summary 2>&1
        """
}

stdout_channel.view { "stdout:\t$it" }


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The results should now be in the same folder as the input data" : "Oops .. something went wrong" )
}

   
