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
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
// params.absence = params.absence ?: { log.error "No absence peak data provided. Make sure you have used the '--absence' option."; exit 1 }()
params.inifile = params.inifile ?: { log.error "No absence peak data provided. Make sure you have used the '--absence' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
//params.loc_in = params.loc_in ?: { log.error "No output directory provided. Will put the results into './results'"; return "./results" }()
params.raw_repo = params.raw_repo ?: { log.error "21"; exit 1 }()


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

// Channel.fromPath( params.absence ).set { input_absence; }

Channel.fromPath( params.inifile ).set { input_ini; }


Channel.fromPath( params.raw_repo ).into { raw_repo; loc_in; }



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
 * STEP 1 - convert raw files to mgf
 */
process moff_all {
    echo true
    publishDir "${params.outdir}"

    container "veitveit/moffworkflow:dev"


    input:
        //file rawfile from input_raw
        file inifile from input_ini
        file rawrepodir from raw_repo
        //file loc_in from loc_in
        //file absencefile from input_absence


    output:
        stdout stdout_channel
        path "absence_peak_data/raw_repo/apex_output" into apex_output
        path "absence_peak_data/mbr_output/" into mbr_output

    script:
        """
        source activate moFF;

        moff_filepath=\$(which moff_all.py;);

        #python3 "\${moff_filepath}" --loc_out "${params.outdir}" --config_file "${inifile}" --loc_in "${rawrepodir}" --ext "txt"; # --raw_repo "${rawrepodir}"
        #python3 "\${moff_filepath}" --loc_out "${params.outdir}" --loc_in "${rawrepodir}" --mbr only;
        #python3.6 "\${moff_filepath}" --loc_out "${params.outdir}"  --config_file "${inifile}";
        python3.6 "\${moff_filepath}" --config_file "${inifile}" --peptide_summary 2>&1;
        """
}

stdout_channel.view { "stdout:\t$it" }


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The results should now be in the same folder as the input data" : "Oops .. something went wrong" )
}

   
