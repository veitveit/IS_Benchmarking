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
        path "${params.loc_out}" into allinput
        //path "absence_peak_data/output/**" into ALLOUTPUT
        //file './**' into ALLOUTPUT
        //file "${rawfile.baseName}.mzML" into (mzml, mzml2)

    script:
        // docker run -it veitveit/moffworkflow:dev bash
        // docker run -v data:/data -v results:/results  -i -t veitveit/moffworkflow:dev python /moFF/moff_all.py --config_file /moFF/configuration_file.ini -i "\$(find "/data/raw/" -type f -name "*.raw" | head -n1;)" --tol 10 -rt_win_peak 1 --xic_length 3 --loc_out /reuslts/ --mbr off

        """
        #find / -name "moff_all.py" -type f 2> /dev/null;
        #find / -name "${inifile}" -type f 2> /dev/null;

        repodir="/opt/conda/envs/nf-core-moff-ms1/share/moff-2.0.3-3/";


        echo -e "import pandas\npandas.show_versions()" | python3;

        #echo "${inifile}";

        #python3 "\${repodir}/moff_all.py" --loc_out "${params.outdir}" --config_file "${inifile}" --loc_in "${rawrepodir}" --ext "txt"; # --raw_repo "${rawrepodir}"
        #python3 "\${repodir}/moff_mbr.py" --loc_out "${params.outdir}" --loc_in "${rawrepodir}" --mbr only;
        python3 "\${repodir}/moff_all.py" --loc_out "${params.outdir}"  --config_file "${inifile}";

        #find "./" -type f 2> /dev/null;
        """
}

stdout_channel.view { "stdout:\t$it" }


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

   
