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
params.absence = params.absence ?: { log.error "No absence peak data provided. Make sure you have used the '--absence' option."; exit 1 }()
params.inifile = params.inifile ?: { log.error "No absence peak data provided. Make sure you have used the '--absence' option."; exit 1 }()
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

Channel.fromPath( params.absence ).set { input_absence; }

Channel.fromPath( params.inifile ).set { input_ini; }


/*
 * STEP 1 - convert raw files to mgf
 */
process moff_all {
    echo true
    publishDir "${params.outdir}"

    container "veitveit/moffworkflow:dev"


    input:
        file rawfile from input_raw
        file inifile from input_ini
        //file absencefile from input_absence


    output:
        file 'docker.log' into ( output )
        //file "${rawfile.baseName}.mzML" into (mzml, mzml2)

    script:
        // docker run -it veitveit/moffworkflow:dev bash
        // docker run -v data:/data -v results:/results  -i -t veitveit/moffworkflow:dev python /moFF/moff_all.py --config_file /moFF/configuration_file.ini -i "\$(find "/data/raw/" -type f -name "*.raw" | head -n1;)" --tol 10 -rt_win_peak 1 --xic_length 3 --loc_out /reuslts/ --mbr off

        """
        python3 /opt/conda/envs/nf-core-moff-ms1/share/moff-2.0.3-3/moff_all.py --config_file ${inifile}  -i ${rawfile} -o ./  -f 2 -z 1> docker.log
        """
}



workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

   
