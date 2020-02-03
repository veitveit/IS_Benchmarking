#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/prottemplate
========================================================================================
 nf-core/proline_workflow Analysis Pipeline.
 #### Homepage / Documentation
TODO https://github.com/nf-core/prottemplate
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    You might need to build the docker image. A Dockerfile is available. Name it "veitveit/searchguiworkflow:dev" 
    or change the name in the nextflow.conf

    The typical command for running the pipeline is as follows:
    nextflow run main.nf  -profile docker

    For a test run with pre-given data, use:
    nextflow run main.nf --raws 'folder/file.raw' --fasta 'folder/file.fasta' -profile docker, test


    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      --fasta                Fasta file for database search
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --peptide_min_length              Minimum peptide length for filtering
      --peptide_max_length              Maximum peptide length for filtering
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (ppm)
      --fragment_bin_offset             Offset of fragment mass bin (Comet specific parameter)
      --fions                          Forward ions for spectral matching
      --rions                          Reverse ions for spectral matching
      --fdr_threshold                   Threshold for FDR filtering
      --fdr_level                       Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --digest_mass_range               Mass range of peptides considered for matching
      --activation_method               Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                          Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --miscleavages            Number of allowed miscleavages
      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods                   Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --num_hits                        Number of reported hits
      --run_centroidisation             Specify whether mzml data is peak picked or not (true, false)
      --pick_ms_levels                  The ms level used for peak picking (eg. 1, 2)
     --min_charge                       Minimal precursor charge 
     --max_charge                       Maximal precursor charge 
      --max_rt_alignment_shift          Maximal retention time shift (sec) resulting from linear alignment      
      --spectrum_batch_size             Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)
      --description_correct_features    Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration) 
      --klammer                         Retention time features are calculated as in Klammer et al. instead of with Elude.
      --skip_decoy_generation           Use a fasta databse that already includes decoy sequences
      --quantification_fdr              Assess and assign ids matched between runs with an additional quantification FDR
      --quantification_min_prob         Specify a minimum probability cut off for quantification
      --run_xtandem                     SearchGui runs xtandem database search
      --run_msgf                        SearchGui runs msgf+ database search
      --run_comet                       SearchGui runs comet database search
      --run_ms_amanda                   SearchGui runs msamanda database search
      --run_myrimatch                   SearchGui runs myrimatch database search
      
    Options for Proline:

    Other options:
      --outdir                          The output directory where the results will be saved
      --email                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                        The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                       The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Validate inputs
params.raws = params.raws ?: { log.error "No read data privided. Make sure you have used the '--raws' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.02
params.precursor_mass_tolerance = 5
params.fions = "b"
params.rions = "y"
params.fragment_bin_offset = 0
params.fdr_threshold = 0.01
params.fdr_level = 'peptide-level-fdrs'
fdr_level = (params.fdr_level == 'psm-level-fdrs') ? '' : '-'+params.fdr_level
params.description_correct_features = 0
params.klammer = false
params.number_mods = 3

params.num_hits = 1
params.digest_mass_range = "800:2500"
params.pick_ms_levels = 2
params.run_centroidisation = false

params.min_charge = 2
params.max_charge = 3
params.activation_method = 'ALL'

params.enzyme = 'Trypsin'
params.miscleavages = 1
params.fixed_mods = 'Carbamidomethylation of C'
params.variable_mods = 'Oxidation of M'
params.spectrum_batch_size = 500
params.run_xtandem = 1
params.run_msgf = 0
params.run_comet = 0
params.run_ms_amanda = 0
params.run_myrimatch = 0
if (params.run_xtandem == 0  && params.run_msgf == 0 && params.run_comet == 0 && params.run_ms_amanda == 0 && params.run_myrimatch == 0) {
           log.error "No database engine defined. Make sure you have set one of the --run_searchengine options to 1 (searchengine can be xtandem, msgf, comet, ms_amanda, myrimatch)."; exit 1 
}

params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.quantification_fdr = false
params.quantification_min_prob = 0
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}


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
 * Create a channel for input raw files
 */
  input_raw =   Channel
        .fromPath( params.raws )

/*
 * Create a channel for fasta file
 */
  Channel
        .fromPath( params.fasta ).into {input_fasta; input_fasta2}



/*
 * STEP 1 - convert raw files to mgf
 */
process convert_raw_mgf {
    publishDir "${params.outdir}"
        input:
      file rawfile from input_raw


    output:
     file "${rawfile.baseName}.mgf" into (mgfs, mgfs2, mgfs3, mgfs4)

    script:
     """
     ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 0 -m 0
     """

}



/*
 * STEP 2 - create  decoy database
 */
process create_decoy_database {
    publishDir "${params.outdir}"
    input:
      file fasta from input_fasta

    output:
      file "${fasta.baseName}_concatenated_target_decoy.fasta" into (fasta_with_decoy)

    when:
      !params.skip_decoy_generation

    script:
     """
         searchgui eu.isas.searchgui.cmd.FastaCLI -in ${fasta} -decoy
         """    
}    



/*
 * STEP 3 - create  searchgui parameter file
 */
process create_searchgui_paramfile {
    publishDir "${params.outdir}"
    input:
      file fasta_decoy from fasta_with_decoy.ifEmpty(input_fasta2)
      

    output:
      file "searchgui.par" into (searchgui_param, searchgui_param2)

    script:
     """
         searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -prec_tol ${params.precursor_mass_tolerance} \\
         -frag_tol ${params.fragment_mass_tolerance} -enzyme ${params.enzyme} -mc ${params.miscleavages}  \\
             -fixed_mods "${params.fixed_mods}" -variable_mods "${params.variable_mods}" -min_charge ${params.min_charge} -max_charge ${params.max_charge} \\
         -fi ${params.fions} -ri ${params.rions} -import_peptide_length_min ${params.peptide_min_length} \\
         -import_peptide_length_max ${params.peptide_max_length} \\
          -db ${fasta_decoy}  -out searchgui.par
         """    
} 



/*
 * STEP 4 - run database search
 */
process run_searchgui_search{
    publishDir "${params.outdir}"
    input:
      each file(mgffile) from mgfs
      file paramfile from searchgui_param
      

    output:
     tuple file("${mgffile.baseName}.zip"), file(mgffile) into (searchgui_out)

    script:
     """
        mkdir tmp
        mkdir log       
        searchgui eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ./tmp -log ./log
         searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files ./  -output_folder ./  -id_params ${paramfile} -threads ${task.cpus} \\
         -xtandem ${params.run_xtandem} -msgf ${params.run_msgf} -comet ${params.run_comet} -ms_amanda ${params.run_ms_amanda} -myrimatch ${params.run_myrimatch}
         mv searchgui_out.zip ${mgffile.baseName}.zip

         """    
}

/*
 * STEP 5 - fdr, ... by PeptideShaker
 */
process run_peptideshaker {
    publishDir "${params.outdir}"
    input:
      tuple file(search_out), file(mgffile) from searchgui_out
      

    output:
      tuple file("${mgffile.baseName}.cpsx"), file(mgffile) into (peptideshaker_file, peptideshaker_file2)

    script:
     """
        mkdir tmp
        mkdir log    
        unzip ${search_out} searchgui.par
        peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files "./${mgffile}"  -identification_files "./${search_out}"  -id_params searchgui.par \\
        -experiment "${params.name}" -sample "${mgffile.baseName}" -out "./${mgffile.baseName}.cpsx" -replicate 1 -threads ${task.cpus}
         """    
}

/*
 * STEP 6 - get PSM report from peptideshaker results (tsv)
 */
process get_peptideshaker_tsv {
    publishDir "${params.outdir}"
    input:
        tuple file(pepshaker), file(mgffile) from peptideshaker_file

    output:
      file "${params.name}_${mgffile.baseName}_1_Default_PSM_Report_with_non-validated_matches.txt" into (peptideshaker_tsv_file)

    script:
     """
        peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in "./${pepshaker}" -out_reports "./" -reports "4"
         """    
}


/*
 * STEP 7 - get PSM report from peptideshaker results (mzid)
*/

process get_peptideshaker_mzid {
    publishDir "${params.outdir}"
    input:
        tuple file(pepshaker), file(mgffile) from peptideshaker_file2

    output:
      file "${pepshaker.baseName}.mzid" into (peptideshaker_mzid_file)

    script:
     """
    peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.MzidCLI -in "./${pepshaker}" -output_file "./${pepshaker.baseName}.mzid" \\
           -contact_first_name Anonyomous -contact_last_name Nextflow -contact_email 'whoami@gmail.com' -contact_address "Greenland" \\
           -organization_name Illuminati -organization_email 'xyz@hell.xyz' -organization_address "California"
         """    
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

   
