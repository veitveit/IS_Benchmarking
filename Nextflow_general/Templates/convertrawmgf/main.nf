#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/prottemplate
========================================================================================
 nf-core/prottemplate Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/prottemplate
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --mzmls '*.mzML' --fasta '*.fasta' -profile standard,docker

    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --peptide_min_length              Minimum peptide length for filtering
      --peptide_max_length              Maximum peptide length for filtering
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (ppm)
      --fragment_bin_offset             Offset of fragment mass bin (Comet specific parameter)
      --use_x_ions                      Use x ions for spectral matching in addition
      --use_z_ions                      Use z ions for spectral matching in addition
      --use_a_ions                      Use a ions for spectral matching in addition
      --use_c_ions                      Use c ions for spectral matching in addition
      --fdr_threshold                   Threshold for FDR filtering
      --fdr_level                       Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --digest_mass_range               Mass range of peptides considered for matching
      --activation_method               Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                          Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods                   Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --num_hits                        Number of reported hits
      --run_centroidisation             Specify whether mzml data is peak picked or not (true, false)
      --pick_ms_levels                  The ms level used for peak picking (eg. 1, 2)
      --prec_charge                     Precursor charge (eg. "2:3")
      --max_rt_alignment_shift          Maximal retention time shift (sec) resulting from linear alignment      
      --spectrum_batch_size             Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)
      --description_correct_features    Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration) 
      --klammer                         Retention time features are calculated as in Klammer et al. instead of with Elude.
      --skip_decoy_generation           Use a fasta databse that already includes decoy sequences
      --quantification_fdr              Assess and assign ids matched between runs with an additional quantification FDR
      --quantification_min_prob         Specify a minimum probability cut off for quantification

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

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


// Validate inputs
params.raws = params.raws ?: { log.error "No read data privided. Make sure you have used the '--raws' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.02
params.precursor_mass_tolerance = 5
params.use_x_ions = false
x_ions = params.use_x_ions ? '-use_X_ions true' : ''
params.use_z_ions = false
z_ions = params.use_z_ions ? '-use_Z_ions true' : ''
params.use_a_ions = false
a_ions = params.use_a_ions ? '-use_A_ions true' : ''
params.use_c_ions = false
c_ions = params.use_c_ions ? '-use_C_ions true' : ''
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

params.prec_charge = '2:3'
params.activation_method = 'ALL'

params.enzyme = 'unspecific cleavage'
params.fixed_mods = ''
params.variable_mods = 'Oxidation (M)'
params.spectrum_batch_size = 500

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
 * STEP 1 - convert raw files to mzML
 */
process convert_raw_mzml {
        input:
	  file rawfile from input_raw


	output:
	 file "${rawfile.baseName}.mzML" into (mzmls)

	script:
	 """
	 ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 1
	 """

}
