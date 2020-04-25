#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/tpp_workflow
========================================================================================
(NOT YET A nf-core!)
 #### Homepage / Documentation
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run main.nf  -profile docker

    For a test run with pre-given data, use:
    nextflow run main.nf -profile docker, test


    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      --fasta                Fasta file for database search
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (da)
      --digest_mass_range               Mass range of peptides considered for matching
      --activation_method               Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                          Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see comet enzymes)
      --miscleavages                    Number of allowed miscleavages
      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      TODO: Not working yet. Fixed modification is always 'Carbamidomethyl (C)'
      --variable_mods                   TODO: Only variable modifications 'Oxidation (M)' and 'Phosphorylation (STY)' allowed. Separated by commas without additional spaces
      --num_hits                        Number of reported hits
     --min_charge                       Minimal precursor charge 
     --max_charge                       Maximal precursor charge 
      --skip_decoy_generation           Use a fasta databse that already includes decoy sequences
      --comet_param_file                Parameter file for comet database search. This will overwrite all other parameters for the comet search
      
      
    Identification validation:
      --peptide_min_length              Minimum accepted peptide length
      --fdr_peptide_threshold           Threshold for FDR filtering
      
      
      
      --quantification_fdr              FDR threshold to accept peptides for quantification
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


// Validate inputs
params.raws = params.raws ?: { log.error "No read data privided. Make sure you have used the '--raws' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
//params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.1
params.precursor_mass_tolerance = 5
params.fragment_bin_offset = 0
params.fdr_peptide_threshold = 0.01
params.number_mods = 3

params.comet_param_file = "none"

params.num_hits = 1
params.pick_ms_levels = 2
params.run_centroidisation = false

params.min_charge = 2
params.max_charge = 3
params.activation_method = 'ALL'

params.enzyme = 'Trypsin'
params.miscleavages = 1
params.fixed_mods = 'Carbamidomethylation of C'
params.variable_mods = 'Oxidation of M'


params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.quantification_fdr = 0.01
params.quantification_min_prob = 0


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
        .fromPath( params.fasta ).into {input_fasta; input_fasta2; input_fasta3}

/*
 * Create a channel for comet parameter file
 */
  input_comet_param =   Channel
        .fromPath( params.comet_param_file )

if (params.comet_param_file == "none") {
} else if( !(file(params.comet_param_file).exists()) ) {
        log.error "Comet parameter file does not exit"; exit 1
    
}



/*
 * STEP 1 - convert raw files to mgf
 */
process convert_raw_mzml {
    publishDir "${params.outdir}"
    input:
      file rawfile from input_raw

    output:
     file "${rawfile.baseName}.mzML" into (mzml, mzml2)

    script:
     """
     ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 2 -z
     """

}

/*
 * STEP 3 - set comet parameter file
 */
process set_comet_configuration {
    publishDir "${params.outdir}"
    input:
      file comet_param from input_comet_param.ifEmpty(file("none"))
     output:
      file "comet.params" into (comet_param_file)
    script: 
    // no file provided
    if (comet_param.getName() == "none") {  
      ttrue = true
      enzymemap = ["Trypsin": 1, "Trypsin/P": 2, "Lys_C": 3, "Lys_N": 4, "Arg_C": 5, "Asp_N": 6, "CNBr": 7, "Glu_C": 8, "PepsinA": 9, "Chymotrypsin": 10, "Unspecified": 0]
      enzyme = enzymemap[params.enzyme]
      if (enzyme == null) enzyme = enzymemap["Unspecified"]
      skip_decoy = params.skip_decoy_generation
      modmap = ["Oxidation of M": "15.9949 M 0 3 -1 0 0 0.0", "Phosphorylation of STY": "79.966331 STY 0 3 -1 0 0 97.976896", "none": "0.0 X 0 3 -1 0"]
      mods = params.variable_mods.split(",")  
      modout = ""
      for (int i=0; i<10; i++) {
        if(mods.size() > i ) {
          tmod =  modmap[mods[i]]
          if (tmod != null) {
              modout += "variable_mod0" + i + " = " + tmod + "\n"
          } else {
              modout += "variable_mod0" + i + " = " + modmap["none"] + "\n"     
          }
        } else {
            modout += "variable_mod0" + i + " = " + modmap["none"] + "\n"     
        }
      }
      """
      comet -p
      mv comet.params.new comet.params
      sed -i '/decoy_search/d' comet.params
      if [ "${skip_decoy}" = "${ttrue}" ]
          then
             echo "decoy_search = 0" >> comet.params
      else
             echo "decoy_search = 1" >> comet.params
      fi
      # given in ppm
      sed -i '/num_threads/d' comet.params
      echo "num_threads = ${task.cpus}" >> comet.params
      sed -i '/peptide_mass_tolerance/d' comet.params
      echo "peptide_mass_tolerance = ${params.precursor_mass_tolerance}" >> comet.params
      sed -i '/search_enzyme_number/d' comet.params
      echo "search_enzyme_number = ${enzyme}" >> comet.params
      # given in da and assuming high resolution
      sed -i '/fragment_bin_tol/d' comet.params
      sed -i '/fragment_bin_offset/d' comet.params
      echo "fragment_bin_tol = ${params.fragment_mass_tolerance}" >> comet.params
      echo "fragment_bin_offset = 0.0" >> comet.params
      # mass range
      sed -i '/activation_method/d' comet.params
      echo "activation_method = ${params.activation_method}" >> comet.params
      sed -i '/allowed_missed_cleavage/d' comet.params
      echo "allowed_missed_cleavage = ${params.miscleavages}" >> comet.params
      sed -i '/max_variable_mods_in_peptide/d' comet.params
      echo "max_variable_mods_in_peptide = ${params.number_mods}" >> comet.params
      sed -i '/num_results/d' comet.params
      echo "num_results = ${params.num_hits}" >> comet.params
      sed -i '/precursor_charge/d' comet.params
      echo "precursor_charge = ${params.min_charge} ${params.max_charge}" >> comet.params  
      sed -i '/variable_mod0/d' comet.params
      echo "${modout}" >> comet.params
      """ 
       // with a given comet file, all other paramter will be discarded
    } else {
          """
          cp ${comet_param} comet.params
          """
    }
}

      
 
/*
 * STEP 3 - run comet
*/ 

process db_search_comet {
    publishDir "${params.outdir}"
    input:
     each mzml_file from mzml
     file fasta from input_fasta
     file comet_param from comet_param_file

    output:
     file "${mzml_file.baseName}.pep.xml" into pepxml

    script:
     """
     comet -P"${comet_param}" -N"${mzml_file.baseName}" -D"${fasta}" "${mzml_file}"
     """
}


/*
 * STEP 4 - run PeptideProphet to validate PSMs
*/ 

process run_peptideprophet {
    publishDir "${params.outdir}"
    input:
     each pepxml_file from pepxml
     file fasta from input_fasta2

    output:
     file "${pepxml_file.baseName}.interact.pep.xml" into interact_pepxml
     
    script:
     enzymemap = ["Trypsin": "", "Trypsin/P": "", "Lys_C": "-eN", "Lys_N": "-eL", "Arg_C": "-eN", "Asp_N": "-eA", "CNBr": "-eM", "Glu_C": "-eG", "PepsinA": "-eN", "Chymotrypsin": "-eC", "Unspecified": "-eN"]
     enzyme = enzymemap[params.enzyme]

    """
    xinteract -N"${pepxml_file.baseName}.interact.pep.xml" -p"${params.fdr_peptide_threshold}" ${enzyme} -l"${params.peptide_min_length}" -THREADS=${task.cpus} -PPM -O -D"${fasta}" "${pepxml_file}"
    """
}

/*
 * STEP 4 - run ProteinProphet to validate on protein level
*/ 

process run_proteinprophet {
    publishDir "${params.outdir}"
    input:
     file pepxml_file from interact_pepxml

    output:
     file "${pepxml_file.baseName}.prot.xml" into (protxml, protxml2)
     
    script:
    """
    ProteinProphet "${pepxml_file}" "${pepxml_file.baseName}.prot.xml"
    """
}

/*
 * STEP 4 - run StPepter for protein quantification (label-free)
*/ 

process run_stpeter {
    publishDir "${params.outdir}"
    input:
     each file(protxml) from protxml 
     file fasta from input_fasta3

    output:
     file "*csv" into protquant
     
    script:
    """
    StPeter -f ${params.quantification_fdr}  "${protxml}"
    """
    
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

   
