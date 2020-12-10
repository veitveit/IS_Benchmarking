#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/proline
========================================================================================
 nf-core/proline_workflow Analysis Pipeline.
#### Homepage / Documentation
TODO https://github.com/nf-core/proline
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run main.nf --raws '*.raw' --fasta '*.fasta' --experiment_design 'test.txt' --lfq_param 'lfq_param_file.txt' -profile docker

    For testing purposes:
    nextflow run main.nf  -profile docker,test 
    (you need to put the file "UPS1_500amol_R1.raw" into the data folder before)
    Link: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/UPS1_500amol_R1.raw
    

    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      --fasta                           Fasta file for database search
      --lfq_param                       Parameter file for Proline
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --peptide_min_length              Minimum peptide length for filtering
      --peptide_max_length              Maximum peptide length for filtering
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (Da)
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
      --skip_decoy_generation           Use a fasta databse that already includes decoy sequences
      --quantification_fdr              Assess and assign ids matched between runs with an additional quantification FDR
      --quantification_min_prob         Specify a minimum probability cut off for quantification
      --run_xtandem                     SearchGui runs xtandem database search
      
    Options for Proline:
      --experiment_design                text-file containing 2 columns: first with mzDB file names and second with names for experimental conditions

    Run statistics
      --run_statistics			Either true or false

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
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1 }()
params.lfq_param = params.lfq_param ?: { log.error "No lfq parameter file for Proline provided. Make sure you have used the '--lfq_param' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.5
params.precursor_mass_tolerance = 30
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
if (params.run_xtandem == 0) {
           log.error "No database engine defined. Make sure you have set one of the --run_searchengine options to 1 (searchengine can be xtandem)."; exit 1 
}

params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.experiment_design = "none"

params.quantification_fdr = false
params.quantification_min_prob = 0
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}


params.run_statistics = true

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
Channel
        .fromPath( params.raws ).set {input_raw}

/*
 * Create a channel for fasta file
 */
  Channel
        .fromPath( params.fasta ).into {input_fasta; input_fasta2}

/* 
 * Create a channel for proline experimental design file
 */
input_exp_design =  Channel.fromPath(params.experiment_design)
if (params.experiment_design == "none") {
    log.warn "No experimental design! All raw files will be considered being from the one and the same experimental condition."
} else if(!(file(params.experiment_design).exists())) {
        log.error "File with experimental design does not exit"; exit 1
    
}

/*
 * Create a channel for Proline LFQ parameter file
 */
  input_proline_param =  Channel
        .fromPath( params.lfq_param)

/*
 * STEP 1.1 - Convert Thermo RAW files to mzDB
 */
process convert_raw_to_mzdb {
    publishDir "${params.outdir}"
    input:
      file rawfile from input_raw

    output:
      file("${rawfile.baseName}.mzDB") into (mzdbs, mzdbs2, mzdbs3)

    script:
     """
     ls -al ./
     thermo2mzdb "${rawfile}"
     mv "${rawfile}.mzDB" "${rawfile.baseName}.mzDB" 
     """
}

/*
 * STEP 1.2 - Convert mzDB files to MGF
 */
process convert_mzdb_to_mgf {
    publishDir "${params.outdir}"
    input:
      file mzdbfile from mzdbs
    
    output:
     file "${mzdbfile.baseName}.mgf" into (mgfs, mgfs2, mgfs3, mgfs4)

    script:
     """
     mzdb2mgf "${mzdbfile}"
     mv "${mzdbfile}.mgf" "${mzdbfile.baseName}.mgf" 
     """
}

/*
 * STEP 1.7 - convert raw files to mzDB
 */
/*process convert_mzml_mzdb {
    publishDir "${params.outdir}"
    input:
      file mzmlfile from mzmls

    output:
     file("${mzmlfile.baseName}.mzDB") into (mzdbs, mzdbs2)

    script:
     """
     /mzdbconvert/mzML2mzDB "${mzmlfile}"
     mv "${mzmlfile}.mzDB" "${mzmlfile.baseName}.mzDB" 
     """

}*/

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
         -import_peptide_length_max ${params.peptide_max_length} -xtandem_quick_acetyl 0 -xtandem_quick_pyro 0 \\
          -db ${fasta_decoy}  -out searchgui.par
         """    
}


/*
 * STEP 4 - run database search
 */
process run_searchgui_search {
    publishDir "${params.outdir}"
    input:
      each file(mgffile) from mgfs
      file paramfile from searchgui_param

    output:
      file "./searchgui_results/${mgffile.baseName}.t.xml" into (search_engine_result_files, search_engine_result_files2)
      //tuple file("${mgffile.baseName}.zip"), file(mgffile) into (searchgui_out)

    script:
     """
        mkdir tmp
        mkdir log

        searchgui eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ./tmp -log ./log
        searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files ./  -output_folder ./  -id_params ${paramfile} -threads ${task.cpus} \\
          -xtandem ${params.run_xtandem} -msgf ${params.run_msgf} -comet ${params.run_comet} -ms_amanda ${params.run_ms_amanda} -myrimatch ${params.run_myrimatch}

        mkdir ./searchgui_results
        unzip searchgui_out.zip -d ./searchgui_results/
        ls -al ./searchgui_results
        """
}

/*
 * STEP 8 - set up Proline config filelist
*/
process set_proline_inputfiles {
    publishDir "${params.outdir}"
    input:
      val result_files from search_engine_result_files.collect()
    
    output:
      file "import_file_list.txt" into import_files

    script:
      def cmd = 'touch import_file_list.txt\n'
      for (int i=0; i<result_files.size(); i++ ) {
        cmd += "echo ./${result_files[i].name} >> import_file_list.txt\n"
      } 
      cmd
}

/*
 * STEP 9 - set up Proline experimental design
*/

process set_proline_exp_design {
    publishDir "${params.outdir}"
    input:
      val mzdb from mzdbs2.collect()
      file expdesign from input_exp_design.ifEmpty(file("none"))
      
      
//      TODO if file available, change raw file names to mzDB and add experimental design
      
       
    output:
      file "quant_exp_design.txt" into (exp_design_file, exp_design_file2)
    script: 
       // no file provided
      if (expdesign.getName() == "none") {
       expdesign_text = "mzdb_file\texp_condition"
       for( int i=0; i<mzdb.size(); i++ ) {
         expdesign_text += "\n./${mzdb[i].baseName}.mzDB\tMain"
         //cmd += "echo ${mzdb[i].baseName}.mzDB  Main >> quant_exp_design.txt\n"
       } 
       // TODOcheck whether all mzdb files are named in the experimental design file
       // For now only copying file
       """
       touch quant_exp_design.txt    
       echo "${expdesign_text}" >> quant_exp_design.txt
       """
      } else {
          """
          cp ${expdesign} quant_exp_design.txt
          """
/*          for ( int i=0; i<mzdb.size(); i++) {
              if (mzdb[i] != "") print "k"
           }     
           log.warn( expdesign.  )*/
//       expdesign.readLines().each { println it }
      }
}

/*
 * STEP 11 - set up Proline config filelist
*/
process run_proline {
    publishDir "${params.outdir}"

    input:
      file rfs from search_engine_result_files2.collect()
      file mzdbs from mzdbs3.collect()
      file lfq_param from input_proline_param
      file import_param from import_files
      file exp_design from exp_design_file
       
    output:
      file "proline_results/*.xlsx" into proline_out

    script:
    mem = " ${task.memory}"
    mem = mem.replaceAll(" ","")
    mem = mem.replaceAll("B","")
    """
    cp -r /proline/* .
    java8 -Xmx${mem} -cp "config:lib/*:proline-cli-0.2.0-SNAPSHOT.jar" -Dlogback.configurationFile=config/logback.xml fr.proline.cli.ProlineCLI \\
          run_lfq_workflow -i="${import_param}"  -ed="${exp_design}" -c="${lfq_param}"
    """
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

       
