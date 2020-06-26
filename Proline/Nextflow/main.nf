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
      --fasta                Fasta file for database search
      --lfq_param                      Parameter file for Proline
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
      --run_msgf                        SearchGui runs msgf+ database search
      --run_comet                       SearchGui runs comet database search
      --run_ms_amanda                   SearchGui runs msamanda database search
      --run_myrimatch                   SearchGui runs myrimatch database search
      
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
if (params.run_xtandem == 0  && params.run_msgf == 0 && params.run_comet == 0 && params.run_ms_amanda == 0 && params.run_myrimatch == 0) {
           log.error "No database engine defined. Make sure you have set one of the --run_searchengine options to 1 (searchengine can be xtandem, msgf, comet, ms_amanda, myrimatch)."; exit 1 
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
        .fromPath( params.raws ).into {input_raw; input_raw2}

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
 * Create a channel for Proline lfq parameter file
 */
  input_proline_param =  Channel
        .fromPath( params.lfq_param)



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
 * STEP 1.5 - convert raw files to mzML
 */
process convert_raw_mzml {
    publishDir "${params.outdir}"
    input:
      file rawfile from input_raw2

    output:
     file("${rawfile.baseName}.mzML") into mzmls

    script:
     """
     ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 2 -z 
     """

}


/*
 * STEP 1.7 - convert raw files to mzDB
 */
process convert_mzml_mzdb {
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
    mem = " ${task.memory}"
    mem = mem.replaceAll(" ","")
    mem = mem.replaceAll("B","")
     """
        mkdir tmp
        mkdir log    
        unzip ${search_out} searchgui.par
        peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files "./${mgffile}"  -identification_files "./${search_out}"  -id_params searchgui.par \\
        -experiment "${params.name}" -sample "${mgffile.baseName}" -out "./${mgffile.baseName}.cpsx" -replicate 1 -threads ${task.cpus} -Xmx${mem}
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
      file("${pepshaker.baseName}.mzid") into (peptideshaker_mzids, peptideshaker_mzids2)

    script:
     """
    peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.MzidCLI -in "./${pepshaker}" -output_file "./${pepshaker.baseName}.mzid" \\
           -contact_first_name Anonyomous -contact_last_name Nextflow -contact_email 'whoami@gmail.com' -contact_address "Greenland" \\
           -organization_name Illuminati -organization_email 'xyz@hell.xyz' -organization_address "California"
         """    
}


/*
 * STEP 8 - set up Proline config filelist
*/

process set_proline_inputfiles {
    publishDir "${params.outdir}"
    input:
      val mzid from peptideshaker_mzids.collect()
       
    output:
      file "import_file_list.txt" into import_files

    script:
      def cmd = 'touch import_file_list.txt\n'
      for( int i=0; i<mzid.size(); i++ ) {
        cmd += "echo ./${mzid[i].baseName}.mzid >> import_file_list.txt\n"
      } 
      cmd
}

/*
 * STEP 9 - set up Proline experimental design
*/

process set_proline_exp_design {
    publishDir "${params.outdir}"
    input:
      val mzdb from mzdbs.collect()
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


/* TODO
 * STEP 10 - set up Proline lfq parameters, for now, just read in file

process set_proline_quant_parameters {
    publishDir "${params.outdir}"
    input:
      val mzid from peptideshaker_mzids.collect()
       
    output:
      file "import_file_list.txt" into filenames

    script:
      def cmd = 'touch import_file_list.txt\n'
      for( int i=0; i<mzid.size(); i++ ) {
        cmd += "echo ${mzid[i].baseName}.mzid >> import_file_list.txt\n"
      } 
      cmd
}
*/

/*
 * STEP 11 - set up Proline config filelist
*/
process run_proline {
    publishDir "${params.outdir}"

    input:
      file mzid from peptideshaker_mzids2.collect()
      file mzdb from mzdbs2.collect()
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
    /usr/lib/jvm/adoptopenjdk-8-hotspot-amd64/bin/java -Xmx${mem} -cp "config:lib/*:proline-cli-0.2.0-SNAPSHOT.jar" -Dlogback.configurationFile=config/logback.xml fr.proline.cli.ProlineCLI \\
            run_lfq_workflow -i="${import_param}"  -ed="${exp_design}" -c="${lfq_param}"
    """
}

/*
 * STEP 11 - run PolySTest for stats
*/
process run_polystest {
    publishDir "${params.outdir}"

    when:
      params.run_statistics

    input:
      file exp_design from exp_design_file2
      file proline_res from proline_out
       
    output:
      file "polystest_prot_res.csv"  into polystest_prot_out
      file "polystest_pep_res.csv" into polystest_pep_out

    script:
    """
    convertFromProline.R ${exp_design} ${proline_res}
    sed -i "s/threads: 2/threads: ${task.cpus}/g" pep_param.yml
    sed -i "s/threads: 2/threads: ${task.cpus}/g" prot_param.yml
    runPolySTestCLI.R pep_param.yml
    runPolySTestCLI.R prot_param.yml    
    """
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

       