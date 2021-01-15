exp_design <-   read.csv("exp_design.txt", sep="\t")
exp_design <- exp_design[order(exp_design[,1]),]
rownames(exp_design) <- exp_design[,1]
# adding 3rd columns with indices for replicate numbers
exp_design <- data.frame(exp_design, replicates=1)
counts <- vector("numeric", length(unique(exp_design[,2])))
names(counts) <- unique(exp_design[,2])
for (exp in 1:nrow(exp_design)) {
  counts[exp_design[exp, 2]] <- counts[exp_design[exp, 2]] + 1
  exp_design$replicates[exp] <- counts[exp_design[exp, 2]]
}  


## Protein table
# read and merge all output files from StPeter output converted to csv
quant_out <- list()
keep_columns_once <- c("protein_name","n_indistinguishable_proteins","indistinguishable_protein","protein_description")
keep_columns_all <- c("probability","unique_stripped_peptides","total_number_peptides", "total_number_distinct_peptides", "pct_spectrum_ids","analysis","SI","SIn")
# to map protein groups to their info
prot_info <- NULL

for (file in exp_design[,1]) {
  t_quant <- read.csv(sub(".raw", ".pep.interact.pep.prot_stpeter.prot_prot.csv", file, fixed = T))
   rownames(t_quant) <-  apply(t_quant[, keep_columns_once], 1, paste, collapse=";")
  t_prot_info <- t_quant[, keep_columns_once]
  quant_out[[file]] <- t_quant[,keep_columns_all]
  prot_info <- rbind(prot_info, t_prot_info)
  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
}
prot_info <- unique(prot_info)
all_quant <- Reduce(function(x, y) merge(x, y, by=0, all=TRUE), quant_out)
all_quant <- cbind(all_quant[,1], prot_info[all_quant[,1],], all_quant[,2:ncol(all_quant)])
write.csv(all_quant, "all_prot_quant_merged.csv", row.names = F)

## Peptide table (left out as tables contain PSMs and thus cannot be easily merged)
#quant_out <- list()
#keep_columns_once <- c("modified_peptide","protein","Quantification_SI","peptide_prev_aa","peptide_next_aa","modifications","peptide","num_missed_cleavages","calc_neutral_pep_mass")

#keep_columns_all <- c("spectrum","spectrumNativeID","assumed_charge","precursor_neutral_mass","retention_time_sec","start_scan","end_scan","xcorr","deltacn","deltacnstar","spscore","sprank","expect","num_matched_ions","tot_num_ions","massdiff","num_matched_peptides","fval","ntt","nmc","massd","peptideprophet_probability","peptideprophet_ntt_prob")

# to map protein groups to their info
#prot_info <- NULL
#for (file in exp_design[,1]) {
#  t_quant <- read.csv(sub(".raw", ".pep.interact.pep.prot_stpeter.prot_pep.csv", file, fixed = T))
#   rownames(t_quant) <-  apply(t_quant[, c("modified_peptide","assumed_charge")], 1 , paste, collapse=";")
#   print(keep_columns_once %in% colnames(t_quant))
#   print(keep_columns_all %in% colnames(t_quant))
#  t_prot_info <- t_quant[, keep_columns_once]
#  quant_out[[file]] <- t_quant[,keep_columns_all]
#  prot_info <- rbind(prot_info, t_prot_info)
#  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
#}
#prot_info <- unique(prot_info)
#all_quant <- Reduce(function(x, y) merge(x, y, by=0, all=TRUE), quant_out)
#all_quant <- cbind(all_quant[,1], prot_info[all_quant[,1],], all_quant[,2:ncol(all_quant)])
#write.csv(all_quant, "all_pep_quant_merged.csv", row.names = F)

