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


# accumulate all descriptions
prot_descr <- NULL
prot_descr <- unique(prot_descr)
rownames(prot_descr) <- prot_descr[,1]

# read and merge all output files from StPeter
quant_out <- list()
for (file in exp_design[,1]) {
  quant_out[[file]] <- read.csv(sub(".raw", ".pep.interact.pep_nsi.csv", file, fixed = T))
  prot_descr <- rbind(prot_descr, quant_out[[file]][,1:2])
  quant_out[[file]] <- quant_out[[file]][,c(1,3:ncol(quant_out[[file]]))]
  colnames(quant_out[[file]]) <- paste(colnames(quant_out[[file]]), exp_design[file,2], exp_design[file,3], sep="_")
}
all_quant <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), quant_out)
all_quant <- cbind(all_quant[,1], prot_descr[all_quant[,1],2], all_quant[,2:ncol(all_quant)])
write.csv(all_quant, "all_quant_merged.csv", row.names = F)

