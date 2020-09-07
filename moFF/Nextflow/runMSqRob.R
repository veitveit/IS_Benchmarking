library(MSnbase)
library(MSqRob)
options(stringsAsFactors=T)
peptides <- MSqRob::import2MSnSet("moff.tab",filetype="moFF")
#head(exprs(peptides))
runs <- colnames(exprs(peptides))
exp_annotation <- read.csv("exp_design.tsv",sep="\t")
exp_annotation$genotype <- make.names(exp_annotation$genotype)
for (c in 1:ncol(exp_annotation)) {
  exp_annotation[,c] <- as.factor(exp_annotation[,c])
}

#print(exp_annotation)
peptides <- preprocess_MSnSet(peptides,accession="prot",split=",", useful_properties="peptide", exp_annotation=exp_annotation)
#head(fData(peptides))
# necessary due to change to R 4.x    
fData(peptides)[,"prot"] <- as.factor(fData(peptides)[,"prot"])
# Set data.frame for experimental design with columns run genotype biorep
proteins <- MSnSet2protdata(peptides, accession="prot")  
      #head(proteins)
system.time(protLM <- fit.model(proteins, response="quant_value", fixed=c("genotype"),  random=c("biorep","run","peptide"), add.intercept=TRUE))
          #create comparisons vs first
contrasts <- NULL
levels <- as.factor(paste0("genotype",make.names(unique(exp_annotation$genotype))))
# levels <- unique(exp_annotation$genotype)
if (length(levels) > 1) {
  for (i in levels[-1]) {
    contrasts <- append(contrasts, paste(i, levels[1], sep="-"))  
  }
} else {
  contrasts <- levels
}
L <- makeContrast(contrasts=contrasts,levels=as.character(levels))
    result <- test.protLMcontrast(protLM, L)
  result <- prot.p.adjust(result, method="fdr")
      write.csv(result, "MSqRobOut.csv")
