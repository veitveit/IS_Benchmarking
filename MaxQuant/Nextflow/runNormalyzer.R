# merging experimental design file from sdrf parser with actual design
A<-read.csv("exp_file.tsv",sep="\t")
B<- read.csv("exp_file2.tsv", sep="\t")

B[,1] <- sub(".raw","",B[,1])
names(B) <- c("file","group")
A$Experiment <- make.names(A$Experiment);
C<-merge(A,B,by.x="Name",by.y="file");
final_exp <- unique(C[,c("Experiment","group")])
num_cond <- nrow(final_exp)
write.table(final_exp, "full_design.txt",sep="\t", row.names=F)

# run Normalyzer
NormalyzerDE::normalyzer(jobName="Project", designPath="full_design.txt", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F)
# comparison set to everything versus first
comps <- paste0(1,"-",2:num_cond)
#NormalyzerDE::normalyzerDE(jobName="Project", comparisons=comps, designPath="full_design.txt", dataPath="./Project/TODONormalyzerMethod-normalized.txt", outputDir="./")'
