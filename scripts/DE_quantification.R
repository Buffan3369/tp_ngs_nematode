library(DESeq2)
library(tximport)

#no need to compute a loop as tximport directly recognizes when there are several files

dir <- "/home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment"

NAMES <- c("SRR5564855",
           "SRR5564856",
           "SRR5564857",
           "SRR5564867",
           "SRR5564868",
           "SRR5564869")


#### we specify the NAMES of the final directories ####

Dir_fin <- c() #list of the 6 final directories
for(name in NAMES){
  dir_fin <- paste(name,"_quant","/quant.sf", sep = '') #name of the dir "name_quant/quant.sf"
  Dir_fin <- c(Dir_fin, dir_fin)
}

files <- file.path(dir, Dir_fin) #the NAMES of the salmon outputs


#### we then run tximport ####

txi_salmon <- tximport(files,
                       type = 'salmon', 
                       txOut = TRUE           #we output the transcript-level
                       )
head(txi_salmon$counts)


#### We create a DESeq object from txi_salmon ####

    #specification of the condition: we attribute 2 different classes for wt mutants
    #this will allow us to know which samples will be opposed in the DE analysis

ColData <- data.frame(NAMES, rep(c("wt", "alg-5(ram2)"), e=3))
colnames(ColData) <- c("samples", "treatment")

dds <- DESeqDataSetFromMatrix(round(txi_salmon$counts), colData = ColData , design = ~ treatment)

quantif <- DESeq(dds, test = "Wald")
res <- results(quantif) #extracts log2(fold_change) and form DESeq output

plot(x = res$log2FoldChange, y = log(res$pvalue)) #raw plot
plotMA(res)