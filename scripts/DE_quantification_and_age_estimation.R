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

path_ref_names <- "~/mydatalocal/tp_ngs_nematode/data/tx2gn.csv"
tx2gn <- read.table(path_ref_names, h=T, sep=',', as.is=T) #nomenclature change (to further be able to extract names that are recognized by WormBase)
txi_salmon <- tximport(files,
                       type = 'salmon', 
                       tx2gene = tx2gn           #we output the transcript-level
                       )
head(txi_salmon$counts)




################### We create a DESeq object from txi_salmon #########################


    #specification of the condition: we attribute 2 different classes for wt mutants
    #this will allow us to know which samples will be opposed in the DE analysis

ColData <- data.frame(NAMES, 
                      factor(rep(c("wt", "alg-5(ram2)"), e=3), 
                             levels = c("wt", "alg-5(ram2)"))) #specifying the levels is CRUCIAL, otherwise, R will set the alphabetic order as default order...

colnames(ColData) <- c("samples", "treatment")

dds <- DESeqDataSetFromMatrix(round(txi_salmon$counts), 
                              colData = ColData , 
                              design = ~ treatment)

quantif <- DESeq(dds, test = "Wald")
res <- results(quantif) #extracts log2(fold_change) and form DESeq output
res

plot(x = res$log2FoldChange, y = -log10(res$pvalue), ylim = c(0,4)) #raw plot
abline(h = -log10(0.05), col = 'red') #significant p-values
abline(v=-1.5, col = 'green')
abline(v = 1.5, col = 'green')

plotMA(res, ylim = c(-12,12))  #nice plot

plot(res$padj)
plot(res$log2FoldChange)

### Extract the output table of DESeq ###

pval_threshold = 0.05
fold_change_threshold = 1.5  #article


  #we split the significantly up and downregulated genes (fold change > thr and < -thr)

ind <- na.omit(rownames(res)[ (res$padj <= pval_threshold) & (abs(res$log2FoldChange) > fold_change_threshold)])
length(ind) #50
cat(ind)

ind_up <- na.omit(rownames(res)[ (res$padj <= pval_threshold) & (res$log2FoldChange > fold_change_threshold)])
length(ind_up) #13
cat(ind_up)

ind_down <- na.omit(rownames(res)[ (res$padj <= pval_threshold) & (res$log2FoldChange < -fold_change_threshold)])
length(ind_down) #37
cat(ind_down)

GO = file("~/mydatalocal/tp_ngs_nematode/results/genes_of_interest.data", "w")
cat(ind, sep='\n', file = GO) #save ind


#################### Age pseudo-estimation ########################

library(RAPToR)
library(limma)

#quantile normalisation and log-transformation of the expression data

r_larv <- prepare_refdata("Cel_larval", "wormRef", n.inter = 600)
ae_algram2 <- ae(samp = txi_salmon$lognorm,                  # input gene expression matrix
                 refdata = r_larv$interpGE,            # reference gene expression matrix
                 ref.time_series = r_larv$time.series)
plot(ae_algram2, group = ColData$treatment, show.boot_estimates = TRUE)
r_larv <- prepare_refdata("Cel_larval_YA", "wormRef", n.inter = 600) #Young Adult Larval reference

#age estimation
ae_algram2 <- ae(samp = txi_salmon$lognorm,                  # input gene expression matrix
                 refdata = r_larv$interpGE,            # reference gene expression matrix
                 ref.time_series = r_larv$time.series)
colors <- c("orange", "green")
plot(ae_algram2, group = ColData$treatment, color = colors[ColData$treatment],show.boot_estimates = TRUE)
