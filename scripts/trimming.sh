#DATA CLEANING (TRIMMING)

#maybe try several parameters, partition in the group
cd /home/rstudio/mydatalocal/tp_ngs_nematode/results
mkdir trimming
cd trimming

trimmomatic_path="/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"  #variable containing the path towards trimmomatic

for file in ls /home/rstudio/mydatalocal/tp_ngs_nematode/data/*.gz
do
  java -jar $trimmomatic_path PE \
  file \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:35 TRAILING:35 SLIDINGWINDOW:15:30 MINLEN:140
done

#WE CHOOSE THE FOLLOWING ARGUMENTS ACCORDING TO THE MULTIQC OUTPUT

#ILLUMINACLIP: cuts the remaining Illumina sequences
#LEADING: remove leading low quality (below qty 35)
#TRAILING: remove trailing low qty (below qty 35)
#SLIDINGWINDOW: Scan the read with a 15-base wide sliding window, cutting when the average quality per base drops below 30
#MINLEN: drops reads (here) below 140 bases long (normally, all reads are 150 bases-long)