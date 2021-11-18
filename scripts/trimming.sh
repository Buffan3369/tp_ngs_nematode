#DATA CLEANING (TRIMMING)

#maybe try several parameters, partition in the group
cd /home/rstudio/mydatalocal/tp_ngs_nematode/results
mkdir -p trimming #ne râle pas si le fichier existe déjà
cd trimming

trimmomatic_path="/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"  #variable containing the path towards trimmomatic
adapter_path="/softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"


for file in ls /home/rstudio/mydatalocal/tp_ngs_nematode/data/*_1.fastq.gz
do
  #on ne veut garder que le nom du préfixe
  prefix=${file/"_1.fastq.gz"/} #on enlève la fin du nom de file
  prefix=${prefix/"/home/rstudio/mydatalocal/tp_ngs_nematode/data/"/} #et son début
  java -jar $trimmomatic_path PE \
  /home/rstudio/mydatalocal/tp_ngs_nematode/data/${prefix}_1.fastq.gz \
  /home/rstudio/mydatalocal/tp_ngs_nematode/data/${prefix}_2.fastq.gz \
  ${prefix}_forward_paired.fq.gz \
  ${prefix}_forward_unpaired.fq.gz \
  ${prefix}_reverse_paired.fq.gz \
  ${prefix}_reverse_unpaired.fq.gz \
  ILLUMINACLIP:${adapter_path}:2:30:10:2:keepBothReads \
  LEADING:8 TRAILING:8 SLIDINGWINDOW:4:15 MINLEN:36
done

#WE CHOOSE THE FOLLOWING ARGUMENTS ACCORDING TO THE MULTIQC OUTPUT

#ILLUMINACLIP: cuts the remaining Illumina sequences
#LEADING: remove leading low quality (below qty 35)
#TRAILING: remove trailing low qty (below qty 35)
#SLIDINGWINDOW: Scan the read with a 15-base wide sliding window, cutting when the average quality per base drops below 30
#MINLEN: drops reads (here) below 140 bases long (normally, all reads are 150 bases-long)
