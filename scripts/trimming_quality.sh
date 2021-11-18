cd /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimming

mkdir -p /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimmed_data
mkdir -p /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimmed_data/multiqc

#we only keep the paired reads

for file in /*_paired.fq.gz
do
  fastqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimmed_data \
  ${file}
done

multiqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimmed_data/multiqc \
/home/rstudio/mydatalocal/tp_ngs_nematode/results/*_paired.fq.gz