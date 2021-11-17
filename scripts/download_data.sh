
#We first list the accession numbers of the sequences we want to download (RAM)

cd /home/rstudio/mydatalocal/tp_ngs_nematode/data   #specification of the path

accession_numbers="SRR5564867
                    SRR5564868
                    SRR5564869
                    SRR5564855
                    SRR5564856
                    SRR5564857"
                    

for SRA in $accession_numbers
do
  fastq-dump -X 1000 --split-3 --gzip $SRA #--gzip wil compresses the output file, --split-3 will put the output of the forward and reverse sequences in 2 different files,
done

# nb: if you want to test your code instead of downloading the whole (heavy) files, -X 1000 will only focus on the first 1000 elements of the files, since they're a bit heavy

#DATA QUALITY ASSESSMENT

fastqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/data/processed /home/rstudio/mydatalocal/tp_ngs_nematode/data/*.gz