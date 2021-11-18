mkdir -p /home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment #directory where we'll store salmon's outputs

cd /home/rstudio/mydatalocal/tp_ngs_nematode/data/reference_transcriptome

ref_transcriptome=Caenorhabditis_elegans.WBcel235.cdna.all.fa

#indexing reference transcriptome

ind=/home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment/ref_idx

salmon index -t $ref_transcriptome \
-i $ind

#pseudo-mapping the transcripts with this ref

Names="SRR5564855
      SRR5564856
      SRR5564857
      SRR5564867
      SRR5564868
      SRR5564869"

for transcript in $Names
do
  salmon quant \
  --libType A \
  --index $ind \
  --mates1 /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimming/${transcript}_forward_paired.fq.gz \
  --mates2 /home/rstudio/mydatalocal/tp_ngs_nematode/results/trimming/${transcript}_reverse_paired.fq.gz \
  --threads 3 \
  --output /home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment/${transcript}_quant
done

# --threads 3 va utiliser 3 coeurs => parallélisation du programme