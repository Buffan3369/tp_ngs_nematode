mkdir -p /home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment/multiqc
final=/home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment/multiqc

multiqc -o $final \
/home/rstudio/mydatalocal/tp_ngs_nematode/results/alignment/*_quant