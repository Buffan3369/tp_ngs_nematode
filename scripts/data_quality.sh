#DATA QUALITY ASSESSMENT

fastqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/processed_data /home/rstudio/mydatalocal/tp_ngs_nematode/data/*.gz

multiqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/processed_data/multiqc /home/rstudio/mydatalocal/tp_ngs_nematode/results/processed_data