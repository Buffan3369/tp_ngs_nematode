#DATA QUALITY ASSESSMENT

fastqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/processed_data /home/rstudio/mydatalocal/tp_ngs_nematode/data/*.gz

multiqc -o /home/rstudio/mydatalocal/tp_ngs_nematode/results/processed_data/multiqc processed_data


#DATA CLEANING

#maybe try several parameters, partition in the group
