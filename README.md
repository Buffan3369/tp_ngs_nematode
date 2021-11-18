# TP NGS 2021 NÉMATODE

## Téléchargement des données de RNAseq

On télécharge les données sur [GEO database](https://www.ncbi.nlm.nih.gov/geo/) à l'aide de `fastqdump`.
<br>
On se focalise sur les mutants ram2.
<br>
Sur le site, on lit que les données sont paired-end, on prend alors garde à bien split les output files, correspondant aux extrémités forward et reverse. (*cf. script "download_data.sh"*).

La qualité des données est testée à l'aide du module `fastqc`. L'analyse comparée entre les différents individus est implémentée à l'aide du module `multiqc`. (*cf. script "data_quality.sh"*)

On performe ensuite un débroussaillage des données téléchargées à l'aide de l'outil `trimmomatic`. (*cf. script "trimming"*) Cet outil va éliminer les primers Illumina résiduels des données de séquençage, les reads de faible qualité et les reads trop petits. 

## Mapping sur génome de référence

On télécharge et unzip le génome de référence de *C.elegans* sur [Ensembl](http://ftp.ensembl.org). (*cf. script "whole_genome*)



