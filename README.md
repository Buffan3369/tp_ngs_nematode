# TP NGS 2021 NÉMATODE

## Téléchargement des données de RNAseq

On télécharge les données sur [GEO database](https://www.ncbi.nlm.nih.gov/geo/) à l'aide de `fastqdump`.On se focalise sur les mutants ram2.
Sur le site, on lit que les données sont paired-end, on prend alors garde à bien split les output files, correspondant aux extrémités forward et reverse. (*cf. script "download_data.sh"*).

La qualité des données est testée à l'aide du module `fastqc`. L'analyse comparée entre les différents individus est implémentée à l'aide du module `multiqc`. (*cf. script "data_quality.sh"*)

On performe ensuite un débroussaillage des données téléchargées à l'aide de l'outil `trimmomatic`. (*cf. script "trimming.sh"*) Cet outil va éliminer les primers Illumina résiduels des données de séquençage, les reads de faible qualité et les reads trop petits. À partir des analyses de qualité menées précédemment, on décide d'éliminer les portions de qualité inférieure à 35 (`LEADING:35` and `TRAILING:35`). Le reste des paramètres de `trimmomatic` sont laissés par défaut.

Enfin, on réitère les analyses de qualité sur l'output de trimmomatic à l'aide de `fastqc`et `multiqc` (*script "trimming_quality.sh"*). 

## Mapping sur transcriptome de référence et quantification

On télécharge et unzip le transcriptome de référence de *C.elegans* sur [Ensembl](http://ftp.ensembl.org). (*cf. script "whole_trasncriptome*).

On procède ensuite au mapping de nos reads sur le transcriptome de référence. Cette étape offre alors la possibilité de quantifier les reads qui associés à une même séquence. Pour cela, on utilise l'outil [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). (*cf script "sushi.sh"*)


