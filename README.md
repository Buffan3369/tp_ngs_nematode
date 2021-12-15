# TP NGS 2021 NÉMATODE

## Téléchargement des données de RNAseq

On télécharge les données sur [GEO database](https://www.ncbi.nlm.nih.gov/geo/) à l'aide de `fastqdump`.On se focalise sur les mutants ram2.
Sur le site, on lit que les données sont paired-end, on prend alors garde à bien split les output files, correspondant aux extrémités forward et reverse. (*cf. script "download_data.sh"*).

La qualité des données est testée à l'aide du module `fastqc`. L'analyse comparée entre les différents individus est implémentée à l'aide du module `multiqc`. (*cf. script "data_quality.sh"*)

On performe ensuite un débroussaillage des données téléchargées à l'aide de l'outil `trimmomatic`. (*cf. script "trimming.sh"*) Cet outil va éliminer les primers Illumina résiduels des données de séquençage, les reads de faible qualité et les reads trop petits. À partir des analyses de qualité menées précédemment, on décide d'éliminer les portions de qualité inférieure à 35 (`LEADING:35` and `TRAILING:35`). Le reste des paramètres de `trimmomatic` sont laissés par défaut.

Enfin, on réitère les analyses de qualité sur l'output de trimmomatic à l'aide de `fastqc`et `multiqc` (*script "trimming_quality.sh"*). 



## Mapping sur transcriptome de référence

On télécharge et unzip le transcriptome de référence de *C.elegans* sur [Ensembl](http://ftp.ensembl.org). (*cf. script "whole_trasncriptome"*).

On procède ensuite au mapping de nos reads sur le transcriptome de référence. Cette étape offre alors la possibilité de quantifier les reads qui associés à une même séquence. Pour cela, on utilise l'outil [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). (*cf script "sushi.sh"*)
Pour vérifier comment s'est passé le mapping, on analyse les données de sortie de `salmon` à l'aide de `multiqc` (*code "mapping_output_ana.sh"*)



## Quantification, analyse d'expression différentielle et estimation d'âges
#Expression différentielle

On commence par charger les output counts de `salmon` dans l'environnement R à l'aide du package `tximport`.

On les analyse à l'aide du package `DESeq2`. On obtient les volcano plots.

On extrait les noms des gènes différetiellement exprîmés entre le WT et notre mutant (*alg-5(ram2)*) et on procède à une Enrichment Analysis sur [Wormbase](https://wormbase.org/tools/enrichment). On trouve, avec un seuil de log2FoldChange de 1.5, 50 gènes dont l'expression diffère significativement entre le mutant et le WT. Parmi eux, 37 sont up-régulés (niveau d'expression supérieur chez le mutant par rapport au sauvage) et 13 sont down-régulés.

Tout est indiqué dans le script (*"DE_quantification_and_age_estimation.R"*)

#Estimation d'âge
Pour estimer l'âge développemental de nos échantillons, on utilise l'outil `RAPToR`.
On copare avec les échantillons avec une référence larvaire/jeune adulte pour couvrir un intervalle de temps large.
Les résultats d'estimation d'âge sont présentés ci-dessous.