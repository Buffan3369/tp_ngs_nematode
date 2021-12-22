# TP NGS 2021 NÉMATODE

## Contexte

[Brown et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587817/) ont travaillé sur les protéines Argonautes de *C.elegans*. Ces protéines s'associent à des miRNA et sont connues pour réguler l'expression d'un grand nombre de gènes et sont assez conservées dans le vivant. Les auteurs se sont plus précisément penchés sur la protéine ALG-5, dont les implications fonctionnelles chez le nématode étaient encore inconnues au temps de l'étude. 

Leur étude comporte deux axes:
 <ol>
  <li> Par RNAseq haut débit, ils ont trouvé qu'ALG-5 se liait à un sous-ensemble de miRNA assez précis, à la différence d'ALG-1 et ALG-2, deux argonautes déjà connues. 
  <li> Par analyse d'expression différentielle chez les mutants *alg-5*, ils identifient les voies dans lesquelles la protéine semble impliquée.
    
 Toutefois, 
      
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
### Expression différentielle

On commence par charger les output counts de `salmon` dans l'environnement R à l'aide du package `tximport`.

On les analyse à l'aide du package `DESeq2`. On obtient les volcano plots.

On extrait les noms des gènes différetiellement exprîmés entre le WT et notre mutant (*alg-5(ram2)*) et on procède à une Enrichment Analysis sur [Wormbase](https://wormbase.org/tools/enrichment). On trouve, avec un seuil de log2FoldChange de 1.5, 50 gènes dont l'expression diffère significativement entre le mutant et le WT. Parmi eux, 37 sont up-régulés (niveau d'expression supérieur chez le mutant par rapport au sauvage) et 13 sont down-régulés.

Tout est indiqué dans le script (*"DE_quantification_and_age_estimation.R"*)

### Estimation d'âge
Pour estimer l'âge développemental de nos échantillons, on utilise l'outil `RAPToR`.
On copare avec les échantillons avec une référence larvaire/jeune adulte pour couvrir un intervalle de temps large.
Les résultats d'estimation d'âge sont présentés ci-dessous.
![Results Age Estimation]()

### Quantification de l'impact du développement sur les résultats obtenus

On cherche à estimer la part du développement expliquant les différences observées entre les 2 groupes (wt et mutant).
On dispose de 2 fonctions:
<ul>
  <li>`getrefTP`: extraire de la référence disponible sur wormRef les estimés de fold change associés au timing développemental
  <li>`refCompare`: comparer nos échantillons à la référence en termes de niveaux d'expression. Elle implémente des modèles linéaires pour expliquer les différences d'expression des gènes par le facteur étudié (traitement: wt/mutant).

  
