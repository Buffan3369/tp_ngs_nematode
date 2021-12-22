# TP NGS 2021 NÉMATODE

## Contexte

[Brown et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587817/) ont travaillé sur les protéines Argonautes de *Caenorhabditis elegans*. Ces protéines, assez conservées au sein du vivant, sont connues pour réguler l'expression d'un grand nombre de gènes en s'associant à des miRNA. Les auteurs se sont plus précisément penchés sur la protéine ALG-5, dont les implications fonctionnelles chez le nématode étaient encore inconnues au temps de l'étude. 

Leur étude comporte deux axes:
 <ol>
  <li> Par RNAseq haut débit, ils ont trouvé qu'ALG-5 se liait à un sous-ensemble de miRNA assez précis, à la différence d'ALG-1 et ALG-2, deux argonautes déjà connues. 
  <li> Par analyse d'expression différentielle chez les mutants *alg-5*, ils identifient les voies dans lesquelles la protéine semble impliquée.
 </ol>    
 
Toutefois, à la lumière des avancées actuelles dans le monde de la recherche en génomique, leurs méthodes, et par conséquent leurs résultats, peuvent être remis en question. En effet, *C.elegans* est connu pour avoir des profils d'expression génétique drastiquement changeants au cours du développement: jusqu'à 10000 gènes/heure peuvent voir leur expression modifiée à certains stades, ce qui représente la moitié du nombre de gènes du nématode! Ainsi, il apparaît crucial, quand on mène des analyses d'expression différentielle, de s'assurer que les échantillons que l'on compare (wild type et mutant) ont bien le même âge développemental. Dans le cas contraire (*i.e* en comparant deux vers d'âge développemental différent), notre expérience a plus d'un paramètre qui varie: les changements observés dans l'expression des gènes peuvent effectivement être dûs à la condition test (la mutation perte de fonction de ALG-5) mais également être un artefect du fait que les gènes de *C.elegans* s'expriment différentiellement au cours du temps dans le développement.
Notons cependant qu'il n'est pas exclu (et même hautement probable) que les changements d'expression observés naturellement chez le nématode au cours du développement soient en partie dûs à la protéine ALG-5. Aussi, les deux conditions qui viennent d'être décrites ne sont pas disjointe: la première est très probablement incluse dans la seconde. Mais cette inclusion n'est sûrement pas réciproque, car il n'y a aucune raison pour qu'ALG-5 soit le seul acteur des dynamiques d'expression des gènes au cours du développement chez *C.elegans*.

Or, dans le papier que nous avons introduit, l'alignement des âges des nématodes comparés pour l'analyse d'expression différentielle est assez critiquable. Les âges développementaux ont été estimés de manière qualitative, à la loupe, et un manque de détermination absolue est à souligner. Leur estimation n'exclut pas une fenêtre non-négligeable de plusieurs heures entre deux nématodes qui seront ensuite comparés, or nous avons pu voir qu'en particulier avec un tel modèle d'étude, chaque minute a son importance.

Le but de notre projet est alors de répéter la deuxième partie de leur étude en y incorporant une détermination absolue de l'âge développemental des échantillons à l'aide d'une méthode récente pour estimer la part de variabilité due au développement et réellement due à la mutation dans les profils d'expression différentielle des vers analysés. Ce fichier renvoit à des scripts qui s'intéressent plus particulièrement à la délétion *alg-5 (ram2)*.

      
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
![Results Age Estimation](https://github.com/Buffan3369/tp_ngs_nematode/blob/master/age_estimation.png)

### Quantification de l'impact du développement sur les résultats obtenus

On cherche à estimer la part du développement expliquant les différences observées entre les 2 groupes (wt et mutant).
On dispose de 2 fonctions:
<ul>
  <li>`getrefTP`: extraire de la référence disponible sur wormRef les estimés de fold change associés au timing développemental
  <li>`refCompare`: comparer nos échantillons à la référence en termes de niveaux d'expression. Elle implémente des modèles linéaires pour expliquer les différences d'expression des gènes par le facteur étudié (traitement: wt/mutant).

  
