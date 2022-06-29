# Stage M2 : Développements de méthodes d’analyse de l’épissage alternatif à partir de données RNA-seq

## Présentation

Ce projet contient contient les résultats de mon stage de M2 portant sur le développement de méthodes d’analyse de l’épissage alternatif à partir de données RNA-seq. Un rapport de stage est disponible sur ce git. L'algorithme principal de ce Git permet de simplifier les graphes de De Bruijn afin de trouver de nouveaux évènements d'épissage alternatif.

Voici une présentation du contenu de ce Git :

### Data :

Ce dossier contient un mini jeu de données pour tester les algorithmes mis en place. Ce jeu de donnée "test" correspond à un graphe de De Bruijn compacté et est décomposé en trois fichiers :

- `test.nodes` : Chaque ligne correspond à un unique sommet qui est défini par un numéro et une séquence de nucléotide.
- `test.edges` : Chaque ligne correspond à une unique arête qui est définie par un sommet de départ, un sommet d'arrivé et un type d'arête (plus d'explications sur les types des arêtes à la structure du graphe sont disponibles dans le rapport).
- `test.abundance` : Chaque ligne correspond à l'abondance des k-mers de chaque sommet. Je n'utilise pas cette métrique, elle ne sert qu'au bon fonctionnement de *KisSplice*.

### Processing :

- `test.nodes.pondere` : Ce fichier correspond au fichier `test.nodes` mais avec en plus pour chaque sommet son poids. Il a était obtenu à la suite de l'étape 1 de l'algorithme.
- `comp0.txt` : Ce fichier contient tous les sommets faisant parties de la composante 0. Il a été généré à la suite de l'étape 2.

### Results :

Ce dossier contient les résultats finaux de mon alogithme. Il est composé de fichiers `clean.nodes`, `clean.edges` et `clean.abundance` qui correspondent au nouveau graphe qui sera par la suite utilisé par *KisSplice*.


### Scr :

Ce dossier contient tous les algorithmes. On peut les décomposer en trois catégories :

Tout d'abord, il y a les programmes de l'algorithme de simplification : 

- `graph.cpp` et `graph.h` : Contiennent toutes les structures et fonctions liées aux graphes que j'ai utilisés.
- `ponderation.cpp` : Etape 1 de l'algorithme permettant de faire la pondération du graphe.
- `gene_comp.cpp` : Etape 2 de l'algorithme permettant de générer toutes les composantes.
- `composantes.cpp` : Etape 3 de l'algorithme permettant de simplifier effectivement le graphe.

Ensuite, il y a les programmes liés à l'analyse de mes données ou de mes méthodes :

- `empty_composantes.cpp` : 
- `filtrage_bulle.py` : 
- `interBulle.py` : 
- `neighborhood.cpp` : 
- `plot.py` : 
- `rapportAgglo` : 
- `reads_to_align.py` : 
- `suppDoublon.py` : 

Finalement, il y a les programmes de pipeline :

- `start.nf` : Permet d'exécuter l'algorithme sur n'importe quel graphe en un ligne. Par défault, s'exécute uniquement sur le graphe `test`.
- `analyse_TE` : Permet d'effectuer l'analyse complète des éléments transposables des jeux de données. Nécéssite d'avoir des fichiers d'annotation pour fonctionner.


## Exemple d'utilisation 

### Etape 1 : la podération du graphe

### Etape 2 : la génération des composantes

### Etape 3 : la simplification du graphe de De Bruijn




## Ancien Protocole utilisé pour la mise en place de la pipeline


### Etape 1: Téléchargement

On a besoin tout d'abord des séquences ref_genome.fna, genome.gtf et TEgenome.gtf. Ensuite, on va utiliser plusieurs programme dans la pipeline :
- STAR
- KisSplice


### Etape 2: Création de la carte

Moustique :

```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ../results/moustique/genome --genomeFastaFiles ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna --sjdbGTFfile ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf -sjdbOverhang 74 --genomeSAindexNbases 13 --genomeSAsparseD 4
```

Path pour le chien : home/sdarmon/Documents/stage-M2/peda/fastq/SRR_Chiens

```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ../results/chien/ --genomeFastaFiles ../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-unmasked.fa --sjdbGTFfile ../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-2021_03-genes.gtf -sjdbOverhang 74 --genomeSAindexNbases 13 --genomeSAsparseD 4 --readFilesCommand zcat
```

### Etape 3: Génération du graphe de De Bruijn

```
./build/bin/kissplice -r read.fa 
```
Attention, il s'agit là d'une solution over-killed, il vaudra mieux repartir de `BCaml` ! (Et attention, être dans le bon répertoire, pour moi data)

### Etape 4: Calcul des poids des sommets

On calcule les poids de chaque sommet (rajouter ici une explication)

Version moustique:
```
g++ graph.cpp main.cpp -o graph.exe
 cd ../../../kissplice_results/kissplice_moustique
../../DmGoth/stage-M2/scr/graph.exe graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges 200 -o ../../DmGoth/data/outputGraphMoustique.txt 
cd ../../DmGoth/stage-M2/scr
```

Version moustique version non optimisée:
```
g++ graph.cpp main.cpp -o graph.exe
 cd ../../../kissplice_results/kissplice_moustique
../../DmGoth/stage-M2/scr/graph.exe graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges 5 -o ../../DmGoth/data/outputGraphMoustiqueNonOpt.txt 
cd ../../DmGoth/stage-M2/scr
```


version chien:

```
g++ graph.cpp main.cpp -o graph.exe
./graph.exe ../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes ../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges 200 -o ../../data/chien/outputGraphChien.txt 
```


Puis on lance l'analyse graphique (dot) et on vérifie que la valeur donnée pour le threshold semble correcte (top10)

```
python3 plot.py ../../data/outputGraphMoustique.txt top10
python3 plot.py ../../data/outputGraphMoustique.txt dot
```

Version moustique version non optimisée:

```
python3 plot.py ../../data/outputGraphMoustiqueNonOpt.txt top10
python3 plot.py ../../data/outputGraphMoustiqueNonOpt.txt dot
```
Version Chien:

```
python3 plot.py ../../data/chien/outputGraphChien.txt top10
python3 plot.py ../../data/chien/outputGraphChien.txt dot
```

Finallement, on génère un fichier de reads à aligner sur le génome de référence. (Ici threshold = 11 pour moustique, 9 pour top 20;  16 pour NonOpt, 11 pour top20, 8 pour chien200, ? pour chien300, 15 moustique300, 20 top1 moustique, 30 top 0.1 moustique)

```
python3 reads_to_align.py ../../data/outputGraphMoustique.txt ../../data/readMoustique.fq 11
```

Version moustique version non optimisée:

```
python3 reads_to_align.py ../../data/outputGraphMoustiqueNonOpt.txt ../../data/readMoustiqueNonOptClean.fq 16
```

version chien:

```
python3 reads_to_align.py ../../data/chien/outputGraphChien.txt ../../data/chien/readChien.fq 8
```

Où the input is the output of graph.exe, output.fq est le format des séquences à générer et threshold est l'output de plot.py top10.

### Etape 5: Alignement de nos reads sur le génome de ref

```
cd ../../data
STAR --genomeDir ../results \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn readMoustique.fq \
--outFileNamePrefix ../results/STAR/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 10000 \
--outReadsUnmapped Fastx
```

Version moustique version non optimisée:
```
cd ../../data
STAR --genomeDir ../results \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn readMoustiqueTop1.fq \
--outFileNamePrefix ../results/moustique/STARClean/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 10000 \
--outReadsUnmapped Fastx
```

Chien:

```
cd ../../data/chien
STAR --genomeDir ../../results/chien \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn readChien.fq \
--outFileNamePrefix ../../results/chien/STAR/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 10000 \
--outReadsUnmapped Fastx
```


### Etape 6: Intersection des reads alignés avec les TE connues

```
bedtools intersect -wa -a AaegL5_TE_repeats.gff -b ../results/STAR/Aligned.sortedByCoord.out.bam > ../results/moustique/intersectionTE.txt
bedtools intersect -wb -a AaegL5_TE_repeats.gff -b ../results/STAR/Aligned.sortedByCoord.out.bam > ../results/moustique/intersectionKiss.txt
cd ../stage-M2/scr/
python3 suppDoublon.py ../../results/moustique/intersectionKiss.txt ../../results/moustique/intersectionKissNoDouble.txt -s 12
python3 suppDoublon.py ../../results/moustique/intersectionTE.txt ../../results/moustique/intersectionTENoDouble.txt -t 8
wc -l ../../results/moustique/intersectionKissNoDouble.txt 
wc -l ../../results/moustique/intersectionTENoDouble.txt 
less ../../results/STAR/Log.final.out
```


Version CLEAN
```
bedtools intersect -wb -b ../../results/STAR/Aligned.sortedByCoord.out.bam -a ../../data/moustique/AaegL5_TE_repeats.gff  > ../../results/moustique/intersectionKissClean.txt
bedtools intersect -wa -b ../../results/STAR/Aligned.sortedByCoord.out.bam -a ../../data/moustique/AaegL5_TE_repeats.gff  > ../../results/moustique/intersectionTEClean.txt
cd ../stage-M2/scr/
python3 suppDoublon.py ../../results/moustique/intersectionKissClean.txt ../../results/moustique/intersectionKissCleanNoDouble.txt -s 12
python3 suppDoublon.py ../../results/moustique/intersectionTEClean.txt ../../results/moustique/intersectionTECleanNoDouble.txt -t 8
wc -l ../../results/moustique/intersectionKissCleanNoDouble.txt 
wc -l ../../results/moustique/intersectionTECleanNoDouble.txt 
less ../../results/moustique/STARClean/Log.final.out
```


Version moustique version non optimisée: A FAIRE DIRECTEMENT SUR SERVEUR DEMANDE 24Go DE RAM!!!
```
bedtools intersect -wa -a AaegL5_TE_repeats.gff -b ../results/moustique/STAR/Aligned.sortedByCoord.out.bam > ../results/moustique/intersectionTENonOpt.txt
bedtools intersect -wb -a AaegL5_TE_repeats.gff -b ../results/moustique/STAR/Aligned.sortedByCoord.out.bam > ../results/moustique/intersectionKissNonOpt.txt
cd ../stage-M2/scr/
python3 suppDoublon.py ../../results/moustique/intersectionKissNonOpt.txt ../../results/moustique/intersectionKissNonOptNoDouble.txt -s 12
python3 suppDoublon.py ../../results/moustique/intersectionTENonOpt.txt ../../results/moustique/intersectionTENonOptNoDouble.txt -t 8
wc -l ../../results/moustique/intersectionKissNonOptNoDouble.txt 
wc -l ../../results/moustique/intersectionTENonOptNoDouble.txt 
less ../../results/moustique/STAR/Log.final.out
```

## Différence entre les deux méthodes de distance :

### intersection Ref vs nonOpt:

```
  cd ../../results/moustique
  python3 ../../stage-M2/scr/interseq.py intersectionTEAllSeq.txt  intersectionTENonOptNoDouble.txt nonOpt
wc -l nonOptintersection.txt 
wc -l nonOptseq1Remaning.txt 
wc -l nonOptseq2Remaning.txt
mv nonOptseq1Remaning.txt nonOptseq1Remaning.gff
bedtools intersect -wb -a nonOptseq1Remaning.gff -b ../STAR/Aligned.sortedByCoord.out.bam > seqKiss.txt
python3 ../../stage-M2/scr/reads_to_align.py ../../data/outputGraphMoustique.txt seqMissing.txt 0 -reverse seqKiss.txt
```

### Intersection 200 vs 5:

```
cd ../../results/moustique
python3 ../../stage-M2/scr/interseq.py intersectionTENoDouble.txt  intersectionTENonOptNoDouble.txt mous
wc -l mousintersection.txt 
wc -l mousseq1Remaning.txt 
wc -l mousseq2Remaning.txt
mv mousseq1Remaning.txt mousseq1Remaning.gff
mv mousseq2Remaning.txt mousseq2Remaning.gff
mv mousintersection.txt mousintersection.gff
bedtools intersect -wb -a mousseq1Remaning.gff -b ../STAR/Aligned.sortedByCoord.out.bam > seqKiss5.txt
bedtools intersect -wb -a mousseq1Remaning.gff -b STAR/Aligned.sortedByCoord.out.bam > seqKiss200.txt
bedtools intersect -wb -a mousintersection.gff -b STAR/Aligned.sortedByCoord.out.bam > seqKissInter.txt
python3 ../../stage-M2/scr/reads_to_align.py ../../data/outputGraphMoustiqueNonOpt.txt seqMousMissing5.txt 0 -reverse seqKiss5.txt
python3 ../../stage-M2/scr/reads_to_align.py ../../data/outputGraphMoustique.txt seqMousMissing200.txt 0 -reverse seqKiss200.txt
python3 ../../stage-M2/scr/reads_to_align.py ../../data/outputGraphMoustique.txt seqMousMissingInter.txt 0 -reverse seqKissInter.txt
```
python3 ../../stage-M2/scr/reads_to_align.py ../../data/outputGraphMoustique.txt seqTop1Inter.txt 30 -reverse intersectionKissTop1NoDouble.txt

## Premiers pas avec awk:

La première ligne permet de récupérer les lignes d'un fichier dont la 3ème item est supérieur à 2000.

Le seconde ligne permet d'obtenir l'abonance de l'unitig 6991943 (attention à bien faire un plus un, car à la 1ère ligne il y a l'unitig 0).
```
awk '$3 > 2000' outputGraph300.txt 
awk -v ligne=6991944 ' NR == ligne { print $0}' ../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.abundance 
sort -n -k 3 ../../data/outputGraphMoustiqueTop1Clean.txt
```

TTTTTTTTTTTGGATTTTTGGTATAATCTTAGAGAAATTTGCGGAAGAAATGTAGAGGGATCTCAGCAAGAATCCCAGAAGTATTATTAGAAGAATCCCAGAGAGTTTTGGATTTCCAGTAGAATTCCACATGAATTCCAGAAGAATCCAAAAGGAATTACCGCAAGTCCCCTAGCACCGCAACTAACATTTTAGTAACATATATAACACACGCTACGTACACAAAAGTTACATGAATGATGTGCAATGCTCAAAACATAAGTGCAAATTGCGGCAAGCTGAG


## Agglomeration des graphes

```
python3 reads_to_align.py ../../data/outputGraphMoustique.txt ./../data/outputClean.txt 0 -clean
g++ graph.cpp agglo.cpp -o agglo.exe
./agglo.exe ../../data/outputGraphMoustique.txt ../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges -c 30 ../../results/moustique/compMous
```


version chien
```
python3 reads_to_align.py ../../data/chien/outputGraph.txt ../../data/outputCleanChien.txt 0 -clean
g++ graph.cpp agglo.cpp -o agglo.exe
./agglo.exe ../../data/outputCleanChien.txt ../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges -c 25 ../../results/moustique/compChien

```
