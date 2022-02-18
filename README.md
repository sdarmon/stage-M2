# stage-M2



## Utiles

### Pour se connecter en sshfs:

```
sshfs sdarmon@pedago-ngs:/localdata/pandata/students/Projet_KS/DmGoth peda/
```
### Pour copier un fichier sur le serveur

```
scp -p ../../Bureau/genome_assemblies_genome_fasta.tar sdarmon@pedago-ngs:/localdata/pandata/students/Projet_KS/DmGoth/data
```


## proto


### Etape 1: Téléchargement

On a besoin tout d'abord des séquences ref_genome.fna, genome.gtf et TEgenome.gtf. Ensuite, on va utiliser plusieurs programme dans la pipeline :
- STAR
- KisSplice


### Etape 2: Création de la carte

Moustique : 

```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ../results/ --genomeFastaFiles ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna --sjdbGTFfile ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf -sjdbOverhang 74 --genomeSAindexNbases 13 --genomeSAsparseD 4
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
cd ../../../kissplice_moustique
../DmGoth/stage-M2/scr/graph.exe graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges 200 -o ../DmGoth/data/outputGraph.txt 
cd ../DmGoth/stage-M2/scr
```

version chien:

```
g++ graph.cpp main.cpp -o graph.exe
./graph.exe ../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes ../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges 200 -o ../../data/chien/outputGraph.txt 
```


Puis on lance l'analyse graphique (dot) et on vérifie que la valeur donnée pour le threshold semble correcte (top10)

```
python3 plot.py ../../data/outputGraph.txt top10
python3 plot.py ../../data/outputGraph.txt dot
```

Iench version:

```
python3 plot.py ../../data/chien/outputGraph.txt top10
python3 plot.py ../../data/chien/outputGraph.txt dot
```

Finallement, on génère un fichier de reads à aligner sur le génome de référence. (Ici threshold = 11 pour moustique, 8 pour chien200, ? pour chien300)

```
python3 reads_to_align.py ../../data/outputGraph.txt ../../data/read.fq 8
```

version chien:

```
python3 reads_to_align.py ../../data/chien/outputGraph300.txt ../../data/chien/read300.fq 8
```

Où the input is the output of graph.exe, output.fq est le format des séquences à générer et threshold est l'output de plot.py top10.

### Etape 5: Alignement de nos reads sur le génome de ref

```
cd ../../data
STAR --genomeDir ../results \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn read.fq \
--outFileNamePrefix ../results/STAR/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 10000 \
--outReadsUnmapped Fastx
```

iench:

```
cd ../../data/chien
STAR --genomeDir ../../results/chien \
--runMode alignReads \
--runThreadN 8 \
--readFilesIn read.fq \
--outFileNamePrefix ../../results/chien/STAR/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMultimapNmax 10000 \
--outReadsUnmapped Fastx
```


### Etape 6: Intersection des reads alignés avec les TE connues

```
bedtools intersect -wa -a AaegL5_TE_repeats.gff -b ../results/STAR/Aligned.sortedByCoord.out.bam > ../results/intersectionTE.txt
bedtools intersect -wb -a AaegL5_TE_repeats.gff -b ../results/STAR/Aligned.sortedByCoord.out.bam > ../results/intersectionKiss.txt
cd ../stage-M2/scr/
python3 suppDoublon.py ../../results/intersectionKiss.txt ../../results/intersectionKissNoDouble.txt -s 12
python3 suppDoublon.py ../../results/intersectionTE.txt ../../results/intersectionTENoDouble.txt -t 8
wc -l ../../results/intersectionKissNoDouble.txt 
wc -l ../../results/intersectionTENoDouble.txt 
less ../../results/STAR/Log.final.out
```


### Etape 7: Plotting des reads concernés

```
python3 reads_to_align.py ../../data/outputGraph.txt ../../results/seq.txt 11 -reverse ../../results/intersectionKissNoDouble.txt
```


## Afficher une séquence sur Cytoscape:

```
/data/home/vincent/TiffanyDelhomme$ ./nbh -n /localdata/pandata/students/Projet_KS/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes -e /localdata/pandata/students/Projet_KS/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges -k 41 -o /localdata/pandata/students/Projet_KS/DmGoth/data/nbh/graph -d 10 -q AGTAAATGTCACAGTTACAATCTCCGGCCATGGAAAAACCAGAAGTAATGAGCCGATCTCCGAGATGCTCAAGTGCAACACC

python3 labellingCytoscape.py ../../data/nbh/graph.edges ../../data/outputGraph.txt 8 ../../data/nbh/graph.label

Cytoscape &
```

Puis File > Import > Import Network from file



## Graphe foward et reverse


