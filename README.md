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

```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ../results/ --genomeFastaFiles ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna --sjdbGTFfile ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf -sjdbOverhang 74 --genomeSAindexNbases 13 --genomeSAsparseD 4
```

### Etape 3: Génération du graphe de De Bruijn

```
./build/bin/kissplice -r read.fa 
```
Attention, il s'agit là d'une solution over-killed, il vaudra mieux repartir de `BCaml` ! (Et attention, être dans le bon répertoire, pour moi data)

### Etape 4: Calcul des poids des sommets

On calcule les poids de chaque sommet (rajouter ici une explication)

```
g++ graph.cpp -o graph.exe
cd ../../../kissplice_moustique
../DmGoth/stage-M2/scr/graph.exe graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges 200 -o ../DmGoth/data/outputGraph.txt 
cd ../DmGoth/stage-M2/scr
```

Puis on lance l'analyse graphique (dot) et on vérifie que la valeur donnée pour le threshold semble correcte (top10)

```
python3 plot.py ../../data/outputGraph.txt top10
python3 plot.py ../../data/outputGraph.txt dot
```

Finallement, on génère un fichier de reads à aligner sur le génome de référence. (Ici threshold = 11)

```
python3 reads_to_align.py ../../data/outputGraph.txt ../../data/read.fq 8
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


### Etape 8: Afficher une séquence sur Cytoscape:

```
python3 labellingCytoscape.py ../../data/nbh.edges ../../data/outputGraph.txt ../../results/seq.txt ../../data/nbh.label

/data/home/vincent/TiffanyDelhomme$ ./nbh -n /localdata/pandata/students/Projet_KS/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes -e /localdata/pandata/students/Projet_KS/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges -k 41 -o /localdata/pandata/students/Projet_KS/DmGoth/data/nbh -d 20 -q CCTCCCCCCCTTGGAGCGTGACGTAATTTGTGCATGACCCC
Cytoscape &
```

Puis File > Import > Import Network from file