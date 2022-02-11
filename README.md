# stage-M2



## Utiles


sshfs sdarmon@pedago-ngs:/localdata/pandata/students/Projet_KS/DmGoth peda/

scp -p ../../Bureau/genome_assemblies_genome_fasta.tar sdarmon@pedago-ngs:/localdata/pandata/students/Projet_KS/DmGoth/data



## proto



download ref_genome.fna et genome.gtf et read.fa
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ../results/ --genomeFastaFiles ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna --sjdbGTFfile ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf -sjdbOverhang 74 --genomeSAindexNbases 13 --genomeSAsparseD 4

./build/bin/kissplice -r read.fa (attention, être dans le bon répertoire)

g++ graph.cpp -o graph.exe

./graph.exe file.nodes file.edges -o output.txt

python3 plot.py output.txt dot
python3 plot.py output.txt top10
python3 reads_to_align.py input.txt output.fq threshold where the input is the output of graph.exe, output.fq est le format des séquences à générer et threshold est l'output de plot.py top10.