#!/bin/bash
#
##Echantillonage of the reads
#Done

##Compute the HC
echo "HC of the reads ..."
python3 homomorphic_compression.py  /beegfs/data/sdarmon/mus/mus_10_1.fq /beegfs/data/sdarmon/mus/hc_mus_1.fq 5
python3 homomorphic_compression.py  /beegfs/data/sdarmon/mus/mus_10_2.fq /beegfs/data/sdarmon/mus/hc_mus_2.fq 5

##Compute the DGB with kissplice
echo "DGB with kissplice ..."
kissplice-binary-ubuntu-2.6.7/bin/kissplice -r /beegfs/data/sdarmon/mus/hc_mus_1.fq -r /beegfs/data/sdarmon/mus/hc_mus_2.fq -k 41

##Compute the weighting
echo "Weighting of the nodes..."
#g++ -g graph.cpp ponderation.cpp -o graph.exe
./graph.exe results/graph_hc_mus_1_hc_mus_2_k41.nodes results/graph_hc_mus_1_hc_mus_2_k41_C0.05.edges 10 -k 41 -o results/outputNodes.txt


## Compute the connexe components
echo "Agglomeration of connexe components..."
mkdir -p /beegfs/data/sdarmon/results/
mkdir -p /beegfs/data/sdarmon/results/mus/
mkdir -p /beegfs/data/sdarmon/results/mus/processing/

g++ -g graph.cpp agglo.cpp -o agglo.exe
./agglo.exe results/outputNodes.txt results/graph_hc_mus_1_hc_mus_2_k41_C0.05.edges -c 5 -d 10  /beegfs/data/sdarmon/results/mus -clean results/graph_hc_mus_1_hc_mus_2_k41.abundance > /beegfs/data/sdarmon/results/mus/rapportAgglo.txt


## Compute the genome STAR ref

mkdir -p /beegfs/data/sdarmon/results/mus
mkdir -p /beegfs/data/sdarmon/results/mus/te_ref
mkdir -p /beegfs/data/sdarmon/results/mus/ref

/beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 8     --runMode genomeGenerate     --genomeDir /beegfs/data/sdarmon/results/mus/ref     --genomeFastaFiles  /beegfs/data/sdarmon/mus/Mus_musculus.GRCm39.dna.primary_assembly.fa --genomeSAindexNbases 13  --genomeSAsparseD 4 --genomeChrBinNbits 16

/beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 8     --runMode genomeGenerate     --genomeDir /beegfs/data/sdarmon/results/mus/te_ref     --genomeFastaFiles  /beegfs/data/sdarmon/mus/TE.fa --genomeSAindexNbases 13  --genomeSAsparseD 4 --genomeChrBinNbits 16


## Compute the intersection
#bash inter_mus.sh


##Fetch the neighborhood of a comp 
#bash neigh_comp.sh ../../results/dog/genome/ ../../results/dog/te_ref/ ../../data/dog/ref.gtf ../../data/dog/TE.gtf ../../data/dog/outputNodes.txt ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41_C0.05.edges ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41.abundance ../../results/dog/processing/comp2.txt ~/Documents/graph/sineCA2 41 10 10


##Build the vizitig file for visualisation
#./visitig ...
