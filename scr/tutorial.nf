#!/usr/bin/env nextflow


println "\tDébut de la Pipeline Nextflow\nA executer dans le dossier scr du serveur pedago-ngs.\n"

moust = []
moust.name = "moustique"
moust.genome = "../../data/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna"
moust.gtf = "../../data/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf"
moust.nodes = "../../../kissplice_results/kissplice_${name}/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes"
moust.edges = "../../../kissplice_results/kissplice_${name}/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges"


donnees = Channel.from(moust)

process creaCarte {
    input:
    env spe from donnes

    output:
    path '../../results/${spe.name}/STAR' into STARDir

    script:
    """
    mkdir ../../results/${spe.name}
    STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ../../results/${spe.name} \
    --genomeFastaFiles ${spe.genome} \
    --sjdbGTFfile ${spe.gtf} \
    -sjdbOverhang 74 \
    --genomeSAindexNbases 13 \
    --genomeSAsparseD 4
    """
}


process calculpoids {
    input:
    env spe from donnes
    path STAR from STARDir

    output:
    path '../../DmGoth/data/outputGraph${spe.name}.txt' into OutputGraph

    script:
    """
    g++ graph.cpp main.cpp -o graph.exe
    graph.exe  10 -o ../../DmGoth/data/outputGraph${spe.name}.txt
    """
}