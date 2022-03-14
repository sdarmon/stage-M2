#!/usr/bin/env nextflow


println "\tDébut de la Pipeline Nextflow\nA executer dans le dossier scr du serveur pedago-ngs.\n"
workDir = '/home/sdarmon/Documents/stage-M2/peda/DmGoth/stage-M2/scr'
moust = ["name":"", "genome":"", "gtf":"", "nodes":"", "edges":""]
moust["name"] = "moustique"
moust["genome"] = "data/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna"
moust["gtf"] = "/home/sdarmon/Documents/stage-M2/peda/DmGoth/data/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf"
moust["nodes"] = "/home/sdarmon/Documents/stage-M2/peda/kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes"
moust["edges"] = "/home/sdarmon/Documents/stage-M2/peda/kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges"


donnees = Channel.from(moust)

process creaCarte {
    input:
    val spe from donnees

    output:
    val spe into STARDir

    script:
    name = spe.name
    genome = spe.genome
    gtf = spe.gtf

    """
    mkdir -p ../../results/${name}
    echo '$PATH'
    STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ../../results/${name} \
    --genomeFastaFiles ${genome} \
    --sjdbGTFfile ${gtf} \
    -sjdbOverhang 74 \
    --genomeSAindexNbases 13 \
    --genomeSAsparseD 4
    """
}


process calculpoids {
    input:
    val spe from STARDir

    output:
    val spe into OutputGraph

    script:
    """
    g++ graph.cpp main.cpp -o graph.exe
    graph.exe  10 -o ../../DmGoth/data/outputGraph${spe.name}.txt
    """
}