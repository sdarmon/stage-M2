#!/usr/bin/env nextflow

workDir = '/home/sdarmon/Documents/stage-M2/peda/DmGoth/stage-M2/scr'

println "\tDébut de la Pipeline Nextflow\nA executer dans le dossier scr du serveur pedago-ngs.\nExecution type: \n\t-nextflow run tutorial.nf --path [YourWorkDirPath]\n\t-pwd | xargs nextflow run tutorial.nf --path : Permet de récupérer directement le WorkDir dedans\n "

if (params.path != null){
workDir = params.path
println "Path loaded\n"
}


moust = ["name":"", "genome":"", "gtf":"", "nodes":"", "edges":"", "TE":""]
moust["name"] = "moustique"
moust["genome"] = "${workDir}/../../data/moustique/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna"
moust["gtf"] = "${workDir}/../../data/moustique/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf"
moust["nodes"] = "${workDir}/../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes"
moust["edges"] = "${workDir}/../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges"
moust["TE"] = "${workDir}/../../data/moustique/AaegL5_TE_repeats.gff"

chien = ["name":"", "genome":"", "gtf":"", "nodes":"", "edges":"", "TE":""]
chien["name"] = "chien"
chien["genome"] = "${workDir}/../../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-unmasked.fa"
chien["gtf"] = "${workDir}/../../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-2021_03-genes.sorted.gtf"
chien["nodes"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes"
chien["edges"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges"
chien["TE"] = "${workDir}/../../data/chien/Cfam_GSD_TE.gtf"


topVal = Channel.from("top10","top10")
topAgglo = Channel.from("top1","top1")

donnees = Channel.from(moust) //moust,chien
// = Channel.from() //moust
//intersecter = Channel.from()  //moust
//agglo = Channel.from() //moust


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
    mkdir -p ${workDir}/../../results/${name}
    mkdir -p ${workDir}/../../results/${name}/genome
    STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ${workDir}/../../results/${name}/genome \
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
    val spe into aligner

    script:
    name = spe.name
    nodes = spe.nodes
    edges = spe.edges
    """
    g++ -g ${workDir}/graph.cpp ${workDir}/main.cpp -o ${workDir}/graph.exe
    ${workDir}/graph.exe  ${nodes} ${edges} 10 -o ${workDir}/../../data/${name}/outputGraph${spe.name}.txt
    python3 ${workDir}/reads_to_align.py ${workDir}/../../data/${name}/outputGraph${spe.name}.txt ${workDir}/../../data/${name}/outputGraph${spe.name}Clean.txt 0 -clean
    """
}



process top {
    input:
    val spe from aligner
    val topV from topVal

    output:
    stdout into top
    val spe into aligner2

    exec:
    name = spe.name

    script:
    """
    python3 ${workDir}/plot.py ${workDir}/../../data/${name}/outputGraph${spe.name}Clean.txt ${topV}
    """
}

process read_to_align {

    input:
    val spe from aligner2
    val value from top

    output:
    val spe into intersecter

    exec:
    name = spe.name
    println value

    script:
    """
    python3 ${workDir}/reads_to_align.py ${workDir}/../../data/${name}/outputGraph${spe.name}Clean.txt ${workDir}/../../data/${name}/read${name}.fq ${value}
    mkdir -p ${workDir}/../../results/${name}/STAR_alignment
    STAR --genomeDir ${workDir}/../../results/${name}/genome \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${workDir}/../../data/${name}/read${name}.fq \
    --outFileNamePrefix ${workDir}/../../results/${name}/STAR_alignment/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --outReadsUnmapped Fastx
    """
}

process intersect {
    input:
    val spe from intersecter

    output:
    val spe into agglo

    exec:
    name = spe.name
    TE = spe.TE

    script:
    """
    bedtools intersect -wa -a ${TE} -b ${workDir}/../../results/${name}/STAR_alignment/Aligned.sortedByCoord.out.bam > ${workDir}/../../results/${name}/intersectionTE.txt
    bedtools intersect -wb -a ${TE} -b ${workDir}/../../results/${name}/STAR_alignment/Aligned.sortedByCoord.out.bam > ${workDir}/../../results/${name}/intersectionKiss.txt
    python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/intersectionKiss.txt ${workDir}/../../results/${name}/intersectionKissNoDouble.txt -s 12
    python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/intersectionTE.txt ${workDir}/../../results/${name}/intersectionTENoDouble.txt -t 8
    echo "Intersections uniques dans KisSplice : " > ${workDir}/../../results/${name}/rapportIntersect.txt
    wc -l ${workDir}/../../results/${name}/intersectionKissNoDouble.txt >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\nIntersections uniques dans les TE : " >> ${workDir}/../../results/${name}/rapportIntersect.txt
    wc -l ${workDir}/../../results/${name}/intersectionTENoDouble.txt >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\n" >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep "Number of input reads" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\n" >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep "Uniquely mapped reads number" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\n" >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep "Number of reads unmapped: too many mismatches" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\n" >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep "Number of reads unmapped: too short" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\n" >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep "Number of reads unmapped: other" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    """
}

process topA {
    input:
    val spe from agglo
    val topV from topAgglo

    output:
    stdout into topA
    val spe into agglo2

    exec:
    name = spe.name

    script:
    """
    python3 ${workDir}/plot.py ${workDir}/../../data/${name}/outputGraph${name}Clean.txt ${topV}
    """
}


process agglomeration {
    input:
    val value from topA
    val spe from agglo2

    output:
    val spe into agglomerate
    exec:
    name = spe.name
    edges = spe.edges

    script:
    """
    mkdir -p ${workDir}/../../results/${name}/processing
    g++ -g ${workDir}/graph.cpp ${workDir}/agglo.cpp -o ${workDir}/agglo.exe
    ${workDir}/agglo.exe ${workDir}/../../data/${name}/outputGraph${name}Clean.txt \
    ${edges} \
    -c ${value.replaceAll(/\n/, "")} \
    -d 100\
    "${workDir}/../../results/${name}" \
    > ${workDir}/../../results/${name}/rapportAgglo.txt
    """
}

process intersectComp {

    input:
    val spe from agglomerate

    exec:
    name = spe.name
    edges = spe.edges
    TE = spe.TE

    script:
    """
        for i in {0..99..1}
        do
            python3 \
            ${workDir}/reads_to_align.py ${workDir}/../../results/${name}/processing/comp\$i.txt \
            ${workDir}/../../results/${name}/processing/comp\$i.fq \
            0
            mkdir -p ${workDir}/../../results/${name}/processing/STAR_alignment
            STAR --genomeDir ${workDir}/../../results/${name}/genome \
            --runMode alignReads \
            --runThreadN 8 \
            --readFilesIn ${workDir}/../../results/${name}/processing/comp\$i.fq  \
            --outFileNamePrefix ${workDir}/../../results/${name}/processing/STAR_alignment/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --outFilterMultimapNmax 10000 \
            --outReadsUnmapped Fastx
            bedtools intersect -wa -a ${TE} \
            -b ${workDir}/../../results/${name}/processing/STAR_alignment/Aligned.sortedByCoord.out.bam \
            > ${workDir}/../../results/${name}/processing/intersectionTE\$i.txt
        done
    """
}
