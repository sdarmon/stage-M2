#!/usr/bin/env nextflow

workDir = '/home/sdarmon/Documents/stage-M2/peda/DmGoth/stage-M2/scr'

println "\tDébut de la Pipeline Nextflow\nA executer dans le dossier scr du serveur pedago-ngs.\nExecution type: \n\t-nextflow run tutorial.nf --path [YourWorkDirPath]\n\t-pwd | xargs nextflow run tutorial.nf --path : Permet de récupérer directement le WorkDir dedans\n "

if (params.path != null){
workDir = params.path
println "Path loaded\n"
}


moust = ["name":"", "genome":"", "gtf":"", "nodes":"","abundance":"",  "edges":"", "TE":""]
moust["name"] = "moustique"
moust["genome"] = "${workDir}/../../data/moustique/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.fna"
moust["gtf"] = "${workDir}/../../data/moustique/ncbi-genomes-2022-02-11/GCF_002204515.2_AaegL5.0_genomic.gtf"
moust["nodes"] = "${workDir}/../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.nodes"
moust["abundance"] = "${workDir}/../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41.abundance"
moust["edges"] = "${workDir}/../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges"
moust["TE"] = "${workDir}/../../data/moustique/AaegL5_TE_repeats.gff"

chien = ["name":"", "genome":"", "gtf":"", "nodes":"","abundance":"",  "edges":"", "TE":""]
chien["name"] = "chien"
chien["genome"] = "${workDir}/../../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-unmasked.fa"
chien["gtf"] = "${workDir}/../../../fastq/SRR_Chiens/Canis_lupus_familiaris-GCA_011100685.1-2021_03-genes.sorted.gtf"
chien["nodes"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes"
chien["abundance"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.abundance"
chien["edges"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges"
chien["TE"] = "${workDir}/../../data/chien/Cfam_GSD_TE.gtf"

chien2 = ["name":"", "genome":"", "gtf":"", "nodes":"","abundance":"",  "edges":"", "TE":""]
chien2["name"] = "chien2"
chien2["genome"] = "${workDir}/../../data/chien2/canFam4.fa"
chien2["gtf"] = "${workDir}/../../data/chien2/refGene.gtf"
chien2["nodes"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes"
chien2["abundance"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.abundance"
chien2["edges"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges"
chien2["TE"] = "${workDir}/../../data/chien2/Cfam_GSD_TE_Americain.gtf"


human = ["name":"", "genome":"", "gtf":"", "nodes":"","abundance":"", "edges":"", "TE":""]
human["name"] = "human"
human["genome"] = "${workDir}/../../../fastq/Human/ref_genome_Homo.fa"
human["gtf"] = "${workDir}/../../../fastq/Human/hg38.ncbiRefSeq.gtf"
human["nodes"] = "${workDir}/../../../kissplice_results/kissplice_human_fastp/graph_fastp_SknshRACellRep1_10M_fastp_SknshRACellRep2_10M_fastp_SknshCellRep3_10M_fastp_SknshCellRep4_10M_k41.nodes"
human["abundance"] = "${workDir}/../../../kissplice_results/kissplice_human_fastp/graph_fastp_SknshRACellRep1_10M_fastp_SknshRACellRep2_10M_fastp_SknshCellRep3_10M_fastp_SknshCellRep4_10M_k41.abundance"
human["edges"] = "${workDir}/../../../kissplice_results/kissplice_human_fastp/graph_fastp_SknshRACellRep1_10M_fastp_SknshRACellRep2_10M_fastp_SknshCellRep3_10M_fastp_SknshCellRep4_10M_k41_C0.05.edges"
human["TE"] = "${workDir}/../../fastq/Human/Human_TE.gtf"


topVal = Channel.from("top10","top10")
topAgglo = Channel.from("top1","top1")

donnees = Channel.from(human,chien) //moust,chien
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
    bedtools intersect  -wa -a ${TE} -b ${workDir}/../../results/${name}/STAR_alignment/Aligned.sortedByCoord.out.bam > ${workDir}/../../results/${name}/intersectionTE.txt
    bedtools intersect  -wb -a ${TE} -b ${workDir}/../../results/${name}/STAR_alignment/Aligned.sortedByCoord.out.bam > ${workDir}/../../results/${name}/intersectionKiss.txt
    python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/intersectionKiss.txt ${workDir}/../../results/${name}/intersectionKissNoDouble.txt -s 12
    python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/intersectionTE.txt ${workDir}/../../results/${name}/intersectionTENoDouble.txt -t 8
    echo "Intersections uniques dans KisSplice : " > ${workDir}/../../results/${name}/rapportIntersect.txt
    wc -l ${workDir}/../../results/${name}/intersectionKissNoDouble.txt >> ${workDir}/../../results/${name}/rapportIntersect.txt
    echo "\nIntersections uniques dans les TE : " >> ${workDir}/../../results/${name}/rapportIntersect.txt
    wc -l ${workDir}/../../results/${name}/intersectionTENoDouble.txt >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep -F "Number of input reads" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    grep -F "Uniquely mapped reads number" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    #grep -F "Number of reads unmapped: too many mismatches" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    #grep -F "Number of reads unmapped: too short" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
    #grep -F "Number of reads unmapped: other" ${workDir}/../../results/${name}/STAR_alignment/Log.final.out >> ${workDir}/../../results/${name}/rapportIntersect.txt
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
    ab = spe.abundance
    script:
    """
    mkdir -p ${workDir}/../../results/${name}/processing
    g++ -g ${workDir}/graph.cpp ${workDir}/agglo.cpp -o ${workDir}/agglo.exe
    ${workDir}/agglo.exe ${workDir}/../../data/${name}/outputGraph${name}Clean.txt \
    ${edges} \
    -c ${value.replaceAll(/\n/, "")} \
    -d 10 \
    "${workDir}/../../results/${name}" \
    -clean ${ab}\
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
        python3 rapportAgglo.py ${TE} ${workDir}/../../results/${name}/processing/intersectionTE 100 -target > ${workDir}/../../results/${name}/rapportHisto.txt
    """
}

"""


