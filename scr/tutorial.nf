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

chien2 = ["name":"", "genome":"", "gtf":"", "nodes":"", "edges":"", "TE":""]
chien2["name"] = "chien2"
chien2["genome"] = "${workDir}/../../data/chien2/canFam4.fa"
chien2["gtf"] = "${workDir}/../../data/chien2/refGene.gtf"
chien2["nodes"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41.nodes"
chien2["edges"] = "${workDir}/../../../kissplice_results/kissplice_Chiens/graph_SRR15254976_1_SRR15254976_2_SRR15254978_1_SRR15254978_2_SRR15254980_1_SRR15254980_2_SRR15254982_1_SRR15254982_2_SRR15254985_1_SRR15254985_2_SRR15254986_1_SRR15254986_2_SRR15254989_1_SRR15254989_2_SRR1k41_C0.05.edges"
chien2["TE"] = "${workDir}/../../data/chien2/Cfam_GSD_TE_Americain.gtf"


topVal = Channel.from("top10","top10")
topAgglo = Channel.from("top1","top1")

donnees = Channel.from(chien,moust) //moust,chien
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
    script:
    """
    mkdir -p ${workDir}/../../results/${name}/processing
    g++ -g ${workDir}/graph.cpp ${workDir}/agglo.cpp -o ${workDir}/agglo.exe
    ${workDir}/agglo.exe ${workDir}/../../data/${name}/outputGraph${name}Clean.txt \
    ${edges} \
    -c ${value.replaceAll(/\n/, "")} \
    -d 10 \
    "${workDir}/../../results/${name}" \
    -clean \
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
        grep "SINEC2A1_CF" ${gtf} > ${workDir}/../../data/${name}/SINEC2A1_CF.gtf
        for i in {0..99..1}
        do
        STAR --genomeDir ${workDir}/../../results/${name}/genome \
                      --runMode alignReads \
                      --runThreadN 8 \
                      --readFilesIn ${workDir}/../../results/${name}/processing/comp$i.fq  \
                      --outFileNamePrefix ${workDir}/../../results/${name}/processing/STAR_alignment/ \
                      --outSAMtype BAM SortedByCoordinate \
                      --outSAMunmapped Within \
                      --outSAMattributes Standard \
                      --outFilterMultimapNmax 10000 \
                      --outReadsUnmapped Fastx
            bedtools intersect -wb -a ${workDir}/../../data/${name}/SINEC2A1_CF.gtf \
            -b ${workDir}/../../results/${name}/processing/STAR_alignment/Aligned.sortedByCoord.out.bam \
            > ${workDir}/../../results/${name}/processing/intersectionSINEC2A1_CF$i.txt
        python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/processing/intersectionSINEC2A1_CF$i.txt \
            ${workDir}/../../results/${name}/processing/intersectionSINEC2A1_CFNoDouble$i.txt -s 12
        done
        for i in {0..99..1}
                do
                echo ${workDir}/../../results/${name}/processing/intersectionSINEC2A1_CFNoDouble$i.txt >> ${workDir}/../../results/${name}/stackSINEC2A1_CF.txt
                done
        python3 recupSeq.py ${workDir}/../../results/${name}/processing/intersectionSINEC2A1_CFNoDouble \
            ${workDir}/../../results/${name}/processing/comp \
            100 \
            12 \
            > ${workDir}/../../results/${name}/seqSINEC2A1_CF.txt

        python3 \
            ${workDir}/reads_to_align.py ${workDir}/../../results/${name}/seqSINEC2A1_CF.txt \
            ${workDir}/../../results/${name}/seqSINEC2A1_CF.fq \
            0
        mkdir -p ${workDir}/../../results/${name}/processing/STAR_alignment
        STAR --genomeDir ${workDir}/../../results/${name}/genome \
            --runMode alignReads \
            --runThreadN 8 \
            --readFilesIn ${workDir}/../../results/${name}/seqSINEC2A1_CF.fq  \
            --outFileNamePrefix ${workDir}/../../results/${name}/processing/STAR_alignment/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --outFilterMultimapNmax 10000 \
            --outReadsUnmapped Fastx
        bedtools intersect -wa -a ${TE} \
            -b ${workDir}/../../results/${name}/processing/STAR_alignment/Aligned.sortedByCoord.out.bam \
            > ${workDir}/../../results/${name}/intersectionTECompSeqSINEC2A1_CF.txt
        python3 ${workDir}/suppDoublon.py ${workDir}/../../results/${name}/intersectionTECompSeqSINEC2A1_CF.txt \
                    ${workDir}/../../results/${name}/intersectionTECompSeqSINEC2A1_CFNoDouble.txt -t 8
    """

#['"L2"', '"L1MDa"', '"HAL1ME"', '"L1M4"', '"L2d"', '"L1MC4"', '"BLACKJACK"', '"MIR"', '"Tigger1"', '"MamRTE1"', '"MLT1C"', '"L1ME1"', '"L1MEc"', '"SINEC2A1_CF"', '"ERVL-B4-int"', '"SINEC_c2"', '"MIRc"', '"L1_Carn3"', '"MER5C"', '"L1MEf"', '"L1_Canis1"', '"L2a"', '"L1_Cf"', '"Tigger2a_Car"', '"MamGypLTR1a"', '"SINEC_a2"', '"MLT1B"', '"L1_Carn5"', '"MER21C"', '"L1ME2"', '"L1MC"', '"MLT2D"', '"L1ME2z"', '"L1ME3"', '"L1MC1"', '"L1MB2"', '"L1MB3"', '"L1M5"', '"SINEC_a1"', '"L1MC3"', '"MLT2C2"', '"HAL1b"', '"SINEC_b2"', '"L1MD1"', '"L1MB5"', '"MLT1J"', '"L1MC2"', '"LTR78B"', '"L1MA9"', '"MLT1J2"', '"L1_Carn2"', '"L1M8"', '"L1_Carn7"', '"SINEC_Cf2"', '"L1_Carni"', '"MER21-int"', '"L1ME3Cz"', '"Charlie8"', '"MER103C"', '"LTR16A"', '"L1MD"', '"Tigger14a"', '"L1MC4a"', '"L1M3c"', '"MLT1E1"', '"SINEC1D_CF"', '"L4_C_Mam"', '"L4_B_Mam"', '"L1_Canid_"', '"MLT1A0"', '"Tigger13a"', '"L1_Canis2"', '"Zaphod"', '"Tigger18a"', '"L1MEi"', '"Tigger17a"', '"SINEC_Cf"', '"L1MB8"', '"L3"', '"SINEC_c1"', '"SINEC_Cf3"', '"ERV3-16A3_I-int"', '"Charlie21a"', '"MamGypLTR2c"', '"L2d2"', '"SINEC_old"', '"MLT1C2"', '"L2c"', '"MER94"', '"Charlie15b"', '"Charlie14a"', '"L1MEd"', '"HAL1"', '"L1MCc"', '"L1MD2"', '"LTR40a"', '"MIRb"', '"MER89"', '"L1MEj"', '"L1M6"', '"Charlie1a"', '"CarLTR11a"', '"CarLTR1A2A"', '"L1MB1"', '"L2b"', '"MLT1H1"', '"L1ME4a"', '"L1MA10"', '"MER34B-int"', '"L1M4c"', '"CarLTR8b"', '"GA-rich"', '"MLT1D"', '"MamSINE1"', '"LTR37A"', '"MER3"', '"X5a_DNA"', '"MER102a"', '"MER33"', '"MLT1K"', '"L1MC5"', '"LTR84a"', '"L1ME3A"', '"MamGypLTR1c"', '"MIR3"', '"L1MB7"', '"SINEC_b1"', '"L1M4a1"', '"CarERVR1_LTR"', '"MER20"', '"tRNA-Lys-AAA"', '"L1MB4"', '"L1MC5a"', '"L1M3b"', '"CarLTR2"', '"LTR33A_"', '"CfERV1z_LTR"', '"CfERV1a_LTR"', '"SINEC1C1_CF"', '"A-rich"', '"L1M3"', '"G-rich"', '"L1_Carn1"', '"SINEC1C2_CF"', '"SINEC2A2_CF"', '"Kanga1b"', '"Kanga1d"', '"SINEC1A_CF"', '"MER63D"', '"Zaphod4a"', '"Tigger2f"', '"MLT1H-int"', '"CfERV1-int"', '"MLT1H"', '"Charlie15a"', '"EuthAT-2"', '"Tigger6a"', '"L1MEg"', '"L1ME3G"', '"MLT1G"', '"LTR33A"', '"CarERV2a_LTR"', '"MER5A1"', '"MLT1A"', '"MER5A"', '"MARNA"', '"CfERV1c_LTR"', '"LTR107_Mam"', '"MER34A_CF"', '"MER34-int"', '"MER31-int"', '"Charlie30b"', '"CarLTR9"', '"L3b"', '"CarLTR10"', '"CarLTR4"', '"MER34A"', '"MER53"', '"Tigger15a"', '"LTR103b_Mam"', '"MamGypLTR2b"', '"MER58A"', '"X6A_LINE"', '"MamRep137"', '"MER91A"', '"Charlie26a"', '"tRNA-Lys-AAG"', '"L1M4a2"', '"Charlie2a"', '"Charlie2b"', '"Tigger17b"', '"Tigger17d"', '"MLT2B1"', '"MIR1_Amn"', '"Cheshire"', '"Kanga1"', '"MER45C"', '"MLT1N2"', '"MLT1G3"', '"L1ME4b"', '"MER89-int"', '"MER21B"', '"MER46C"', '"MER74A"', '"LTR51_CF"', '"L1ME3E"', '"hAT-5_Mam"', '"MLT1E"', '"L1ME4c"', '"L1M4b"', '"Charlie13a"', '"LTR16C"', '"OldhAT1"', '"Plat_L3"', '"L1ME3B"', '"MER113A"', '"MamTip2b"', '"AmnL2-1"', '"MER63B"', '"MER112"', '"MER115"', '"MamRep434"', '"MER121"', '"MER77"', '"SINEC1B1_CF"', '"MLT1-int"', '"L1MCa"', '"MLT1I"', '"Charlie10"', '"MER68-int"', '"MER68"', '"LTR67B"', '"MER58C"', '"MLT1J1"', '"MLT1F2"', '"Tigger9a"', '"SINEC1B2_CF"', '"L1MDb"', '"LTR37B"', '"Penelope1_Vert"', '"Charlie9"', '"UCON105"', '"X6b_DNA"', '"SINEC1_CF"', '"Charlie4a"', '"ORSL-2b"', '"X17_LINE"', '"X30_DNA"', '"MER58B"', '"MLT2B2"', '"Charlie17b"', '"L1M7"', '"MamRep605"', '"LTR78"', '"LTR81C"', '"Tigger10"', '"MLT2C1"', '"Tigger19b"', '"CfERVF2_LTR"', '"L4_A_Mam"', '"MER131"', '"AmnSINE1"', '"Eutr2"', '"MER96B"', '"MLT1F-int"', '"LTR104_Mam"', '"MER5B"', '"LTR10_EC"', '"L1MD3"', '"Kanga1a"', '"LTR40b"', '"MamGypLTR3"', '"Kanga11a"', '"Tigger6"', '"X9a_DNA"', '"MER77B"', '"MLT1L"', '"Tigger17"', '"DNA1_Mam"', '"UCON66"', '"HAL1M8"', '"CarLTR1A2B"', '"MER119"', '"CR1_Mam"', '"LFSINE_Vert"', '"MADE2"', '"U5"', '"MamTip2"', '"MLT2B4"', '"Charlie23a"', '"MLT1E2"', '"LTR41"', '"LTR91A"', '"Charlie17a"', '"LTR79"', '"MLT1E1A"', '"Charlie20a"', '"CarERV2c1_LTR"', '"Charlie1b"', '"Kanga1c"', '"Charlie4"', '"CarERV2b1_LTR"', '"MER68B"', '"LTR82B"', '"MamGypLTR1b"', '"L1ME3C"', '"HERV16-int"', '"MamGypLTR3a"', '"L1MCb"', '"MER67C"', '"Arthur1A"', '"MER92-int"', '"MER92D"', '"Tigger22N1"', '"Charlie10a"', '"Tigger12A"', '"MER102c"', '"Helitron2Na_Mam"', '"X26_DNA"', '"L1ME3D"', '"U6"', '"LTR16E2"', '"Charlie7"', '"LTR90A"', '"MLT1F"', '"LTR91"', '"LTR53"', '"L1MEh"', '"L1ME5"', '"hAT-N1a_Mam"', '"Tigger16a"', '"X12_DNA"', '"MER135"', '"MARE4"', '"Tigger1a_Car"', '"Tigger19a"', '"UCON41"', '"LTR37-int"', '"MLT2F"', '"MamGypsy2-I"', '"MLT1H2"', '"CR1-3_Croc"', '"MER41_CF"', '"Eutr18"', '"MER58D"', '"MER45A"', '"MamRep1527"', '"MLT1O"', '"LTR68"', '"CfERV1b_LTR"', '"Charlie18a"', '"X7A_LINE"', '"UCON26"', '"X7B_LINE"', '"MER91B"', '"LTR16A2"', '"Tigger16b"', '"MamGypsy2-LTR"', '"MamRep605b"', '"X4b_DNA"', '"LTR85a"', '"UCON4"', '"UCON40"', '"L1MEb"', '"AmnSINE2"', '"MER91C"', '"EUTREP14"', '"Arthur1B"', '"Arthur1"', '"X6a_DNA"', '"MLT1F1"', '"LTR33C"', '"Charlie5"', '"Tigger12"', '"MER20B"', '"MER126"', '"LTR16"', '"Charlie22a"', '"MLT1G1"', '"MLT1M"', '"LTR109A1"', '"HERVL40-int"', '"LTR106_Mam"', '"Charlie19a"', '"ERV54-EC_LTR"', '"CarERV4-int"', '"Charlie4z"', '"LTR16E1"', '"Charlie16"', '"Charlie16a"', '"LTR89B"', '"LTR88a"', '"LTR40c"', '"UCON48"', '"MER94B"', '"MER34C_CF"', '"CarERV4a_LTR"', '"MamRep38"', '"LTR105_Mam"', '"LTR89"', '"Arthur2"', '"MER117"', '"LTR72A_CF"', '"MLT1J2-int"', '"L1M3a"', '"MER97b"', '"MLT1E3"', '"MER90a"', '"CarLTR11b"', '"LTR33"', '"Charlie29a"', '"Eutr15"', '"MER34A1"', '"MER90"', '"MamRep1894"', '"L1M"', '"MER110A"', '"MER110"', '"CarLTR1A1A"', '"LTR85c"', '"MER34B_CF"', '"LTR16B1"', '"LTR83"', '"X1_DNA"', '"L1MEg1"', '"MLT1J-int"', '"LTR16B2"', '"CarLTR1B1"', '"Kanga2_a"', '"MER67D"', '"X24_DNA"', '"ERVL-E-int"', '"Eutr10"', '"LTR50"', '"LTR80B"', '"MLT1D-int"', '"MLT1G1-int"', '"Charlie1"', '"LTR33B"', '"LTR86C"', '"CarLTR7"', '"MamRep4096"', '"Tigger17c"', '"MER31B"', '"MER57LA"', '"MER63C"', '"MER63A"', '"LTR85b"', '"X24_LINE"', '"LTR55"', '"LTR88b"', '"LTR16A1"', '"Zaphod3"', '"LTR16D"', '"MLT1H1-int"', '"X2_LINE"', '"LTR16B"', '"LTR52"', '"LTR81"', '"tRNA-Gly-GGA"', '"CarLTR3"', '"CarERV2c2_LTR"', '"MLT1A0-int"', '"L1_Canid2"', '"5S"', '"UCON1"', '"MER95"', '"LTR109A2"', '"X9b_DNA"', '"MER110-int"', '"CarLTR6c"', '"Zaphod2"', '"LTR86B1"', '"MLT1F2-int"', '"LTR102_Mam"', '"LTR82A"', '"U2"', '"CarLTR8y"', '"Charlie7a"', '"Charlie25"', '"Charlie24"', '"Charlie13b"', '"Tigger9b"', '"L1ME3F"', '"LTR87"', '"LTR75"', '"MER106B"', '"MER81"', '"CR1-L3A_Croc"', '"MamRep1879"', '"SINEC_Fc2"', '"U4"', '"MER124"', '"MER104"', '"LTR81A"', '"X9_LINE"', '"MLT1I-int"', '"MLT1E2-int"', '"LTR84b"', '"CarERV4b_LTR"', '"MER74B"', '"ERVL-int"', '"CarERV2b2_LTR"', '"Charlie29b"', '"LTR41C"', '"L1MEg2"', '"MLT1H2-int"', '"MER92A"', '"MER106A"', '"MamRep1151"', '"MER76-int"', '"L1M2"', '"UCON31"', '"MamGypLTR4"', '"MERX"', '"7SK"', '"MamRep564"', '"LTR101_Mam"', '"MER45B"', '"MER92C"', '"CarLTR1B2"', '"LTR108e_Mam"', '"X15_LINE"', '"CarERVR1-int"', '"CarLTR1-int"', '"MER102b"', '"MamTip3"', '"MamTip1"', '"Tigger20a"', '"LTR41B"', '"Tigger2"', '"MER45R"', '"Looper"', '"ERV3-16A3_LTR"', '"MER31A"', '"MER92B"', '"CarLTR6b"', '"SINEC_Fc"', '"MER127"', '"Mam_R4"', '"Charlie31a"', '"MER96"', '"Tigger6b"', '"EUTREP7"', '"X11_DNA"', '"Helitron1Na_Mam"', '"MER105"', '"LTR88c"', '"UCON34"', '"CarLTR5"', '"MER113"', '"X8_LINE"', '"CarERV2-int"', '"MER99"', '"ORSL"', '"L2-2_Mam"', '"MLT1A-int"', '"X20_LINE"', '"Tigger12c"', '"MER5C1"', '"Charlie30a"', '"Helitron3Na_Mam"', '"UCON6"', '"MLT2E"', '"UCON44"', '"ORSL-2a"', '"hAT-N1_Mam"', '"Charlie10b"', '"MamGypLTR1d"', '"LTR9A_EC"', '"LSU-rRNA_Hsa"', '"LTR40A1"', '"Eutr5"', '"MamTip1b"', '"Helitron1Nb_Mam"', '"LTR9B_EC"', '"hAT-1_Mam"', '"Chap1a_Mam"', '"X32_DNA"', '"MER67B"', '"L1M6B"', '"EUTREP16"', '"X7C_LINE"', '"MamRep488"', '"FordPrefect"', '"UCON33"', '"LTR35B_CF"', '"X10b_DNA"', '"LTR86A1"', '"L1M3e"', '"MLT1B-int"', '"UCON60"', '"LTR103_Mam"', '"Eutr2B"', '"Eutr11"', '"Eulor12"', '"Tigger8"', '"LTR72B_CF"', '"Tigger11a"', '"Eulor1"', '"U13_"', '"UCON39"', '"LTR65"', '"MER70-int"', '"CarLTR8a"', '"CarERV3_LTR"', '"MER97c"', '"X7D_LINE"', '"hAT-4b_Ther"', '"LTR16D2"', '"L2-1_AMi"', '"LTR81B"', '"Eulor6E"', '"Charlie17"', '"MamTip3B"', '"CarLTR6a"', '"Tigger23a"', '"Chap1_Mam"', '"Tigger21a"', '"L1M3d"', '"MER134"', '"UCON9"', '"UCON58"', '"FordPrefect_a"', '"LTR75B"', '"MLT1F1-int"', '"MER88"', '"CarLTR1A1C"', '"UCON47"', '"CarLTR1A1B"', '"Eulor6D"', '"MLT1K-int"', '"UCON99"', '"Charlie32a"', '"EutTc1-N1"', '"Carsat2"', '"MER101B_CF"', '"U7"', '"MARE10"', '"MLT1E1A-int"', '"tRNA-Met_"', '"UCON2"', '"tRNA-Phe-TTY"', '"U17"', '"U13"', '"tRNA-Arg-AGG"']
