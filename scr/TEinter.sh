#!/bin/bash
MAXI=$(ls ../../results/droso/TE/processing/comp*.txt | wc -l)
echo $MAXI
for ((i=0;i<$(($MAXI-1));i++));
do
        echo $i    
	python3 reads_to_align.py ../../results/droso/TE/processing/comp$i.txt \
                ../../results/droso/TE/processing/comp$i.fq 0
            mkdir -p ../../results/droso/TE/processing/STAR_alignment
            STAR --genomeDir ../../results/droso/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/droso/TE/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/droso/TE/processing/STAR_alignment/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
     bedtools intersect -wa -a ../../data/droso/TE.gtf  -b ../../results/droso/TE/processing/STAR_alignment/Aligned.sortedByCoord.out.bam -split  > ../../results/droso/TE/processing/intersectionTE$i.txt


     bedtools intersect -wb -a ../../data/droso/TE.gtf  -b ../../results/droso/TE/processing/STAR_alignment/Aligned.sortedByCoord.out.bam -split  > ../../results/droso/TE/processing/seq_intersectionTE$i.txt

     python3 add_ref_TE.py ../../results/droso/TE/processing/comp$i.txt ../../results/droso/TE/processing/seq_intersectionTE$i.txt

done
