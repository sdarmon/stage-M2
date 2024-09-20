#!/bin/bash

        M=$(ls /beegfs/data/sdarmon/results/mus/processing/comp*.txt | wc -l)
	N=$(ls /beegfs/data/sdarmon/results/mus/processing/comp*_TE.txt | wc -l)
	MAXI=$((${M}-${N}))
	echo "Number of comps : ${MAXI}"
        for ((i=0; i<$MAXI; i++))
        do
		echo $i
            	python3 \
                reads_to_align.py /beegfs/data/sdarmon/results/mus/processing/comp$i.txt \
                 /beegfs/data/sdarmon/results/mus/processing/comp$i.fq \
                0
            mkdir -p /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i
            mkdir -p /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_cons_$i

            /beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STARlong \
	        --genomeDir /beegfs/data/sdarmon/results/mus/ref \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn /beegfs/data/sdarmon/results/mus/processing/comp$i.fq  \
                --outFileNamePrefix /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx

            /beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STARlong \
	        --genomeDir /beegfs/data/sdarmon/results/mus/te_ref \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn /beegfs/data/sdarmon/results/mus/processing/comp$i.fq  \
                --outFileNamePrefix /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_cons_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
            
	    /beegfs/home/sdarmon/Documents/samtools-1.9/samtools view -h -F 256 /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam | /beegfs/home/sdarmon/Documents/samtools-1.9/samtools view -b -o /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/output_filtered.bam
	    

	/beegfs/home/sdarmon/Documents/samtools-1.9/samtools view -F 256 /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_cons_$i/Aligned.sortedByCoord.out.bam  > /beegfs/data/sdarmon/results/mus/processing/seq_intersectionCons$i.txt

	   # bedtools intersect -wa -a /beegfs/data/sdarmon/mus/TE.gtf -b /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/output_filtered.bam -split > /beegfs/data/sdarmon/results/mus/processing/intersectionTE$i.txt
	   echo "" > /beegfs/data/sdarmon/results/mus/processing/intersectionTE$i.txt

	/beegfs/home/sdarmon/Documents/bedtools2/bin/bedtools intersect -wa -a /beegfs/data/sdarmon/mus/ref.gtf -b /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/output_filtered.bam -split > /beegfs/data/sdarmon/results/mus/processing/intersectionRef$i.txt
	#bedtools intersect -wb -a /beegfs/data/sdarmon/mus/TE.gtf -b /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/output_filtered.bam  -split > /beegfs/data/sdarmon/results/mus/processing/seq_intersectionTE$i.txt
	   echo "" > /beegfs/data/sdarmon/results/mus/processing/seq_intersectionTE$i.txt

/beegfs/home/sdarmon/Documents/bedtools2/bin/bedtools intersect -wb -a /beegfs/data/sdarmon/mus/ref.gtf -b /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/output_filtered.bam  -split > /beegfs/data/sdarmon/results/mus/processing/seq_intersectionRef$i.txt

	python3 analysis_comp.py /beegfs/data/sdarmon/results/mus/processing/comp$i.fq /beegfs/data/sdarmon/results/mus/processing/intersectionRef$i.txt /beegfs/data/sdarmon/results/mus/processing/seq_intersectionCons$i.txt /beegfs/data/sdarmon/results/mus/processing/analysis$i.txt --custom-sep "split('\t')[8].split(';')[0][8:-1]"

python3 analysis_comp.py /beegfs/data/sdarmon/results/mus/processing/comp$i.fq /beegfs/data/sdarmon/results/mus/processing/intersectionRef$i.txt /beegfs/data/sdarmon/results/mus/processing/seq_intersectionCons$i.txt /beegfs/data/sdarmon/results/mus/processing/cr_analysis$i.txt --custom-sep "split('\t')[8].split(';')[0][8:-1]" --short 
	
	python3 add_ref_TE.py /beegfs/data/sdarmon/results/mus/processing/comp$i.txt /beegfs/data/sdarmon/results/mus/processing/seq_intersectionTE$i.txt /beegfs/data/sdarmon/results/mus/processing/seq_intersectionRef$i.txt /beegfs/data/sdarmon/results/mus/processing/seq_intersectionCons$i.txt  results/graph_hc_mus_1_hc_mus_2_k41.abundance

	/beegfs/home/sdarmon/Documents/samtools-1.9/samtools index /beegfs/data/sdarmon/results/mus/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam 
	
done


