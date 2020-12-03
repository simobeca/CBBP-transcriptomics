# Becattini et al, CHM 2021, 'Rapid transcriptional and metabolic adaptation of intestinal microbes to host immune activation'
# the following bash script was used to generate sam files to be inputted into FeatureCounts;
# MyFolder should contain one folde per each pair of reads (= 1 folder per sample) in .fq.gz format

#!/bin/bash
cd /MyFolder

# make a new directory calles 'SAM_BAM_ALIGNMENTS'
mkdir SAM_BAM_ALIGNMENTS

# the following for loop will trim reads based on quality, align to mouse genome to get rid of mouse derived reads,
# then align to CBBP genomes to get 4 .BAM files that can be used for RNAseq analysis, containing only aligned sequences

#run trimmomatic 
for i in */ ; 
do 
	java -jar PATH_TO_TRIMMOMATIC/Trimmomatic-0.36/trimmomatic-0.36.jar  PE -phred33 -trimlog trimLogFile  $i/*_R1_001.fastq.gz $i/*_R2_001.fastq.gz $i/${i%/}output_forward_paired.fq.gz $i/${i%/}output_forward_unpaired.fq.gz $i/${i%/}output_reverse_paired.fq.gz $i/${i%/}output_reverse_unpaired.fq.gz ILLUMINACLIP:PATH_TO_TRIMMOMATIC/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; 

# cancel original files 
	rm  $i/*_001.fastq.gz ; 

# gunzip files for bowtie2 alignment
	gunzip $i/*_paired.fq.gz ; 

# align to mouse genome, spare unaligned reads (called no_mouse_paired.fq)
	bowtie2 -p 4 -k 1 --dovetail -x PATH_TO_INDEXES/mouse_genome/mouse_genome -1 $i/*output_forward_paired.fq -2 $i/*output_reverse_paired.fq -S $i/${i%/}_alignment_to_mouse.sam --al-conc $i/${i%/}_al_to_mouse_paired.fq --un-conc $i/${i%/}_nomouse_paired.fq; 
# cancel to_mouse alignments 
	rm  $i/*al_to_mouse_paired.*.fq ; 
	rm  $i/*.sam

# cancel trimmomatic output
	rm  $i/*output_*.fq ; 
	rm  $i/*output_*.fq.gz ; 

# align to each of the CBBP genomes, and then transform the .sam files into .bam files containing only aligned reads using samtools
	bowtie2 -p 4 --dovetail -x PATH_TO_INDEXES/BP/BP -1 $i/*nomouse_paired.1.fq -2 $i/*nomouse_paired.2.fq | samtools view -bSF 4 > $i/${i%/}.BP.only_aligned.bam ;
	bowtie2 -p 4 --dovetail -x PATH_TO_INDEXES/BS/BS -1 $i/*nomouse_paired.1.fq -2 $i/*nomouse_paired.2.fq | samtools view -bSF 4 > $i/${i%/}.BS.only_aligned.bam ;
	bowtie2 -p 4 --dovetail -x PATH_TO_INDEXES/CB/CB -1 $i/*nomouse_paired.1.fq -2 $i/*nomouse_paired.2.fq | samtools view -bSF 4 > $i/${i%/}.CB.only_aligned.bam ;
	bowtie2 -p 4 --dovetail -x PATH_TO_INDEXES/PD/PD -1 $i/*nomouse_paired.1.fq -2 $i/*nomouse_paired.2.fq | samtools view -bSF 4 > $i/${i%/}.PD.only_aligned.bam ;

#remove .fq files to save space	and move .bam files to SAM_BAM_ALIGNMENTS folder 
	rm $i/*.fq ;
	mv $i/*.bam /MyFolder/SAM_BAM_ALIGNMENTS ;

done


