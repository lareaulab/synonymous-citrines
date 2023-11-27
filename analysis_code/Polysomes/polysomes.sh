#!/bin/bash

# cd Index
# bowtie2-build yeast_transcriptome_plus_citrine_cds.fa citrine_cds_transcriptome
# cd ..

# # SeqData/LLLD001B0_S1_R1.fastq.gz
# # SeqData/LLLD001B1_S3_R1.fastq.gz
# # SeqData/LLLD001B2_S5_R1.fastq.gz
# # SeqData/LLLD001B3_S7_R1.fastq.gz
# # SeqData/LLLD001B4_S9_R1.fastq.gz
# # SeqData/LLLD001B5_S11_R1.fastq.gz
# # SeqData/LLLD001C0_S2_R1.fastq.gz
# # SeqData/LLLD001C1_S4_R1.fastq.gz
# # SeqData/LLLD001C2_S6_R1.fastq.gz
# # SeqData/LLLD001C3_S8_R1.fastq.gz
# # SeqData/LLLD001C4_S10_R1.fastq.gz
# # SeqData/LLLD001C5_S12_R1.fastq.gz

# allow multi mapped reads, sort out later if any of the citrine reads are affected
# bowtie2 -p 16 -k 4 --un Alignments/B0_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B0_S1_R1.fastq.gz -S Alignments/B0_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/B1_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B1_S3_R1.fastq.gz -S Alignments/B1_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/B2_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B2_S5_R1.fastq.gz -S Alignments/B2_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/B3_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B3_S7_R1.fastq.gz -S Alignments/B3_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/B4_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B4_S9_R1.fastq.gz -S Alignments/B4_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/B5_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001B5_S11_R1.fastq.gz -S Alignments/B5_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C0_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C0_S2_R1.fastq.gz -S Alignments/C0_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C1_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C1_S4_R1.fastq.gz -S Alignments/C1_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C2_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C2_S6_R1.fastq.gz -S Alignments/C2_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C3_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C3_S8_R1.fastq.gz -S Alignments/C3_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C4_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C4_S10_R1.fastq.gz -S Alignments/C4_yeast.sam
# bowtie2 -p 32 -k 4 --un Alignments/C5_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x Index/citrine_cds_transcriptome -U SeqData/LLLD001C5_S12_R1.fastq.gz -S Alignments/C5_yeast.sam

cd Alignments
# human transcriptome index copied from ~/BioE290/Index/bowtie2/ ...
# run in the mode where it returns 0 or 1 hits per read, no multi map - should give a "# of mapped reads" we can use to normalize
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B0_unaligned_to_yeast.fastq.gz -S B0_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B1_unaligned_to_yeast.fastq.gz -S B1_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B2_unaligned_to_yeast.fastq.gz -S B2_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B3_unaligned_to_yeast.fastq.gz -S B3_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B4_unaligned_to_yeast.fastq.gz -S B4_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U B5_unaligned_to_yeast.fastq.gz -S B5_human.sam

# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C0_unaligned_to_yeast.fastq.gz -S C0_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C1_unaligned_to_yeast.fastq.gz -S C1_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C2_unaligned_to_yeast.fastq.gz -S C2_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C3_unaligned_to_yeast.fastq.gz -S C3_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C4_unaligned_to_yeast.fastq.gz -S C4_human.sam
# bowtie2 -p 32 --no-unal --no-sq -x ../Index/grch38_transcriptome -U C5_unaligned_to_yeast.fastq.gz -S C5_human.sam

grep "NM:i:0" C0_human.sam > C0_human_perfect.sam
grep "NM:i:0" C1_human.sam > C1_human_perfect.sam
grep "NM:i:0" C2_human.sam > C2_human_perfect.sam
grep "NM:i:0" C3_human.sam > C3_human_perfect.sam
grep "NM:i:0" C4_human.sam > C4_human_perfect.sam
grep "NM:i:0" C5_human.sam > C5_human_perfect.sam

grep "NM:i:0" B0_human.sam > B0_human_perfect.sam
grep "NM:i:0" B1_human.sam > B1_human_perfect.sam
grep "NM:i:0" B2_human.sam > B2_human_perfect.sam
grep "NM:i:0" B3_human.sam > B3_human_perfect.sam
grep "NM:i:0" B4_human.sam > B4_human_perfect.sam
grep "NM:i:0" B5_human.sam > B5_human_perfect.sam

grep cit9 B0_citrines.sam | grep "NM:i:0" > B0_cit9.sam
grep cit9 B1_citrines.sam | grep "NM:i:0" > B1_cit9.sam
grep cit9 B2_citrines.sam | grep "NM:i:0" > B2_cit9.sam
grep cit9 B3_citrines.sam | grep "NM:i:0" > B3_cit9.sam
grep cit9 B4_citrines.sam | grep "NM:i:0" > B4_cit9.sam
grep cit9 B5_citrines.sam | grep "NM:i:0" > B5_cit9.sam

grep cit9 C0_citrines.sam | grep "NM:i:0" > C0_cit9.sam
grep cit9 C1_citrines.sam | grep "NM:i:0" > C1_cit9.sam
grep cit9 C2_citrines.sam | grep "NM:i:0" > C2_cit9.sam
grep cit9 C3_citrines.sam | grep "NM:i:0" > C3_cit9.sam
grep cit9 C4_citrines.sam | grep "NM:i:0" > C4_cit9.sam
grep cit9 C5_citrines.sam | grep "NM:i:0" > C5_cit9.sam

grep citmin B0_citrines.sam | grep "NM:i:0" > B0_citmin.sam
grep citmin B1_citrines.sam | grep "NM:i:0" > B1_citmin.sam
grep citmin B2_citrines.sam | grep "NM:i:0" > B2_citmin.sam
grep citmin B3_citrines.sam | grep "NM:i:0" > B3_citmin.sam
grep citmin B4_citrines.sam | grep "NM:i:0" > B4_citmin.sam
grep citmin B5_citrines.sam | grep "NM:i:0" > B5_citmin.sam

grep citmin C0_citrines.sam | grep "NM:i:0" > C0_citmin.sam
grep citmin C1_citrines.sam | grep "NM:i:0" > C1_citmin.sam
grep citmin C2_citrines.sam | grep "NM:i:0" > C2_citmin.sam
grep citmin C3_citrines.sam | grep "NM:i:0" > C3_citmin.sam
grep citmin C4_citrines.sam | grep "NM:i:0" > C4_citmin.sam
grep citmin C5_citrines.sam | grep "NM:i:0" > C5_citmin.sam

wc -l *cit9.sam|grep -v total |cut -f 1 -d "_" |perl -ne '($num, $frac) = /(\d+) ([BC]\d)/; print $frac . "\t" . $num . "\n"' > cit9_counts.txt
wc -l *citmin.sam|grep -v total |cut -f 1 -d "_" |perl -ne '($num, $frac) = /(\d+) ([BC]\d)/; print $frac . "\t" . $num . "\n"' > citmin_counts.txt

wc -l *human_perfect.sam |grep -v total |cut -f 1 -d "_" |perl -ne '($num, $frac) = /(\d+) ([BC]\d)/; print $frac . "\t" . $num . "\n"' > human_perfect_totals.txt
