#!/bin/bash

# yeast transcriptome index created with:
#  bowtie2-build yeast_transcriptome_plus_citrine_cds.fa citrine_cds_transcriptome

yeast_idx="../PolysomeRNASeq20230914/Index/citrine_cds_transcriptome"
human_idx="../PolysomeRNASeq20230914/Index/grch38_transcriptome"

aln_dir="Alignments"
seq_dir="SeqData"

# yeast

# allow multi mapped reads, then check later if any of the citrine reads are affected (they aren't)
# example:
# bowtie2 -p 32 -k 4 --un-gz Alignments/B0_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x ../PolysomeRNASeq20230914/Index/citrine_cds_transcriptome -U SeqData/ix_B0_S1_L002_R1_001.fastq.gz -S Alignments/B0_yeast.sam

# B replicate
for i in {1..8}
do
    s=$(($i + 1))
    bowtie2 -p 32 -k 4 --un-gz ${aln_dir}/B${i}_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x ${yeast_idx} -U ${seq_dir}/ix_B${i}_S${s}_L002_R1_001.fastq.gz -S ${aln_dir}/B${i}_yeast.sam
done

# C replicate
for i in {0..8}
do
    s=$(($i + 10))
    bowtie2 -p 32 -k 4 --un-gz ${aln_dir}/C${i}_unaligned_to_yeast.fastq.gz --no-unal --no-sq -x ${yeast_idx} -U ${seq_dir}/ix_C${i}_S${s}_L002_R1_001.fastq.gz -S ${aln_dir}/C${i}_yeast.sam
done

# human

# run in the mode where it returns 0 or 1 hits per read, no multi map - should give a "# of mapped reads" we can use to normalize
# example:
# bowtie2 -p 32 --no-unal --no-sq -x ../../PolysomeRNASeq20230914/Index/grch38_transcriptome -U B0_unaligned_to_yeast.fastq.gz -S B0_human.sam

for j in {0..8}
do
    bowtie2 -p 32 --no-unal --no-sq -x ${human_idx} -U ${aln_dir}/B${j}_unaligned_to_yeast.fastq.gz -S ${aln_dir}/B${j}_human.sam
    bowtie2 -p 32 --no-unal --no-sq -x ${human_idx} -U ${aln_dir}/C${j}_unaligned_to_yeast.fastq.gz -S ${aln_dir}/C${j}_human.sam
done


# yeast analysis
# party like it's 1999!
grep cit ${aln_dir}/B*_yeast.sam | grep "NM:i:0" | perl -ne 's/_/\t/g; print;' | cut -f 1,4 | sort | uniq -c > ${aln_dir}/citrine_counts_B.txt
grep cit ${aln_dir}/C*_yeast.sam | grep "NM:i:0" | perl -ne 's/_/\t/g; print;' | cut -f 1,4 | sort | uniq -c > ${aln_dir}/citrine_counts_C.txt

# human analysis
grep "NM:i:0" ${aln_dir}/B*_human.sam | cut -f 1 -d '_' | uniq -c > ${aln_dir}/human_counts_B.txt
grep "NM:i:0" ${aln_dir}/C*_human.sam | cut -f 1 -d '_' | uniq -c > ${aln_dir}/human_counts_C.txt
