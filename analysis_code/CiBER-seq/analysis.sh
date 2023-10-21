# some variables
SLOW='CCTTG$'
FAST='TCAGA$'

TURB="L R"
TYPE="DNA_ RNA_"
SAMPLES="post_ pre_"
SF="slow_ fast_"

RAWDIR="RawSeqs/"
TRIMDIR="TrimmedSeqs/"
COUNTDIR="Counts/"

# rename
mv LLEL001_A*.fastq.gz ${RAWDIR}DNA_Pre_L.fastq.gz
mv LLEL001_B*.fastq.gz ${RAWDIR}DNA_Post_L.fastq.gz
mv LLEL001_C*.fastq.gz ${RAWDIR}DNA_Pre_R.fastq.gz
mv LLEL001_D*.fastq.gz ${RAWDIR}DNA_Post_R.fastq.gz

mv LLEL001_E*.fastq.gz ${RAWDIR}RNA_Pre_L.fastq.gz
mv LLEL001_F*.fastq.gz ${RAWDIR}RNA_Post_L.fastq.gz
mv LLEL001_G*.fastq.gz ${RAWDIR}RNA_Pre_R.fastq.gz
mv LLEL001_H*.fastq.gz ${RAWDIR}RNA_Post_R.fastq.gz

# trim
for a in ${TURB}
do
    for b in ${SAMPLES}
    do
        for c in ${TYPE}
        do
	    
	    zcat $c$b$a.fastq.gz | \
		cutadapt -a slow="${SLOW}" -a fast="${FAST}" -u -4 -l 25 -o ${TRIMDIR}{name}_$c$b$a.fastq -
	    
        done
    done
done

count

for a in ${TURB}
do
    for b in ${SAMPLES}
    do
        for c in ${TYPE}
        do
	    for d in ${SF}
	    do

	        python bc-count.py -i ${TRIMDIR}$d$c$b$a.fastq -o ${COUNTDIR}$d$c$b$a.count.txt

	    done
	done
    done
done

python bc-tabulate.py -o counts_all.txt ${COUNTDIR}slow* ${COUNTDIR}fast*
