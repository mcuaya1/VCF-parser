#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 1 -c 16 --mem 16gb
#SBATCH --job-name="fasterq_job"
#SBATCH --error="job_dump_error.txt"


while getopts "f:s:a:d:h" flag; do
        case $flag in
                f)
                INPUT="$OPTARG"
                ;;
                s)
                SAMPLE_COL="$OPTARG"
                ;;
                a)
                ACC_COL="$OPTARG"
                ;;
                d)
                OUTPUT_DIR="$OPTARG"
                ;;
                h)
                echo "BAM-TO-VCF"
                echo "DESCRIPTION"
                echo "          Script that downloads SRA data from an SRA metadata file"
                echo "OPTIONS"
                echo "	-f"
                echo "		SRA metadata file"
                echo "	-s"
                echo "		Column where sample names are located"
                echo "	-a"
                echo "		Column where accession runs are located"
		echo "	-d"
		echo "		Output directory name"
		exit 1
                ;;
                \?)
                echo "Unknown error"
                exit 1
                ;;
        esac
done



module load sratoolkit
echo "Creating necessary files..."
awk -v sample_col="${SAMPLE_COL}" 'BEGIN {FS=","}; NR>1 {print $sample_col}' ${INPUT} | sed 's/ /_/g' > SAMPLE_NAMES.txt
awk -v acc_col="${ACC_COL}" 'BEGIN {FS=","}; NR > 1 {print $acc_col}' ${INPUT} > SRA_ACC.txt

PREFETCH_DIR=prefetch_dump
LIST_OF_SAMPLE_NAMES=SAMPLE_NAMES.txt
LIST_OF_SRAS=SRA_ACC.txt
FASTQ_DIR=$OUTPUT_DIR

#echo "Creating necessary directories..."
#if [ ! -d $FASTQ_DIR ]; then
	mkdir ${FASTQ_DIR}
#fi

#if [ ! -d $PREFETCH_DIR ]; then
	mkdir ${PREFETCH_DIR}
#fi

#for sra in $(cat  $LIST_OF_SRAS); do
#	prefetch $sra -O $PREFETCH_DIR
#done

echo "Downloaded prefetch..."
echo "Unpacking FASTQ data..."
while read -r -u 3 sra_read && read -r -u 4 sample_name; do
	fasterq-dump --skip-technical -o ${sra_read}_${sample_name} -O $FASTQ_DIR -t $SCRATCH $PREFETCH_DIR/$sra_read/$sra_read.sra
	gzip --fast $FASTQ_DIR/${sra_read}_${sample_name}_1.fastq
	gzip --fast $FASTQ_DIR/${sra_read}_${sample_name}_2.fastq
done  3<$LIST_OF_SRAS 4<$LIST_OF_SAMPLE_NAMES
echo "Unpacked FASTQ data..."
