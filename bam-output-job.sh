#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 1 -c 64 --mem 32gb
#SBATCH --job-name="bam-output"
#SBATCH --error="bam-job-error.txt"

while getopts "hd:f:" flag; do
	case $flag in
		f)
		GIVEN_DIR="$OPTARG"
		;;
		d)
		OUTPUT_DIR="$OPTARG"
		;;
		h)
		echo "BAM-TO-VCF"
		echo "DESCRIPTION"
		echo "		Script that produces vcf file from bam file(s)"
		echo "OPTIONS"
		echo "	-f"
		echo "		Directory where fastq files are stored."
		echo "	-d"
		echo "		Custom directory where alignment files will be stored."
		exit 1
		;;
		\?)
		echo "Unknown error"
		exit 1
		;;
	esac
done


module load bwa
module load samtools
CPU=64
GENOME=S288C_reference_sequence_R64-4-1_20230823.fsa

SAMPLES=$(ls ${GIVEN_DIR}/ | grep -v "_2" | sed 's/_1.fastq.gz//g' | sort)


if [ ! -f $GENOME.sa ]; then
   bwa index $GENOME
fi

if [ ! -d ${OUTPUT_DIR} ];then
	mkdir ${OUTPUT_DIR}
fi

for acc in ${SAMPLES}
do
	FWDREAD=${GIVEN_DIR}/${acc}_1.fastq.gz
	REVREAD=${GIVEN_DIR}/${acc}_2.fastq.gz

	bwa mem -t $CPU $GENOME $FWDREAD $REVREAD > ${OUTPUT_DIR}/${acc}.sam
	samtools fixmate -O bam ${OUTPUT_DIR}/${acc}.sam ${OUTPUT_DIR}/${acc}_fixmate.bam
	samtools sort --threads $CPU -O BAM -o ${OUTPUT_DIR}/${acc}.bam ${OUTPUT_DIR}/${acc}_fixmate.bam
	samtools index ${OUTPUT_DIR}/${acc}.bam
	rm ${OUTPUT_DIR}/${acc}_fixmate.bam ${OUTPUT_DIR}/${acc}.sam
done
