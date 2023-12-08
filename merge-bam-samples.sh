#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 1 -c 64 --mem 32gb
#SBATCH --job-name="bam-merge"
#SBATCH --error="merge_error.txt"


while getopts "i:ho:b:" flag; do
	case $flag in
		i)
		ID="$OPTARG"
		;;
		b)
		BAM_DIR="$OPTARG"
		;;
		o)
		OUTPUT="$OPTARG"
		;;
		h)
		echo "MERGE-CONTROL"
		echo "DESCRIPTION"
		echo "		Script that merges bam files"
		echo "OPTIONS"
		echo "	-i"
		echo "		Sample identifer."
		echo "	-o"
		echo "		Custom output directory."
		echo "	-b"
		echo "		Directory where bam files are stored. This is necessary if you want all bam files in a directory to be used."
		exit 1
		;;
		\?)
		echo "Unknown error"
		exit 1
		;;
	esac
done


module load samtools
CPUs=64


if [ ! -z $OUTPUT ]; then
	echo "Creating user defined directory..."
	mkdir ${OUTPUT}
else
	echo "Using default directory..."
	OUTPUT=MERGED_BAM
fi


INPUT=$(ls ${BAM_DIR}/ | grep "${ID}" | grep -v ".bai" | awk -v dir=${BAM_DIR} '{ print dir"/" $0}')

echo "Merging samples..."
samtools merge -f -@ "$CPUs" -o "$OUTPUT/${ID}.bam" ${INPUT}
samtools sort --threads $CPU -T $SCRATCH -O BAM -o ${OUTPUT}/${ID}.bam ${OUTPUT}/${ID}_fixmate.bam
samtools index ${OUTPUT}/${ID}.bam
rm ${OUTPUT}/${ID}_fixmate.bam
echo "Finished..."

