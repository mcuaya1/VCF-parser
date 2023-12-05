#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 1 -c 64 --mem 32gb
#SBATCH --job-name="vcf-output"
#SBATCH --error="vcf_error.txt"

while getopts "f:e:ho:d:" flag; do
	case $flag in
		f)
		INPUT="$OPTARG"
		;;
		e)
		ALL_FILES="$OPTARG"
		;;
		d)
		GIVEN_DIR="$OPTARG"
		;;
		o)
		OUTPUT="$OPTARG"
		;;
		h)
		echo "BAM-TO-VCF"
		echo "DESCRIPTION"
		echo "		Script that produces vcf file from bam file(s)"
		echo "OPTIONS"
		echo "	-f"
		echo "		Takes in a string of files to convert"
		echo "	-e"
		echo "		Takes in a 0 or 1 input where 0 tells the program not to parse through every bam file in a the directory and 1 does parse through every bam file in the directory"
		echo "	-o"
		echo "		Output file name"
		echo "	-d"
		echo "		Directory where bam files are stored. This is necessary if you want all bam files in a directory to be used"
		exit 1
		;;
		\?)
		echo "Unknown error"
		exit 1
		;;
	esac
done



if [ $ALL_FILES == 0 ]
then
	echo "Given files will be used..."
	BAM_DIR=$INPUT
else
	echo "All files in given directory will be used..."
	BAM_DIR=${GIVEN_DIR}/*
fi


if [ ! -d $PWD/VCF_OUTPUT ];
then
	echo "Creating output directory..."
	mkdir VCF_OUTPUT
fi

module load samtools
module load bcftools
GENOME=S288C_reference_sequence_R64-4-1_20230823.fsa
VCF_DIR=VCF_OUTPUT
# need to make a string which is all the bam files you want to process
# but if we do *.bam it will catch the intermediate bam files that are in the folder
VCF=${OUTPUT}.vcf
VCFFILTER=${OUTPUT}.vcf.filtered

bcftools mpileup -Ou -f $GENOME $BAM_DIR | bcftools call -vmO z -o $VCF_DIR/$VCF --ploidy 1
tabix -p vcf $VCF_DIR/$VCF
gzip $VCF_DIR/$VCF.tabi
bcftools stats -F $GENOME -s - $VCF_DIR/$VCF > $VCF_DIR/$VCF.stats
mkdir -p $VCF_DIR/plots
plot-vcfstats -p $VCF_DIR/plots/$VCF.stats
bcftools filter -O z -o $VCF_DIR/$VCFFILTER -s LOWQUAL -i'%QUAL>10' $VCF_DIR/$VCF
gzip $VCF_DIR/$VCF
gzip $VCF_DIR/$VCF.stats
gzip $VCF_DIR/$VCFFILTER.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%TGT]\n' VCF_OUTPUT/UV_sample.vcf.gz > $VCF_DIR/VCF_OUTPUT_INFO.tsv
echo "Finished..."
