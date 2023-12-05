#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 1 -c 64 --mem 32gb
#SBATCH --job-name="vcf-output"
#SBATCH --error="vcf_error.txt"

while getopts "c:ho:b:d:i:" flag; do
	case $flag in
		c)
		CONTROL="$OPTARG"
		;;
		b)
		GIVEN_DIR="$OPTARG"
		;;
		o)
		OUTPUT="$OPTARG"
		;;
		d)
		OUTPUT_DIR="$OPTARG"
		;;
		i)
		SPECIFIC_FILES="$OPTARG"
		;;
		h)
		echo "BAM-TO-VCF"
		echo "DESCRIPTION"
		echo "		Script that produces vcf file from bam file(s)"
		echo "OPTIONS"
		echo "	-c"
		echo "		Include control sample when creating vcf file."
		echo "	-o"
		echo "		Output file name."
		echo "	-b"
		echo "		Directory where bam files are stored. This is necessary if you want all bam files in a directory to be used."
		echo "	-d"
		echo "		Custom directory where vcf outfiles will be stored."
		echo "	-i"
		echo "		Specific bam files to look for when creating vcf file, should be used with -c parameter"
		exit 1
		;;
		\?)
		echo "Unknown error"
		exit 1
		;;
	esac
done


USING_CTRL=0
BAM_DIR=
if [ ! -z ${CONTROL} ]
then
	echo "Control sample given..."
	echo "Include control sample in VCF creation..."
	ALL_BAM_FILES=$(ls ${GIVEN_DIR}/ | grep -v ".bai" |grep "${SPECIFIC_FILES}.bam" | awk -v dir=${GIVEN_DIR} '{ print dir"/" $0}')
	BAM_DIR=${ALL_BAM_FILES}
	USING_CTRL=1
else
	echo "No control sample given..."
	echo "Preform standard VCF creation..."
	BAM_DIR=$(ls ${GIVEN_DIR}/ | grep "${SPECIFIC_FILES}.bam" | grep -v ".bai" | awk -v dir=${GIVEN_DIR} '{ print dir"/" $0}')
fi

if [ -z ${OUTPUT_DIR} ]
then
	echo "Using default directory..."
	if [ ! -d $PWD/VCF_OUTPUT ]
	then
    		echo "Creating output directory..."
        	mkdir VCF_OUTPUT
		VCF_DIR=VCF_OUTPUT
	fi

else
	echo "Directory ${OUTPUT_DIR}... will be used..."
	if [ ! -d $PWD/${OUTPUT_DIR} ]
	then
		echo "Creating user defined output directory..."
		mkdir ${OUTPUT_DIR}
	fi

	VCF_DIR=${OUTPUT_DIR}
fi


module load samtools
module load bcftools
GENOME=S288C_reference_sequence_R64-4-1_20230823.fsa

VCF=${OUTPUT}.vcf
VCFFILTER=${OUTPUT}.vcf.filtered

if [ ${USING_CTRL} == 0 ]
then
	echo "Only using bam files and not control sample..."
	bcftools mpileup -Ou -f $GENOME $BAM_DIR | bcftools call -vmO z -o $VCF_DIR/$VCF --ploidy 1
else
	echo "Using control sample and bam files..."
	bcftools mpileup -Ou -f $GENOME $CONTROL $BAM_DIR | bcftools call -vmO z -o $VCF_DIR/$VCF --ploidy 1
fi

tabix -p vcf $VCF_DIR/$VCF
gzip $VCF_DIR/$VCF.tabi
bcftools stats -F $GENOME -s - $VCF_DIR/$VCF > $VCF_DIR/$VCF.stats
mkdir -p $VCF_DIR/plots
plot-vcfstats -p $VCF_DIR/plots/$VCF.stats
bcftools filter -O z -o $VCF_DIR/$VCFFILTER -s LOWQUAL -i'%QUAL>10' $VCF_DIR/$VCF
echo "Creating custom stats file..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%TGT]\n' ${VCF_DIR}/$VCF > $VCF_DIR/VCF_OUTPUT_INFO.tsv

echo "Compressing files..."
gzip $VCF_DIR/$VCF
gzip $VCF_DIR/$VCF.stats
gzip $VCF_DIR/$VCFFILTER

echo "Finished..."
