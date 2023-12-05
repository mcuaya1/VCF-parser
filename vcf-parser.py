import os
import gzip
import csv
import argparse
from datetime import datetime
"""
TODO:
	Add Frequency count of mutation spotted, perhaps using another dictonary as a value with
	it's key being the mutation and the value being the amount of times it has been seen
		*CURRENT Wokring approach is simply having the control sample be the first sample in VCF file with all other samples coming before the control sample.
		This is done when creating the vcf file which takes in a control sample if present and ensures this sample is the first sample
	Work with new VCF file to be generated from new BAM file
		- Will include other samples as well

	*Add user defined output argument
	Add more comments to code
	Add output to readme
	Add description to readme
	Add option if no control same is given
	*Add pandas data frame to sort Frequency column
"""
parser = argparse.ArgumentParser(add_help=False, prog="vcf-parser.py", description="Program to find the frequency of mutations in a vcf file")

parser.add_argument('-i',"--input", required=True, help="Input generated tsv stats file",type=str)
parser.add_argument('-a', "--added-info", action="store_true", help="add extra info such as sample name, chromosomes, and more")
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program')
args = parser.parse_args()


#Store arguments in variables
input_file = args.input
added_info = args.added_info


sample_column=[]
print("Read file in...")
mutation_list = {}

with open(input_file, "r") as vcf_file:
	vcf_parser=csv.reader(vcf_file,delimiter='\t')
	mutation=()
	mutation_cnt=0

	for row in vcf_parser:
		for i in range(5, len(row)):
			sample_column.append(i)
			print(i)
		print(f"There are {len(sample_column)} samples")

		break
	vcf_file.seek(0)
	#next(vcf_parser)

	for column in vcf_parser:
		mutation_cnt+=1
		control=column[4].split('=')[1].strip()
		chromosome_location=column[0]
		#For loop that loops through each sample and compares that samples mutations with the control sample
		for i in range(0, len(sample_column)):

			#The next 4 lines simply remove the directory path of the sample, keeping only the sample name
			cur_sample=column[sample_column[i]].split('=')[1].strip()
			sample_name=column[sample_column[i]].split('=')[0].strip()
			if '/' in sample_name:
				sample_name=sample_name.split('/')[-1].strip().split('.')[0].strip()


			mutation = (control, cur_sample)
			if mutation not in mutation_list and (control != cur_sample) and (control != "." and cur_sample != ".") and (len(control)==1 and len(cur_sample)==1):
				mutation_list[mutation]={}
				mutation_list[mutation]["NS"]=1
				mutation_list[mutation]["FQ"]=1
				mutation_list[mutation]["SAMPLE"]={}
				mutation_list[mutation]["SAMPLE"][sample_name]=1
				mutation_list[mutation]["LOCATION"]=[]
				mutation_list[mutation]["LOCATION"].append(chromosome_location)
			elif mutation in mutation_list:
				if mutation_list[mutation]["NS"] != len(sample_column):
					mutation_list[mutation]["NS"]+=1

				if chromosome_location not in mutation_list[mutation]["LOCATION"]:
					mutation_list[mutation]["LOCATION"].append(chromosome_location)

				if sample_name in mutation_list[mutation]["SAMPLE"]:
					mutation_list[mutation]["SAMPLE"][sample_name]+=1
				else:
					mutation_list[mutation]["SAMPLE"][sample_name]=1
				mutation_list[mutation]["FQ"]+=1


output_file="VCF_PARSER_OUTPUT.tsv"
time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")
with open(output_file, "w") as f:
	print(f"#File name:{output_file}", file=f)
	print("#Column Descriptions:", file=f)
	print("#NS='Number of samples mutation was seen in'",file=f)
	print("#FQ='Frequenecy of mutation seen'",file=f)
	print(f"#Date file was generated:{time_generated}",file=f)
	if added_info == True:
		print("#CTRL\tALT\tSEEN\tFREQUENCY\tSAMPLES\tLOCATION",file=f)
		for mutation,info in mutation_list.items():
			print(f"{mutation[0]}\t{mutation[1]}",end="\t",file=f)
			for key,value in info.items():
				if(key=='FQ'):
					print(f"{key}={((value/mutation_cnt)*100):.2f}",end="\t",file=f)
				else:
					print(f"{key}={value}",end="\t",file=f)
			f.write("\n")
	else:
		print("#CTRL\tALT\tSEEN\tFREQUENCY",file=f)
		for mutation,info in mutation_list.items():
			print(f"{mutation[0]}\t{mutation[1]}",end="\t",file=f)
			for key,value in info.items():
				if(key=='FQ'):
					print(f"{key}={((value/mutation_cnt))*100:.2f}",end="\t",file=f)
				elif(key != 'SAMPLE' and key !='LOCATION'):
					print(f"{key}={value}",end="\t",file=f)
			f.write("\n")
