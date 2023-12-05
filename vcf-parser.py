import os
import gzip
import csv
import argparse
from datetime import datetime
"""
TODO:
	Add Frequency count of mutation spotted, perhaps using another dictonary as a value with
	it's key being the mutation and the value being the amount of times it has been seen
		-Currently working based off of sample_name being the key to another dictonary
		where the mutation seen in that sample is the mutation
		-A lot of this is hardcoded with column numbers being based off of where I am assuming a control sample will be located
		which I am currently assuming will be the first column after the ALT column
		-A better solution to approach this problem might be to use a conditional or argument that asks "Does the file contain a control sample?"
		and if so "Which column is the control sample located?" assuming it's not in the first column after ALT
	Work with new VCF file to be generated from new BAM file
		- Will include other samples as well

	Added new parameter for amount of sampes(Done)
"""
parser = argparse.ArgumentParser(add_help=False, prog="vcf-parser.py", description="Program to find the frequency of mutations in a vcf file")

parser.add_argument('-i',"--input", required=True, help="Input generated tsv stats file",type=str)
parser.add_argument('-n', "--number-of-samples", type=int, help="Number of samples in generated tsv stats file")
parser.add_argument('-a', "--added-info", action="store_true", help="add extra info such as sample name, chromosomes, and more")
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program')
args = parser.parse_args()


#Store arguments in variables
input_file = args.input
num_of_samples = args.number_of_samples
added_info = args.added_info
"""
Create a list containing the location of samples within the input file
This is currently hardcoded to assume that the first column after the ALT column(column[3]) is the column containing the
Control sample(column[4]) so every other sample is after this column
"""
sample_column=[]
for i in range(5,5+num_of_samples):
	sample_column.append(i)


print("Read file in...\n")
mutation_list = {}

with open(input_file, "r") as vcf_file:
	vcf_parser=csv.reader(vcf_file,delimiter='\t')
	mutation=()
	mutation_cnt=0
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
			if mutation not in mutation_list and (control != cur_sample) and (control != "." and cur_sample != "."):
				mutation_list[mutation]={}
				mutation_list[mutation]["NS"]=1
				mutation_list[mutation]["FQ"]=1
				mutation_list[mutation]["SAMPLE"]={}
				mutation_list[mutation]["SAMPLE"][sample_name]=1
				mutation_list[mutation]["LOCATION"]=[]
				mutation_list[mutation]["LOCATION"].append(chromosome_location)
			elif mutation in mutation_list:
				if mutation_list[mutation]["NS"] != num_of_samples:
					mutation_list[mutation]["NS"]+=1

				if chromosome_location not in mutation_list[mutation]["LOCATION"]:
					mutation_list[mutation]["LOCATION"].append(chromosome_location)

				if sample_name in mutation_list[mutation]["SAMPLE"]:
					mutation_list[mutation]["SAMPLE"][sample_name]+=1
				else:
					mutation_list[mutation]["SAMPLE"][sample_name]=1
				mutation_list[mutation]["FQ"]+=1
			"""
			If statement that checks if a sample_name is in the mutation list.
			Statement is not need after all the samples have been put in a list so I might remove this later and
			have a for loop that does this step once so this check is not done every loop
			"""
#			if sample_name not in mutation_list:
#				mutation_list[sample_name] = {}

			"""
			Based on the sample name, this bit of code checks if the mutation is in the current  sample
			and that the mutation is indeed a mutation as well as ensuring uncertain mutations aren't passing through.
			If all the checks pass then it is put into the sample's dictonary of mutations and added
			"""
#			if (mutation not in mutation_list[sample_name]) and (control != cur_sample) and (control != '.' and cur_sample != '.'):
#				mutation = (control, cur_sample)
#				mutation_list[sample_name][mutation] = 1
#			elif mutation in mutation_list[sample_name]:
#				mutation_list[sample_name][mutation] += 1


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
					print(f"{key}={((value/mutation_cnt)*100):.2f}",end="\t",file=f)
				elif(key != 'SAMPLE' and key !='LOCATION'):
					print(f"{key}={value}",end="\t",file=f)
			f.write("\n")
