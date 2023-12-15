import os
import csv
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime


parser = argparse.ArgumentParser(add_help=False, prog="vcf-parser.py", description="Program to find the frequency of mutations in a vcf file.")

parser.add_argument('-i',"--input-file", required=True, help="Input file to be used by program",type=str)
parser.add_argument('-d', "--input-dir", help="Input file directory")
parser.add_argument('-a', "--added-info", 
                    action="store_true", 
                    help="Add extra info such as the following\n\tChromosome location: Chromosome location where mutation appears as well as frequency of mutation in locaiton.")
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Display commands possible with this program.')
parser.add_argument('-c', "--no_control", action="store_true", help="If the input file contains a control sample then enable this flag.")
parser.add_argument('-o', "--output-file", required=True, help="Name of output file", type=str)
args = parser.parse_args()


#Store arguments in variables
input_file = args.input_file
added_info = args.added_info
output_file = args.output_file
no_control = args.no_control
input_dir = args.input_dir

sample_row=[]
mutation_list = {}


#Function to determine SNPs type
def type_of_mutation(ref_base: str, alt_base: str) -> str:
    type= ""
    if(len(alt_base) == 1):
        if(ref_base == "A" or ref_base == "G"):
            if(alt_base == "G" or alt_base == "A"):
                type = "Transition"
            else:
                type = "Transversion"
        elif(ref_base == "T" or ref_base =="C"):
            if(alt_base == "C" or alt_base =="T"):
                type = "Transition"
            else:
                type = "Transversion"
    else:
        if len(alt_base) > len(ref_base):
            type = "Insertion"
        else:
            type = "Deletion"
    return type


working_file = input_file
if '/' in input_file:
    working_file=os.path.join(input_file)    


with open(working_file, "r") as vcf_file:
    vcf_parser=csv.reader(vcf_file,delimiter='\t')
    mutation=()
    mutation_cnt=0

    if no_control == False:
        for row in vcf_parser:
            for i in range(4, len(row)):
                sample_row.append(i)
            break
        print(f"There are {len(sample_row)} samples")
        vcf_file.seek(0)

        for row in vcf_parser:
            ref_alle=row[4].split('=')[1].strip()
            chromosome_location=row[0]
        
            for i in range(0, len(sample_row)):
                alt_alle=row[sample_row[i]].split('=')[1].strip()
        
                mutation = (ref_alle, alt_alle)
                
                if (ref_alle != alt_alle) and (ref_alle != "." and alt_alle != "."):
                    mutation_cnt+=1
                    if mutation not in mutation_list and (len(ref_alle)==1 and len(alt_alle)==1):
                        mutation_list[mutation]={}
                        mutation_list[mutation]["TYP"] = type_of_mutation(ref_alle, alt_alle)
                        mutation_list[mutation]["FQ"]=1
                        mutation_list[mutation]["LOCATION"]={}
                        mutation_list[mutation]["LOCATION"][chromosome_location]=1
                    elif mutation in mutation_list:
                        if chromosome_location not in mutation_list[mutation]["LOCATION"]:
                            mutation_list[mutation]["LOCATION"][chromosome_location]=1
                        else:
                            mutation_list[mutation]["LOCATION"][chromosome_location]+=1
                        mutation_list[mutation]["FQ"]+=1
    else:
        for row in vcf_parser:
            ref_alle=row[3]
            chromosome_location=row[0]
            alt_alle=row[4]
        
            mutation = (ref_alle, alt_alle) 
            if (ref_alle != alt_alle) and (ref_alle != "." and alt_alle != "."):
                mutation_cnt+=1
            
                if mutation not in mutation_list and (len(ref_alle)==1 and len(alt_alle)==1):
                    mutation_list[mutation]={}
                    mutation_list[mutation]["TYP"] = type_of_mutation(ref_alle, alt_alle)
                    mutation_list[mutation]["FQ"]=1
                    mutation_list[mutation]["LOCATION"]={}
                    mutation_list[mutation]["LOCATION"][chromosome_location]=1
                elif mutation in mutation_list:
                    if chromosome_location not in mutation_list[mutation]["LOCATION"]:
                        mutation_list[mutation]["LOCATION"][chromosome_location]=1
                    else:
                        mutation_list[mutation]["LOCATION"][chromosome_location]+=1
                    mutation_list[mutation]["FQ"]+=1


time_generated=datetime.now().strftime("%d/%m/%y %H:%M:%S")

for mutation,info in mutation_list.items():
    mutation_list[mutation]["LOCATION"] = dict(sorted(mutation_list[mutation]["LOCATION"].items(), key=lambda x: x[1], reverse=True))
    
df = pd.DataFrame.from_dict(mutation_list, orient="index")
if added_info == False:
    df = df.drop(columns=["LOCATION"])
df["FQ"] = round((df["FQ"]/mutation_cnt)*100, 2)

df = df.sort_values(by=["FQ"], ascending=False)


with open(output_file, "w") as f:
	print(f"#File name:{output_file}", file=f)
	print(f"#Number of mutations present: {mutation_cnt}",file=f)
	print("#Column Descriptions:", file=f)
	print("#TYP='Type of mutation'",file=f)
	print("#FQ='Frequenecy of mutation'",file=f)
	if added_info == True:
		print("#LOCATION=List of chromosomes mutation was seen in",file=f)
	print(f"#Date file was generated:{time_generated}",file=f)

fig = df.plot(kind='barh', color='r', title=(f"{output_file}").split('/')[1].split('.')[0]).get_figure()
fig.savefig(f"{output_file.split('.')[0]}.png")
df.to_csv(output_file, sep='\t', mode='a')
