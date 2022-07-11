#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Run fungal functional profiling, Input is the diamond blastx output.\n")

parser.add_argument('-i', '--input',
                        dest = "input_path",
                        action = "store",
                        default = None,
                        help = "Input blastx output path.\n",
                        required = True)


parser.add_argument('-o1', '--output1',
                        dest = "output_path1",
                        action = "store",
                        default = "",
                        help = "Output path of the fungal functional profiling.\n",
                        required = True)
parser.add_argument('-o2', '--output2',
                        dest = "output_path2",
                        action = "store",
                        default = "",
                        help = "Output path of the fungal EC profiling.\n",
                        required = True)

option=parser.parse_args()

sample_to_process = option.input_path
output1=option.output_path1
output2=option.output_path2
id_catalog=pd.read_csv("/mnt/synology/RAW_DATA/FUNGAL_GENOMES/id_to_clean.txt",delimiter="\t",names=["pid","rid"])
annotation=pd.read_csv("/mnt/synology/RAW_DATA/FUNGAL_GENOMES/nju_id_protein_ec.txt",delimiter='\t',names=["pid","pname","ECnumber"])
blastx_out=pd.read_csv(sample_to_process,delimiter='\t',names=["qid","rid","id","lenght","mismatch","gapopne","qstart","qend","sstart","send","evalue","bscore"])
if len(blastx_out) >=1:
	raw_counts=blastx_out['rid'].value_counts()
	raw_counts=blastx_out['rid'].value_counts().rename_axis("rid").reset_index(name='counts')
	raw_counts_clean=pd.merge(raw_counts,id_catalog)
	protein_counts=pd.merge(raw_counts_clean,annotation)
	protein_counts["counts"]=protein_counts["counts"].astype(str).astype(int)
	ec=protein_counts.groupby(["ECnumber"])["counts"].sum().reset_index()
	pname=protein_counts.groupby(["pname"])["counts"].sum().reset_index()
	pname.to_csv(output1,index=False,sep='\t')
	ec.to_csv(output2,index=False,sep='\t')

else:
	print('no protein found\n')
	pass