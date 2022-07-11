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
						help = "Output path of the fungal pathway profiling.\n",
						required = True)
parser.add_argument('-o2', '--output2',
						dest = "output_path2",
						action = "store",
						default = "",
						help = "Output path of the fungal pathway_class profiling.\n",
						required = True)
parser.add_argument('-o3', '--output3',
						dest = "output_path3",
						action = "store",
						default = "",
						help = "Output path of the fungal pathway_type profiling.\n",
						required = True)
parser.add_argument('-o4', '--output4',
						dest = "output_path4",
						action = "store",
						default = "",
						help = "Output path of the fungal full annotation.\n",
						required = True)

option=parser.parse_args()

sample_to_process = option.input_path
output1=option.output_path1
output2=option.output_path2
output3=option.output_path3
output4=option.output_path4

id_catalog=pd.read_csv("/mnt/synology/RAW_DATA/FUNGAL_GENOMES/id_to_clean.txt",delimiter="\t",names=["pid","rid"])
jgi_ann=pd.read_csv("/mnt/synology/RAW_DATA/FUNGAL_GENOMES/JGI_CDS/version_Nov2021/annotation_all.tab",delimiter='\t',names=["proteinId","ecNum","definition","catalyticActivity","cofactors","associatedDiseases","pathway","pathway_class","pathway_type"])
blastx_out=pd.read_csv(sample_to_process,delimiter='\t',names=["qid","rid","id","lenght","mismatch","gapopne","qstart","qend","sstart","send","evalue","bscore"])
if len(blastx_out) >=1:
	raw_counts=blastx_out['rid'].value_counts()
	raw_counts=blastx_out['rid'].value_counts().rename_axis("rid").reset_index(name='counts')
	raw_counts_clean=pd.merge(raw_counts,id_catalog)
	column_names = ["proteinId","ecNum","definition","catalyticActivity","cofactors","associatedDiseases","pathway","pathway_class","pathway_type"]
	pwy_counts=pd.DataFrame(columns=column_names, dtype=object)
	for i in range(len(raw_counts_clean)):
		inter=jgi_ann.loc[jgi_ann['proteinId']==raw_counts_clean.iloc[i]['pid']]
		inter=inter.assign(counts=raw_counts_clean.iloc[i]['counts']/len(inter))
		pwy_counts=pwy_counts.append(inter,ignore_index = True)
	pwy=pwy_counts.groupby(["pathway"])["counts"].sum().reset_index()
	pwy_class=pwy_counts.groupby(["pathway_class"])["counts"].sum().reset_index()
	pwy_type=pwy_counts.groupby(["pathway_type"])["counts"].sum().reset_index()
	pwy.to_csv(output1,index=False,sep='\t')
	pwy_class.to_csv(output2,index=False,sep='\t')
	pwy_type.to_csv(output3,index=False,sep='\t')
	pwy_counts.to_csv(output4,index=False,sep='\t')
else:
	print('no protein found\n')
	pass