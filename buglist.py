#!/usr/bin/env python

import pandas as pd
import argparse

#### Import from arguments
parser = argparse.ArgumentParser(description="process number of reads to buglist")

parser.add_argument('-i', '--input',
			dest = "input_path",
			action = "store",
			default = None,
			help = "Sample to import in this script \n")

parser.add_argument('-o', '--output',
			dest = "output_path",
			action = "store",
			default = None,
			help = "Sample to output in this script \n")

parser.add_argument('-t', '--taxa',
			dest = "taxanomy_path",
			action = "store",
			default = None,
			help = "Path of the taxanomy file \n")

option = parser.parse_args()

sample_to_process = option.input_path
output_path = option.output_path
taxa_path = option.taxanomy_path

bt_counts=pd.read_csv(sample_to_process, delimiter=" ",names=["Id","Abundance"])
if len(bt_counts) >= 1:
	raw_counts=bt_counts.groupby(bt_counts.columns[0])[bt_counts.columns[-1]].median().reset_index()
	# taxa=pd.read_csv('/mnt/storage5TB/zixuan/all_busco_taxa-humannFormat.txt', delimiter="\t",names=["Id","taxa"]) ## taxa of version1
	taxa=pd.read_csv(taxa_path, delimiter="\t",names=["Id","taxa"],index_col=False)
	raw_counts.columns=["Id","Abundance"]
	taxa.columns=["Id","taxa"]
	buglist=pd.merge(taxa,raw_counts)
	column_names=['kingdom',"phylum",'class','order','family','genus','species','strain']
	stratified_counts=pd.DataFrame(buglist.taxa.str.split('|').tolist(),columns=column_names, dtype=object)
	results=pd.concat([buglist,stratified_counts],axis=1)
	kingdom_counts=results.groupby(["kingdom"]).sum().reset_index()
	phylum_counts=results.groupby(["kingdom","phylum"]).sum().reset_index()
	phylum_counts['combined']=phylum_counts['kingdom'].astype(str)+'|'+phylum_counts['phylum']
	class_counts=results.groupby(["kingdom","phylum","class"]).sum().reset_index()
	class_counts['combined']=class_counts['kingdom'].astype(str)+'|'+class_counts['phylum']+'|'+class_counts['class']
	order_counts=results.groupby(["kingdom","phylum","class","order"]).sum().reset_index()
	order_counts['combined']=order_counts['kingdom'].astype(str)+'|'+order_counts['phylum']+'|'+order_counts['class']+'|'+order_counts['order']
	family_counts=results.groupby(["kingdom","phylum","class","order","family"]).sum().reset_index()
	family_counts['combined']=family_counts['kingdom'].astype(str)+'|'+family_counts['phylum']+'|'+family_counts['class']+'|'+family_counts['order']+'|'+family_counts['family']
	genus_counts=results.groupby(["kingdom","phylum","class","order","family","genus"]).sum().reset_index()
	genus_counts['combined']=genus_counts['kingdom'].astype(str)+'|'+genus_counts['phylum']+'|'+genus_counts['class']+'|'+genus_counts['order']+'|'+genus_counts['family']+'|'+genus_counts['genus']
	species_counts=results.groupby(["kingdom","phylum","class","order","family","genus","species"]).sum().reset_index()
	species_counts['combined']=species_counts['kingdom'].astype(str)+'|'+species_counts['phylum']+'|'+species_counts['class']+'|'+species_counts['order']+'|'+species_counts['family']+'|'+species_counts['genus']+'|'+species_counts['species']
	strain_counts=results.groupby(["kingdom","phylum","class","order","family","genus","species","strain"]).sum().reset_index()
	strain_counts['combined']=strain_counts['kingdom'].astype(str)+'|'+strain_counts['phylum']+'|'+strain_counts['class']+'|'+strain_counts['order']+'|'+strain_counts['family']+'|'+strain_counts['genus']+'|'+strain_counts['species']+'|'+strain_counts['strain']
	#### l8
	#new=pd.concat([kingdom_counts.kingdom.append(phylum_counts.combined).append(class_counts.combined).append(order_counts.combined).append(family_counts.combined).append(genus_counts.combined).append(species_counts.combined).append(strain_counts.combined),kingdom_counts.Abundance.append(phylum_counts.Abundance).append(class_counts.Abundance).append(order_counts.Abundance).append(family_counts.Abundance).append(genus_counts.Abundance).append(species_counts.Abundance).append(strain_counts.Abundance)],axis=1)
	#### l7
	new=pd.concat([kingdom_counts.kingdom.append(phylum_counts.combined).append(class_counts.combined).append(order_counts.combined).append(family_counts.combined).append(genus_counts.combined).append(species_counts.combined),kingdom_counts.Abundance.append(phylum_counts.Abundance).append(class_counts.Abundance).append(order_counts.Abundance).append(family_counts.Abundance).append(genus_counts.Abundance).append(species_counts.Abundance)],axis=1)
	new.columns=["taxa","Abundance"]
	new.to_csv(output_path, sep='\t',index=False)
else:
	print('no fungi found\n')
	pass
