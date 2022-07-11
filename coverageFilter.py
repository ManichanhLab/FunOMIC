#!/usr/bin/env python

import pandas as pd
import argparse
import re

#### Import from arguments
parser = argparse.ArgumentParser(description="get percentage of coverage of each reads")

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

option = parser.parse_args()

q30reads = option.input_path
coverage80 = option.output_path

q30sam=pd.read_csv(q30reads,delimiter='\t', error_bad_lines=False,names=["QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","BT1","BT2","BT3","BT4","BT5","BT6","BT7","BT8","BT9","BT10","BT11"])
l=[]
len=[]
for cigar in q30sam['CIGAR']:
	countM=0
	countAll=0
	match=re.findall(r"(\d+)M",cigar)
	all=re.findall(r"(\d+)",cigar)
	for m in match:
		countM+=int(m)
	for a in all:
		countAll+=int(a)
	coverage=countM/countAll*100
	rlen=countAll
	len.append(rlen)
	l.append(coverage)
q30sam['coverage']=l
q30sam['RLEN']=len
q30cov80sam=q30sam[(q30sam['coverage']>=80) & (q30sam['RLEN']>=60)]
qnames=q30cov80sam['QNAME']
qnames.to_csv(coverage80,sep='\t',index=False)
