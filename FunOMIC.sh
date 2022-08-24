#!/usr/bin/env bash


################################################################################
Help()
{
   # Display Help
   echo
   echo "FunOMIC: mycobiome taxonomic and functional profiling pipeline"
   echo
   echo "Syntax: FunOMIC -1 -2 -p -o -a -b -c -t"
   echo "options:"
   echo "1     path of paired end read1 (we recommend to apply kneaddata or other quality control first)"
   echo "2     path of paired end read2 (we recommend to apply kneaddata or other quality control first)"
   echo "p     output prefix"
   echo "o     path of output directory"
   echo "a     path of UHGG bacterial database directory"
   echo "b     path of FunOMIC-T databases directory"
   echo "c     path of FunOMIC-P databases directory"
   echo "t     number of threads"
   echo
}



######################################
########### SET VARIABLES ############
######################################

while getopts "1:2:p:o:a:b:c:t:h" flag
do
    case "$flag" in
        1) read1=$OPTARG;; # path of paired end read1 (we recommend to apply kneaddata or other quality control first)
        2) read2=$OPTARG;; # path of paired end read2 (we recommend to apply kneaddata or other quality control first)
        p) prefix=$OPTARG;; #output prefix
        o) outdir=$OPTARG;; #path of output directory
        a) bactdb=$OPTARG;; #path of UHGG bacterial database
        b) taxadb=$OPTARG;; #path of FunOMIC-T databases
        c) protdb=$OPTARG;; #path of FunOMIC-P databases
        t) threads=$OPTARG;; #number of threads
        h) Help
           exit;; #display help
    esac
done



# ##############################################
# ######### DECONTAMINATION WITH UHGG! #########
# ##############################################

printf "Start removing bacterial reads for ${prefix}\n"

if [[ ! -f  $outdir ]]
    then   # When the pipeline has been executed properly and the script recognizes the output
      mkdir $outdir
      if [[ ! -f  $outdir/tmp ]]
         then mkdir $outdir/tmp
      fi 
      if [[ ! -f $outdir/taxonomic_profiling ]]
         then mkdir $outdir/taxonomic_profiling
   fi
   if [[ ! -f $outdir/functional_profiling ]]
         then mkdir $outdir/functional_profiling
   fi 
fi &>/dev/null


bowtie2 -p $threads -x $bactdb/uhgg \
-1 $read1 \
-2 $read2 \
-S $outdir/tmp/${prefix}_Bact.sam 2> ${outdir}/bact_decontam.log

samtools view -b -f 4 $outdir/tmp/${prefix}_Bact.sam > $outdir/tmp/${prefix}_noBact.bam 

samtools sort -n $outdir/tmp/${prefix}_noBact.bam -o $outdir/tmp/${prefix}_noBact_sorted.bam &>/dev/null
samtools fastq -@ 8 $outdir/tmp/${prefix}_noBact_sorted.bam \
   -1 $outdir/${prefix}_noBact_1.fastq.gz \
   -2 $outdir/${prefix}_noBact_2.fastq.gz \
   -0 /dev/null -s /dev/null -n &>/dev/null



# # ##############################################
# # ########### TAXONOMICAL PROFILING! ###########
# # ##############################################

printf "Start taxonomic annotation for ${prefix}\n"

bowtie2 -p $threads -x $taxadb/FunOMIC.T.v1 \
-1 $outdir/${prefix}_noBact_1.fastq.gz \
-2 $outdir/${prefix}_noBact_2.fastq.gz \
-S $outdir/tmp/${prefix}.sam 2> $outdir/taxonomic_profiling/log

# filter hits with q-score over 30 and coverage over 80
samtools view -Sq 30 $outdir/tmp/${prefix}.sam > $outdir/tmp/${prefix}.30.sam 
samtools view -bSq 30 $outdir/tmp/${prefix}.sam > $outdir/tmp/${prefix}.30.bam
samtools sort -o $outdir/tmp/${prefix}.30.sorted.bam $outdir/tmp/${prefix}.30.bam  &>/dev/null
coverageFilter.py -i $outdir/tmp/${prefix}.30.sam -o $outdir/tmp/${prefix}.30.c80.list &>/dev/null
samtools view $outdir/tmp/${prefix}.30.sorted.bam | grep -f $outdir/tmp/${prefix}.30.c80.list > $outdir/tmp/${prefix}.30.c80.sam
samtools view -bt $taxadb/FunOMIC.T.v1.fasta.fai -o $outdir/tmp/${prefix}.30.c80.bam $outdir/tmp/${prefix}.30.c80.sam &>/dev/null
samtools index $outdir/tmp/${prefix}.30.c80.bam &>/dev/null
samtools idxstats --threads $threads $outdir/tmp/${prefix}.30.c80.bam > $outdir/tmp/${prefix}.idxstats.txt


awk '$3 > 0' $outdir/tmp/${prefix}.idxstats.txt | sed 's/|/\t/g' | awk '{print $1,$4}' > $outdir/tmp/${prefix}"_assembly_hits.txt"
buglist.py -i $outdir/tmp/${prefix}"_assembly_hits.txt" -t $taxadb/taxonomy.txt -o $outdir"/taxonomic_profiling"/${prefix}"_buglist_stratified.txt" &>/dev/null
  

#### rm temporary files
rm -r $outdir/tmp/

##############################################
########### FUNCTIONAL PROFILING! ############
##############################################


printf "Starting functional annotation for ${prefix}\n"

cd ${outdir}/functional_profiling
mkdir $outdir/functional_profiling/tmp/

flash2 $outdir/${prefix}_noBact_1.fastq.gz $outdir/${prefix}_noBact_2.fastq.gz -o tmp/${prefix} -q -t $threads &>/dev/null
cat $outdir/functional_profiling/tmp/${prefix}.extendedFrags.fastq $outdir/functional_profiling/tmp/${prefix}.notCombined_1.fastq $outdir/functional_profiling/tmp/${prefix}.notCombined_2.fastq > $outdir/functional_profiling/joined.fastq


printf "merged clean reads stored in: $outdir/functional_profiling/joined.fastq \n"

  diamond blastx -b5 -c1 --query-cover 95 \
  --id 99 -d $protdb"/FunOMIC.P.v1.dmnd" \
  -q $outdir"/functional_profiling"/joined.fastq --outfmt 6 \
  -o $outdir/"functional_profiling"/${prefix}.func.out 2> $outdir/functional_profiling/log

# Filter alignment hits: coverage>80%; query length>100bp; sort by perc_identity, length, bit score, remove redundancy
awk '$4 >= 20' $outdir/"functional_profiling"/${prefix}.func.out | sort -k3,3nr -k4,4nr -k12,12nr | awk '!a[$1]++' > $outdir/"functional_profiling"/${prefix}.func.filtered.out

    if [[ -f $outdir/"functional_profiling"/${prefix}.func.filtered.out ]]
    then   # When the pipeline has been executed properly and the script recognizes the output
         printf "removing secondary files files for ${p}\n\n"
         rm $outdir/"functional_profiling"/${prefix}.func.out
         rm $outdir/functional_profiling/*.fastq
         rm -r $outdir/functional_profiling/${prefix}.hist*
    else echo "No output generated"
    fi &>/dev/null

functional_profiling.py -i $outdir/"functional_profiling"/${prefix}.func.filtered.out\
 -c $protdb/id_to_clean.txt -a $protdb/nju_id_protein_ec.txt \
 -o1 $outdir/"functional_profiling"/${prefix}"_protein.csv"\
 -o2 $outdir/"functional_profiling"/${prefix}"_EC.csv" &>/dev/null

get_pwy_abundance.py -i $outdir/"functional_profiling"/${prefix}.func.filtered.out\
 -c $protdb/id_to_clean.txt -a $protdb/jgi_annotation_Nov2021.tab \
 -o1 $outdir/"functional_profiling"/${prefix}"_pwy_abundance.csv"\
 -o2 $outdir/"functional_profiling"/${prefix}"_pwyClass_abundance.csv"\
 -o3 $outdir/"functional_profiling"/${prefix}"_pwyType_abundance.csv"\
 -o4 $outdir/"functional_profiling"/${prefix}"_annotation.csv" &>/dev/null
