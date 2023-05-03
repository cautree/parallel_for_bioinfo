
## this is for batch download stuff from s3 bucket, first create a temp file, then loop the temp file

aws s3 ls s3://s3_bucket/ | grep ecoli_REL606 | grep -v _ecoli_ |\
    awk '{$1=$1}1' OFS="," |\
    cut -f 4 -d "," > temp


#!/bin/bash
set -e
while read line;
do
  aws s3 cp s3://s3_bucket/$line  ref/
done <temp


## regex for well location of 96 well plate
cat info_email | tail -10 | grep -o -E  "[A-H]{1}[0-9]{2}" > well.txt

## show the files that match the patten in a file, then move that file to a fold called keep
ls -1 | grep -f ../../meta_data/20230322_20ng.txt > file
cat file | parallel --eta --verbose "mv {} keep"


## counts total reads of all the R1 file  in one folder
zcat *R1* | paste - - - - | wc -l 


## get the well position for the fastq files
aws s3 ls s3://s3_bucket/220616-UDI-PlexC_FASTQ/ |  awk '{$1=$1}1' OFS="," | cut -f 4 -d "," | cut -f 2  -d "_" | sort | uniq > well


## find the position of forwrd and reverse primer 
cat lambda_seq/lambda.fa | dreg -filter -pattern 'GTGCAGCCGGTCTTAAAC' | grep -E "regex:" > frag_250_f.txt
echo GCCTCCTGGGCAGTC | tr ACGTacgt TGCAtgca | rev
cat lambda-seq/lambda.fa | dreg -filter -pattern 'GACTGCCCAGGAGGC' | grep -E "regex:" > frag_250_r.txt


## create a correct formated bed file
## [according to https://www.biostars.org/p/404859/]
awk 'OFS=" " {print $1"\t", $2"\t", $3"\t"}' regions_file.bed | tr -d " " > outputs/lambda_amplicon_4.bed



