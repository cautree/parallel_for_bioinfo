
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

## select fasta from reference file based on bed file
bedtools getfasta -fi ../lambda_seq/lambda.fa -bed lambda_amplicon_4.bed 



#contain both hunger & and on the same line
awk '/hunger/&&/and/' resources/gutenberg/pg74.txt

#contain either hunger or hunger
awk '/hunger|together/' resources/gutenberg/pg74.txt

## counts lines has digits
cat test.txt | awk -e ' /[[:digit:]]/ { acount +=1 }  END { print acount}'

## get everything except the first column
cat test.txt | cut -f2- -d " " 

## bcftools view options, -H; supress header, -i: filter, -c1: at least one allele difference from ref, -V: exclude 
cat 220623-PB184-diploid-8x.vcf | bcftools view -H -i "QUAL>20" -c1 -V indels - | cut -f 2,4,5,6

## grep space, get ecoli.fa, not kk_ecoli.fa
aws s3 ls s3://s3_bucket | grep "[[:space:]]ecoli.fa"| awk '{$1=$1}1' OFS="," | cut -f 4 -d ","  


