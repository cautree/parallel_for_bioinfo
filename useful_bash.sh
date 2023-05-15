
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
aws s3 ls s3://s3_bucket/ | grep "[[:space:]]ecoli.fa"| awk '{$1=$1}1' OFS="," | cut -f 4 -d ","  


## awk filter based on one column, very slow so better use bed file during using samtools depth
cat coverage_E10_E11.coverage | awk '/chr22/' > chr22.coverage
samtools depth -a 230426_E09.md.bam 230426_E10.md.bam 230426_E11.md.bam -b chr22.bed > chr22_E091011.coverage


## a quick summary of dataframe
summarytools::view(dfSummary(as.data.frame(mtcars), 
                             style = 'grid',
                             max.distinct.values = 10, 
                             plain.ascii =   F, 
                             valid.col = F, 
                             headings = F), method = "render")

## corrplot
cormat <- cor(X, use="complete.obs")
corrplot(cormat, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"),
         title = "",
         addCoef.col = "black",
         tl.cex=.5, number.cex=.5)
         
         
## read sample info from a file and do something to all the files, from bio carpentry
#!/bin/bash

set -e
set -u
set -o pipefail

# specify the input samples file, where the third column is the path to each sample FASTQ file
sample_info=samples.txt

# create a Bash array from the third column of $sample_info
sample_files=($(cut -f 3 "$sample_info"))

for fastq_file in ${sample_files[@]} 
do
    # strip .fastq from each file, and add suffix "-stats.txt" to create an output filename
    results_file="$(basename $fastq_file .fastq)-stats.txt"
    
    # run fastq_stat on a file, writing results to the filename we've # above
    fastq_stat $fastq_file > stats/$results_file
done
