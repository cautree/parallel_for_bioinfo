
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



## sed and regex
### (chr[^:]+) is to get chr1 etc, [^:]+ means get any char or number but that is not :
### :([0-9]+)-  get the numbers, which are in btw : and -
### \1\t\2\t\3 means the three groups are sepearated by tabs
### original format is sed -E 's/a/b/g' 
echo "chr1:28427874-28425431" | sed -E 's/^(chr[^:]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/'

## using[] in sed, change ":" and "-" in to "\t"
echo "chr1:28427874-28425431" | sed 's/[:-]/\t/g'

## using sed twice, with pipe or without pipes
echo "chr1:28427874-28425431" | sed 's/:/\t/' | sed 's/-/\t/' 
echo "chr1:28427874-28425431" | sed  -e 's/:/\t/' -e 's/-/\t/'

## using tr
echo "chr1:28427874-28425431" | tr ':-' '\t'


## subshell, the end result is both commdn will be changed to step, if no (), only one will be changed
(echo "this command"; echo "that command") | sed 's/command/step/'


## awk filter rows, two conditions connected by &&, ~ here is search patten of /chr1/
awk '$1 ~ /chr1/ && $3 - $2 > 10' example.bed


## wak BEGIN, END, notice {} is connected by ;
awk 'BEGIN{ s = 0 }; { s += ($3-$2) }; END{ print "mean: " s/NR };'


## filter to see if the row contains a specific val, if yes, that row is returned
awk '/chr3/' example.bed 


## awk associative array
awk '/Lypla1/ { feature[$3] += 1 }; END { for (k in feature) print k "\t" feature[k] }' ../data/Mus_musculus.GRCm38.75_chr1.gtf

## linux way to do the same thing as above
grep "Lypla1" ../data/Mus_musculus.GRCm38.75_chr1.gtf | cut -f 3 | sort | uniq -c


## Commands:"create,"extract,"gzip,"file,"list,"verbose"
tar cvf $myfile.tar. $sample1.fq. $sample2.fq


## regex
####AT{3}, A followed by 3 T
####(AT}{3}, 3 AT in a row
####(ATG)+C{2}, one or more ATG followed by 2 C


## complicated awk, file name is split.awk
$3 == "gene" {
# split the 9th column by ;
split($9, x, ";")
split(x[1], y, " ")
## remove the double quotes around the gene name
name = y[2]
## substitue of " with empty space
gsbu("\"", "", name)
print $3, name, $5-$4 +1

}

cat NC.gff | awk -f split.awk | head -5
