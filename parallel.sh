parallel --dry-run -j 4 -k echo \"{} "<-- a number"\" ::: `seq 1 5`

#cat seq
#1
#2
#3
#4
#5
cat seq | parallel --dry-run -j 4 -k echo \"{} "<-- a number"\"

##real example, get the coverage
ls -1 bam/20230420_file/ | sed 's/.bam//g' > file
cat file|sed  s/\ //g | parallel "samtools depth -a bam/{}.bam > bam_coverage/{}.coverage"

cat file | parallel -j 4 "bwa mem  ref/hg38.fa fastq/{}_R1_001.fastq.gz fastq/{}_R2_001.fastq.gz | samtools view -bh -F2048 - | samtools sort > bam/{}.bam"


## bwa index
ls -1 | sed 's/.final.fasta//g' > file
cat file| parallel "bwa index DF825_SO11626/{}.final.fasta"


## one line of code do one action to all the files in one folder
find . -name "*.fastq" | xargs basename -s ".fastq" | xargs -I{} fastq_stat --in {}.fastq --out ../summaries/{}.txt


## xargs, find, basename and parallel
find . -name "*.fastq" | xargs basename -s ".fastq" | xargs -P 6 -I{} fastq_stat --in {}.fastq --out ../summaries/{}.txt 

## use awk to create scripts for automation
cat sraids.txt | awk '{print "fastq-dump -X 2000 --split-files -0 sra" $1}' > get_data.sh
bash get-data.sh
