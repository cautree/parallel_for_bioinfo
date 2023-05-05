parallel --dry-run -j 4 -k echo \"{} "<-- a number"\" ::: `seq 1 5`

#cat seq
#1
#2
#3
#4
#5
cat seq | parallel --dry-run -j 4 -k echo \"{} "<-- a number"\"

##real example, get the coverage
ls -1 bam/20230420_file/ | tr ".bam" " " > file
cat file|sed  s/\ //g | parallel "samtools depth -a bam/{}.bam > bam_coverage/{}.coverage"

cat file | parallel -j 4 "bwa mem  ref/hg38.fa fastq/{}_R1_001.fastq.gz fastq/{}_R2_001.fastq.gz | samtools view -bh -F2048 - | samtools sort > bam/{}.bam"


