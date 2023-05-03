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
cat file|sed  s/\ //g | parallel "samtools depth bam/{}.bam > bam_coverage/{}.coverage"


