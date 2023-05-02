parallel --dry-run -j 4 -k echo \"{} "<-- a number"\" ::: `seq 1 5`

#cat seq
#1
#2
#3
#4
$5
cat seq | parallel --dry-run -j 4 -k echo \"{} "<-- a number"\"
