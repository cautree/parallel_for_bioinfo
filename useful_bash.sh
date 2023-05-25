
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



## good to have alias for awk
alias awk = "awk -F '\t' =v OFS='\t' "


## create resuable function
i() { (head -n 2; tail -n 2) < "$1" | column -t}
i a.txt




## from https://mp.weixin.qq.com/s?__biz=MzIxMjQxMDYxNA==&mid=2247484065&idx=1&sn=580783ce880008b2809b65d060460d6f&chksm=9747cb38a030422ecea3dcfaed25b7b7d15bf6119e66e919baeab5a601826eb3e5cb3312868a&token=1106388583&lang=zh_CN&scene=21#wechat_redirect
## 作者：简佐义（jianzuoyi@qq.com）
ls `pwd`/file
cat R1.fq.gz R2.fq.gz  # by row
paste -d ' ' file1 file2 #by column
grep '^>' test.fa | cut -c 2-  ## -c is by character position
less -SN file                # 显示文件的行号，并且截断太长的行
tail -n +2 file	             # 跳过第1行，显示从第2行开始的所有行，可用于跳过文件的标题行
touch {file1,file2,file3}	 # 同时创建3个文件
tar czvf file.tar.gz files	 # 打包并压缩
tar xvf file.tar.gz			 # 解包，解压缩
gzip file					 # 压缩
gunzip file.gz				 # 解压
chmod a+x file	             # 增加[所有人]可执行权限
chmod a-x file	             # 取消[所有人]可执行权限
sort -k2,2 -k3,3 file	     # 先按第2列排序，第2列相同，再按第3列排序
sort -k2,2n file		     # 按第2列排序，且第2列是数字，升序
sort -k2,2nr file		     # 按第2列排序，且第2列是数字，降序
sort -u file			     # 先排序文件，然后去除相邻的重复行，只保留一条记录
sort file | uniq		     # 去除相信的重复行，只保留一条记录，相当于： sort -u file
sort a b | uniq			     # 并集
sort a b | uniq -d > c	     # 交集
sort a c | uniq -u 		     # 补集
df -h		                 # 查看磁盘使用情况，-h表示以人类可读的方式显示容量大小
du -sh		                 # 查看当前目录使用了多少磁盘空间
du -sh *	                 # 查看当前目录下各文件或文件夹使用的磁盘空间
free -h		                 # 查看内存使用情况
top -c		                 # 查看CPU，内存的使用情况
ps aut		                 # 查看后台任务运行情况，第2列是任务的PID号
kill -9 PID                  # 删除编号为PID的任务
nohup ./run.sh &> run.sh.o & # 远程任务管理
##########################################################################
./run.sh	                 # 假如任务是直接这样开始跑的
ctrl + z	                 # 按ctrl + z，将任务放到后台
jobs                         # 输入jobs命令，回车，可以看到任务是暂停的： [1]+  Stopped(SIGTSTP)        bash run.sh
bg                           # 让后台暂停的任务开始运行
jobs                         # 再次运行jobs，可以看到任务已经跑起来了：   [1]+  Running                 bash run.sh &
disown -r 	                 # 从当前shell中移除运行中的作业，至此，可以关掉终端回家了
##########################################################################
./run.sh > run.sh.o		     # 标准输出到run.sh.o日志文件
./run.sh 2> run.sh.e	     # 标准错误输出到run.sh.e错误日志文件
./run.sh &> run.sh.log	     # 标准输出和标准错误都输出到定一个文件
./run.sh &> /dev/null	     # 丢弃标准输出和标准错误信息
###################################################################
#command < file			     # 将file的内容作为command的输入
#command << END			     # 从标准输入（键盘）中读取数据，直到遇到分界符END时停止（分界符用户可以自定义）
#command <file1 > file2	     # 将file1作为command的输入，并将处理结果输出到file2
### one application
whileread line
do
    do something
done < file.txt > result.txt
##################################################################
find -name file					# 在当前目录查找名为file的文件
find dir/ -name file			# 在dir/目录下查找名为file的文件
find dir/ -name '*file*'# 在dir/目录下查找包含file关键词的文件，-name参数支持正则表达式
find dir/ -name file -delete	# 查找文件并删除

#####################################################################
cat file | xargs		                                            # 将file的内容显示成一行
find /ifs/result -name '*.fq.gz' | xargs -n1 -I{} cp {} /ifs/data/	# 查找fq.gz文件并复制到/ifs/data目录下
find /ifs/result -name '*.fq.gz' | xargs tar czvf all.fq.gz			# 查找fq.gz文件并打包在一起
find . -type f -name '*.log' -print0 | xargs -0 rm -f				# 当rm文件过多时，可以这样删除
find . -type f -name '*.py' -print0 | xargs -0 wc -l				# 统计一个目录中所有python文件的行数
#####################################################################
find *.fq | parallel -j 12 "fastqc {} --outdir ."# 同时执行12个Fastqc任务
#####################################################################
grep pattern files			                # 搜索文件中包含pattern的行
grep -v pattern files		                # 搜索文件中不包含pattern的行
grep -f pattern.txt files	                # 搜索的pattern来自于文件中
grep -i pattern files		                # 不区分大小写。默认搜索是区分大小写的
grep -w pattern files		                # 只匹配整个单词，而不是字符串的一部分（如搜索hello，不会匹配到helloworld）
grep -n pattern files		                # 显示行号信息
grep -c pattern files		                # 显示匹配的行数
grep -l pattern files		                # 只显示匹配的文件名
grep -L pattern files		                # 显示不匹配的文件名
grep -C number pattern files                # 额外显示匹配行的上下[number]行
grep pattern1 | grep pattern2 files         # 显示既匹配pattern1，又匹配pattern2的行
grep -E "pattern1|pattern2" files	        # 显示匹配pattern1或者pattern2的行, grep -E相当于egrep

# 用于搜索的特殊字符
^: 表示行前
$: 表示行尾

grep '^#' result.vcf		                # 显示VCF文件的表头信息
grep '^hello$' files		                # 显示只包含hello的行
grep -v '^\s*$' file		                # 删除空白行

########################################################################
##sed
#d：删除该行
#p：打印该行
#i：在行的前面插入新行
#a：在行的后面插入新行
#r：读取指定文件的内容。
#w：写入指定文件。
sed -n '10p' file		                # 显示第10行
sed -n '10,20p' file	                # 显示第10到20之间的行
sed -n '/pattern/p' file                # 显示含有pattern的行
sed -n '/pattern1/,/pattern2/p' file    # 显示patter1与pattern2之间的行

sed '10d' file			                # 删除第10行
sed '10,20d' file		                # 删除第10到20之间的行
sed '/pattern/d'                        # 删除匹配pattern的行
sed '/^\s*$/d' file		                # 删除空白行
sed 's/^\s*//' file		                # 删除行前的空白：空格，制表符
sed 's/\s*$//' file		                # 删除行尾的空白：空格，制表符
sed 's/^\s*//;s/\s*$//' file            # 删除行首和行尾的空白：空格，制表符

sed 's/AA/BB/' file		                # 将文件中的AA替换成BB，只替换一行中第一次出现的AA，替换后的结果输出到屏幕
sed 's/AA/BB/g' file	                # 将文件中的所有AA都替换成BB，替换后的结果输出到屏幕
sed -i 's/AA/BB/g' file                 # 将文件中的所有AA都替换成BB，直接更改文件的内容
sed '/CC/s/AA/BB/g' file                # 只替换那些含有CC的行
sed 's/pattern/&XXXX/' file	            # 在pattern之后加上XXXX。&表示之前被匹配的内容
sed 's/pattern.*/&XXXX' file            # 在匹配pattern的行尾加上XXXX。pattern.*表示包含pattern的整行内容

sed -n '1~4s/^@/>/p;2~4p' file.fq > file.fa	# Fastq文件转Fasta文件
sed -n '2~4p' file.fq		            # 提取Fastq文件的序列

sed 'y/ABC/XYZ/' file	                # 将ABC逐字替换成XYZ

sed '1i\hello' file		                # 在第1行前面插入一行，内容为hello，通常用来为文件增加标题
sed '1a\hello' file		                # 在第1行后面插入一行，内容为hello
sed '1r file2' file1	                # 在第1行后面读入file2的内容
sed '/pattern/w file2' file1 # 将匹配的行写入file2中

############################### awk
seq 20 | xargs -n5 > file
# cat file
1 2 3 4 5
6 7 8 9 10
11 12 13 14 15
16 17 18 19 20

awk '$5 ~ /10/' file
awk '$5 ~ "10"' file
awk '$5 ~ 10' file



#####
## getline							# 读取下一条记录到 $0，更新NF，NR和FNR
seq 10 | awk '{print $0;getline}'      # 显示奇数行
seq 10 | awk '{getline; print $0}'     # 显示偶数行
seq 10 | awk '{getline tmp; print tmp; print $0}'   # 奇偶行对调
# fastq转换成fasta
awk '{getline seq; getline comment; getline quality; sub("@", ">", $0); print $0"\n"seq}' file
awk -v n=$number '{print n, $0}' file
awk '{print $0}' file	# 打印整行
awk '{print $1}' file	# 打印第一列
awk '{print $2}' file	# 打印第二列
awk '{print $NF}' file	# 打印最后一列
awk '{print $(NF-1)}' file#打印倒数第二列
awk -F ';' -v OFS='\t''{print $1,$2,$NF}' file	# 读入的文件以逗号;分隔列，打印第1列，第2列和最后一列，并且打印时以制表符作为列的分隔符
awk '$2 > 100' file		# 打印第2列大于100的行
awk 'NR>1 && NR<4' file # 打印第2~3行

awk '/EGFR/' file		# 打印含有EGFR的行，相当于grep EGFR file
awk '$1 ~ /EGFR/' file	# 打印第1列含有EGFR的列


# 按指定列去除重复行
# cat file
#1 2 3 4 5
#6 2 8 9 10
#11 12 13 14 15
#16 17 18 19 20
awk '!a[$2]++' file		# 第二列出现两次2，只保留第一次出现的那一行，结果如下：
#1 2 3 4 5
#11 12 13 14 15
#16 17 18 19 20

awk '{sum+=$1} END {print sum}' file	# 累加文件的第一列
awk '{sum+=$1} END {print sum/NR}' file	# 求第一列的平均数

time command# 显示命令执行时间
seq 10			# 产生1到10的整数
md5sum			# 生成，或验证文件的MD5值


#### notice how replacement is happening
FILE=my_picture.jpg
echo $FILE
echo ${FILE/jpg/png}
echo ${FILE/my/your}
## my_picture.jpg
## my_picture.png
## your_picture.jpg


## array in bash
ArrayVariable=(words or things 'or stuff' separated 'by whitespace')
##single-quoted groups of words are treated as single values that will be assigned to an array element as a group. 
##(Though the same is not true of double-quoted values…)
echo $ArrayVariable
## below are the two ways to get all the values at once as a single string
echo ${ArrayVariable[*]} 
echo ${ArrayVariable[@]}
## check length
echo ${#ArrayVariable[@]}
for i in {0..5}; do
  echo $i: ${ArrayVariable[$i]}
done

# $(command) ,it means take the output of the command and insert it into the command line.

##grouped command, subshell
(cat FileA; echo xxxxx; cat FileB; NewVar=15) > Both
#the shell knows nothing about it FileB

{ cat FileA; echo xxxxx; cat FileB; NewVar=15; } > Both
#the shell does know nothing about it FileB, notice the blank space after { and before }, those are required


# here the cycled variable is "i"
for i in some things separated by whitespace; do
  commands involving $i
done
for i in oranges bananas apples; do
  echo "I like $i"
done


# this
test $VAR = small_and_sweet

# is equivalent to:
[ $VAR = "small_and_sweet" ];

#A[A-D]
#img.{png,jpg,svg}
#underscore, _, or a dash, -

#alias ls='ls -GFh'
#PATH=$PATH:/a-new/program/directory

## echo -e, use as regex, -E, use as it is
echo -e Hi "\t\t\t" There "\n\n"
echo -E Bye "\t\t\t" For now "\n\n"
## -n    Do not print the trailing newline character
echo -n Good to see you "\n\n"


## notice eval is need as otherwise = in the qutation is not recognized
# create variable that is a line with assignments
ASSIGNMENTS="D=four; E=five; F=six"

# evaluate that line after variable substitution:
eval $ASSIGNMENTS

# see that the variables D, E, and F have values assigned:
echo $D, $E, $F

##
# same as above, but quit processing the file as soon as you hit
# an @SQ line with a sequence name starting with "NW_"

##

awk '
  /^@SQ/ && /SN:NC_/ {print}
  /^@SQ/ && /SN:NW_/ {exit}  
' data/DPCh_plate1_F12_S72.sam


## notice how gzcat is used
# print only lines 101, and 103 from the fastq file
gzcat data/DPCh_plate1_F12_S72.R1.fq.gz | awk 'NR==101 || NR==102 {print}'


## [ ]	Match any single character from a character class (see below)


# Match either "big mess" or "huge mess"
/(big|huge) mess/
