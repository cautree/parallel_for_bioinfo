
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





