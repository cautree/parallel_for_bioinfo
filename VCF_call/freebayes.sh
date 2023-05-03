#!/bin/bash
set -e
while read line;
do
name=$(basename $line)
  freebayes -f ref/Cannabis_sativa.fa bam/$line > freebayes_out/${name}.vcf
done < bam/bamfile

