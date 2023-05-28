# who is in this file
bcftools query -l all.vcf.gz

#How many variants are in this file?
bcftools stats all.vcf.gz | less

#where are these variants
bcftools query -f '%CHROM\t%POS\n' all.vcf.gz | head 

#how many variants were on each of the chromosomes/scaffolds, sorted by number of variants
bcftools query -f '%CHROM\t%POS\n' all.vcf.gz | awk '{n[$1]++} END {for(i in n) print i, n[i]}' | sort -nbr -k 2 | less

# show the whole file from the top
bcftools view all.bcf | less

# of course, this works with either bcf or vcf or vcf.gz
bcftools view all.vcf.gz | less

# show just the header with -h.  Here we look at just the last 10 lines of the header
bcftools view -h all.bcf  | tail

# show the variants themselves (no header) with -H
bcftools view -H all.vcf.gz | head

# -O z: bgzipped VCF (vcf.gz)
# -O v: uncompressed VCF (the default)
# -O u: uncompressed BCF
# -O b: compressed BCF

## rename the samples in the file
bcftools reheader -s sample-renames.txt all.bcf  > renamed-all.bcf
bcftools reheader -s sample-renames.txt all.vcf.gz  > renamed-all.vcf.gz

##Extract keyed values from the INFO field
# extract CHROM POS and BaseQRankSum, separated by TABs
bcftools query -f '%CHROM\t%POS\t%INFO/BaseQRankSum\n' all.vcf.gz | less

# extract CHROM POS and total read depth DP
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' all.bcf | less

##extract information from each of the genotype columns
bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' all.bcf | less

##View data from specified regions, You can also specify those regions in a file with the -R option.
bcftools view -H -r CM031199.1:1-10000,CM031200.1:1000000-1005000 all.vcf.gz | less

## select samples
bcftools view -H -s s001,s002,s003 all.vcf.gz | less

## deselect samples
bcftools view -H -s ^s001,s002,s003 all.vcf.gz | less

##Combine VCF files in various ways
#Catenate VCF files
# make two files from different regions
bcftools view -O z -r CM031199.1:1-10000 all.vcf.gz  > A.vcf.gz
bcftools view -O z -r CM031200.1:1000000-1005000 all.vcf.gz  > B.vcf.gz

# how many variants in each of those?
bcftools stats A.vcf.gz | awk '/^SN/'
bcftools stats B.vcf.gz | awk '/^SN/'

# catenate the back together
bcftools concat -Oz  A.vcf.gz B.vcf.gz > CAT.vcf.gz

# how many variants in that?
bcftools stats CAT.vcf.gz | awk '/^SN/'

##Merge VCF files
##If you have files with different samples in them you can easily combine them:
# make file with first three samples
bcftools view -Oz -s s001,s002,s003 all.vcf.gz > first3.vcf.gz

# make another with the last three samples
bcftools view -Oz -s s006,s007,s008 all.bcf > last3.vcf.gz

# merging requires that the files be indexed
bcftools index first3.vcf.gz
bcftools index last3.vcf.gz

# merge those into a file with 6 samples
bcftools merge -Oz first3.vcf.gz last3.vcf.gz > 6-samples.vcf.gz


#Filter out variants for a variety of reasons
#Just the biallelic SNPs please Get things with no more than 2 alleles and no fewer than two alleles, and of a type = SNP:
bcftools view -Ou -m 2 -M 2 --types=snps all.bcf | bcftools stats - | awk '/^SN/'

#Just the biallelic indels please
# do it and see the result all in one line:
bcftools view -Ou -m 2 -M 2 --types=indels all.vcf.gz | bcftools stats - | awk '/^SN/'

#Fraction of missing sites less than X
#If you want to make sure that 60% of your individuals have at least one read at the genotype, you can do this:
bcftools view -i 'F_MISSING < 0.4' all.vcf.gz | bcftools stats - | awk '/^SN/'

#Exclude based on various features of the data
#You can use the -e option to bcftools view or bcftools filter to exclude sites that meet certain criteria.
#(You can use -i to include those sites and no others).
#only keep things with a maximum-likelihood-estimated allele frequency between 0.4 and 0.6
bcftools view -i 'INFO/MLEAF >= 0.4 && INFO/MLEAF <= 0.6' all.bcf | bcftools query -f '%INFO/MLEAF\n' | less

#retain only those sites at which the MLEAF value for the first alternate allele is between 0.4 and 0.6
bcftools view -i 'INFO/MLEAF[0] >= 0.4 && INFO/MLEAF[0] <= 0.6' all.bcf | bcftools query -f '%INFO/MLEAF\n' | less

#excluding those sites in which any individual had a DP less than 5
bcftools view -H -e 'FMT/DP < 5' all.bcf | less

#make it easier to see what the DPs are there, let’s print them:
bcftools view -e 'FMT/DP < 5' all.bcf | bcftools query -f '%CHROM\t%POS\t[%DP\t]\n' | less

#filter our data set down to variant sites at which at least two of the eight individuals had at least 10 reads of the first alternate allele
bcftools view -i 'COUNT(FMT/AD[:1] > 10) > 2' all.bcf | bcftools query -f '%CHROM\t%POS\t[%AD\t]\n' | less
