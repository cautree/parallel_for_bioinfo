bcftools view -Oz -s 221219-BioEcho-Solanum-lycopersicum-C03-7600000,\
221219-BioEcho-Solanum-lycopersicum-D03-7600000,\
221219-BioEcho-Solanum-lycopersicum-E03-7600000,\
221219-BioEcho-Solanum-lycopersicum-F03-7600000,\
221219-BioEcho-Solanum-lycopersicum-G03-7600000,\
221219-BioEcho-Solanum-lycopersicum-H03-7600000,\
221219-BioEcho-Solanum-lycopersicum-A04-7600000,\
221219-BioEcho-Solanum-lycopersicum-B04-7600000,\
221219-BioEcho-Solanum-lycopersicum-C04-7600000 tomato_impute-vcf-merged.vcf.gz --force-samples > filtered_tomato.vcf.gz



bcftools query -l filtered_tomato.vcf.gz
tabix -fp vcf filtered_tomato.vcf.gz



VCF=filtered_tomato.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--indep-pairwise 50 10 0.1 --out tomato


plink --vcf $VCF --double-id --allow-extra-chr \
--extract tomato.prune.in \
--make-bed --pca --out tomato





##this is for get the relatedness using plink genome
VCF=mpileup.vcf.gz
/Users/yanyan/Documents/software/plink_mac_20221210/plink --vcf $VCF --double-id --allow-extra-chr --indep-pairwise 50 10 0.1 --out cannabis_mpileup
/Users/yanyan/Documents/software/plink_mac_20221210/plink --vcf $VCF --double-id --allow-extra-chr \
--extract cannabis_mpileup.prune.in \
--recode transpose \
--out cannabis_mpileup 

/Users/yanyan/Documents/software/plink_mac_20221210/plink --tfile cannabis_mpileup --genome --out cannabis_mpileup --allow-extra-chr


