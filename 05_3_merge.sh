#!/bin/bash
# For set of subjects, concatenate the chromosome chunks together and index the file
filename=full_imputation
chr=1 # 2..22, X

for j in {00..06}; do
        ls imputation | grep "$j".chr"$chr" | grep .dose.vcf.gz | awk '{print "imputation/"$0}'> imputation/chr"$chr"."$j".imputation.vcf.txt
        bcftools concat -f imputation/chr"$chr"."$j".imputation.vcf.txt --threads 32 -Oz -o imputation/"$filename"."$j".chr"$chr".vcf.gz
        tabix -f imputation/"$filename"."$j".chr"$chr".vcf.gz
done

# Merge all the separate subjects sets into a full file stored in the merge folder
mkdir merged

ls imputation | grep .chr"$chr".vcf.gz | grep -v .vcf.gz.tbi | awk '{print "imputation/" $0}' > merged/chr"$chr".imputation.vcf.txt
bcftools merge -l merged/chr"$chr".imputation.vcf.txt --threads 32 -Oz -o merged/"$filename".chr"$chr".vcf.gz