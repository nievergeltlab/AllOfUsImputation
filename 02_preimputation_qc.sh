### PREIMPUTATION QC - this was done on all subjects regardless of WGS availability

### Environment parameters
#Set up the environment with 8 CPUs, 7.2 GB RAM and standard persistent disk (500 GB works for largest chromosomes)

# Set environmental variables on start

# Set working directory 
 directory=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/full_imputation
 cd $directory

# Create starting data directory
 mkdir starting_data

#Imputation project name
 filename=full_imputation

# Download the arrays PLINK files if working from new environment
 gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/microarray/plink/arrays.* $directory/starting_data

# Update sex in pheno file with self-reported sex at birth 
## Generate the file "arrays_updated.fam" by running the Jupyter notebook "fam_file_construction" (note this was originally run under v6 phenotypic data)

# ## Move the file into starting_data folder
 # mv ../arrays_updated.fam starting_data/

# ## Create the update sex file and the pheno file
 # cat starting_data/arrays_updated.fam | awk '{print $1, $2, $5}' > starting_data/"$filename".update_sex

# ## Create pheno file
 # cat starting_data/arrays_updated.fam | awk '{print $1, $2, $6}' > starting_data/"$filename".pheno

## Update the phenotypic information
 plink2 --bfile starting_data/arrays --update-name starting_data/AllOfUs_SNPs_switched.annotated.complete --update-sex starting_data/"$filename".update_sex --pheno starting_data/"$filename".pheno --make-just-fam --out starting_data/arrays_upd."$filename"
 cat starting_data/arrays_upd."$filename".fam > starting_data/arrays.fam

# Create folder for pre-imputation QC
 mkdir preimputation_qc

# Create the initial list of SNPs to get the original count and write to the QC log
 plink2 --bfile starting_data/arrays --write-snplist --out preimputation_qc/"$filename"_qc0

 echo "Number of variants at start:" > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc0.snplist | wc -l) > preimputation_qc/"$filename"_qc.log 

# Using the fam file, get the number of subjects, cases, and controls
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of subjects at start:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat starting_data/arrays.fam | wc -l) > preimputation_qc/"$filename"_qc.log

 cat preimputation_qc/"$filename"_qc.log <(echo "Number of cases at start:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat starting_data/arrays.fam | awk '{if ($6==2) print $0}' | wc -l) > preimputation_qc/"$filename"_qc.log

 cat preimputation_qc/"$filename"_qc.log <(echo "Number of controls at start:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat starting_data/arrays.fam | awk '{if ($6==1) print $0}' | wc -l) > preimputation_qc/"$filename"_qc.log

# STEP 1: Filter to SNPs with a call rate above 0.95 for an initial cutoff
## Get missing rates from PLINK
 plink2 --bfile starting_data/arrays --missing variant-only --threads 8 --out preimputation_qc/"$filename"_qc1

## Create included and excluded SNP lists (.info has the specific call rate for reference) 
 cat preimputation_qc/"$filename"_qc1.vmiss | awk '{if ($5 > 0.05) print $2}' | tail -n +2 > preimputation_qc/"$filename"_qc1.excluded_snplist
 cat preimputation_qc/"$filename"_qc1.vmiss | awk '{if ($5 <= 0.05) print $2}' | tail -n +2 > preimputation_qc/"$filename"_qc1.included_snplist
 cat preimputation_qc/"$filename"_qc1.vmiss | awk '{if ($5 > 0.05) print $0}' > preimputation_qc/"$filename"_qc1.excluded_info

## Update the log wtih the number of SNPs after the step 1 filter
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of variants after 0.95 SNP call rate filter:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc1.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

# STEP 2: Filter to subjects with a call rate above 0.98
## Get missing rates for subjects from PLINK
 plink2 --bfile starting_data/arrays --exclude preimputation_qc/"$filename"_qc1.excluded_snplist --missing sample-only --threads 8 --out preimputation_qc/"$filename"_qc2

## Create included and excluded subject lists (.info has call rate for reference)
 cat preimputation_qc/"$filename"_qc2.smiss | awk '{if ($6 > 0.02) print $1,$2}' | tail -n +2 > preimputation_qc/"$filename"_qc2.excluded_samples
 cat preimputation_qc/"$filename"_qc2.smiss | awk '{if ($6 <= 0.02) print $1,$2}' | tail -n +2 > preimputation_qc/"$filename"_qc2.included_samples
 cat preimputation_qc/"$filename"_qc2.smiss | awk '{if ($6 > 0.02) print $0}' > preimputation_qc/"$filename"_qc2.excluded_info

## Update the log with the number of subjects removed after step 2 filter
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of subjects after 0.98 call rate filter:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc2.included_samples | wc -l) > preimputation_qc/"$filename"_qc.log

# STEP 3: Filter out subjects with FHET outside +/- 0.20
## Get heterozygosity of subjects
 plink --bfile starting_data/arrays --het --remove preimputation_qc/"$filename"_qc2.excluded_samples --exclude preimputation_qc/"$filename"_qc1.excluded_snplist --out preimputation_qc/"$filename"_qc3

## Create included and excluded subject lists (.info has heterozygosity values for reference)
 cat preimputation_qc/"$filename"_qc3.het | tail -n +2 | awk '{if ($6>0.2) print $1,$2; else if ($6<-0.2) print $1,$2}' > preimputation_qc/"$filename"_qc3.excluded_samples
 cat preimputation_qc/"$filename"_qc3.het | tail -n +2 | awk '{if ($6<=0.2 && $6>=-0.2) print $1,$2}' > preimputation_qc/"$filename"_qc3.included_samples
 cat preimputation_qc/"$filename"_qc3.het | awk '{if ($6>0.2) print $0; else if ($6<-0.2) print $0}' > preimputation_qc/"$filename"_qc3.excluded_info

## Update the log with the number of subjects removed after step 3 filter
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of subjects after FHET filter:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc3.included_samples | wc -l) > preimputation_qc/"$filename"_qc.log

## Concatenate the step 2 and step 3 excluded subjects lists to get a growing list of excluded subjects
 cat preimputation_qc/"$filename"_qc2.excluded_samples preimputation_qc/"$filename"_qc3.excluded_samples > preimputation_qc/"$filename"_qc3.excluded_samples_combined

# STEP 4: Exclude subjects with a mismatch between self-reported sex at birth and genotypic sex
## Run sex check on PLINK
 plink --bfile starting_data/arrays --remove preimputation_qc/"$filename"_qc3.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc1.excluded_snplist --check-sex --out preimputation_qc/"$filename"_qc4

## Create list of included and excluded subjects (.info has values for reference)
 cat preimputation_qc/"$filename"_qc4.sexcheck | awk '{if($3==1 && $4==2) print $1,$2; else if ($3==2 && $4==1) print $1,$2}' > preimputation_qc/"$filename"_qc4.excluded_samples
 cat preimputation_qc/"$filename"_qc4.sexcheck | awk '{if($5=="OK") print $1,$2; else if ($3==0) print $1,$2; else if ($4==0) print $1,$2}' > preimputation_qc/"$filename"_qc4.included_samples
 cat preimputation_qc/"$filename"_qc4.sexcheck | awk '{if($3==1 && $4==2) print $0; else if ($3==2 && $4==1) print $0}' > preimputation_qc/"$filename"_qc4.sexcheck.excluded_info

## Update the log with the number of subjects excluded for sex mismatch
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of subjects after sex mismatch exclusion:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc4.included_samples | wc -l) > preimputation_qc/"$filename"_qc.log

## Concatenate the step 4 excluded subjects to the growing list of subjects
 cat preimputation_qc/"$filename"_qc3.excluded_samples_combined preimputation_qc/"$filename"_qc4.excluded_samples > preimputation_qc/"$filename"_qc4.excluded_samples_combined

# STEP 5: Create list of subjects with sex warnings (undefined phenotypic or genotypic sex)
## Create list of subjects with sex warning, including type (phenotypic, genotypic, or both)
 cat preimputation_qc/"$filename"_qc4.sexcheck | awk '{if ($3==0 && $4==0) print $0,"phenotype_genotype_ambiguity"; else if ($3==0) print $0,"phenotype_ambiguity"; else if ($4==0) print $0,"genotype_ambiguity"}' > preimputation_qc/"$filename"_qc5.sexcheck_warnings

## Add the number of subjects with sex ambiguity
cat preimputation_qc/"$filename"_qc.log <(echo "Number of subjects with sex warning:") > preimputation_qc/"$filename"_qc_temp.log
cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc5.sexcheck_warnings | wc -l) > preimputation_qc/"$filename"_qc.log

# STEP 6: Exclude SNPs with call rate below 0.98
## Get list of call rates for variant-only missing values
 plink2 --bfile starting_data/arrays --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc1.excluded_snplist --missing variant-only --out preimputation_qc/"$filename"_qc6

## Create lists of included and excluded variants
 cat preimputation_qc/"$filename"_qc6.vmiss | awk '{if ($5 > 0.02) print $2}' | tail -n +2 > preimputation_qc/"$filename"_qc6.excluded_snplist
 cat preimputation_qc/"$filename"_qc6.vmiss | awk '{if ($5 <= 0.02) print $2}' | tail -n +2 > preimputation_qc/"$filename"_qc6.included_snplist
 cat preimputation_qc/"$filename"_qc6.vmiss | awk '{if ($5 > 0.02) print $0}' > preimputation_qc/"$filename"_qc6.excluded_info

## Add number of variants after filter to log
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of variants after 0.98 SNP call rate filter:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc6.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

## Combine lists of excluded variants from step 1 and step 6 to create gorwing list of excluded variants
 cat preimputation_qc/"$filename"_qc1.excluded_snplist preimputation_qc/"$filename"_qc6.excluded_snplist > preimputation_qc/"$filename"_qc6.excluded_snplist_combined

# STEP 7: Exclude subjects with a missing rate difference between cases and controls >= 0.02
## Get missing rates for cases
 cat starting_data/arrays.fam | awk '{if ($6 == 2) print $1, $2}' > preimputation_qc/"$filename"_qc7_cases
 plink2 --bfile starting_data/arrays --keep preimputation_qc/"$filename"_qc7_cases --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc6.excluded_snplist_combined --missing variant-only --out preimputation_qc/"$filename"_qc7_cases

## Get missing rates for controls
 cat starting_data/arrays.fam | awk '{if ($6 == 1) print $1, $2}' > preimputation_qc/"$filename"_qc7_cntls
 plink2 --bfile starting_data/arrays --keep preimputation_qc/"$filename"_qc7_cntls --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc6.excluded_snplist_combined --missing variant-only --out preimputation_qc/"$filename"_qc7_cntls

## Sort the lists
 cat preimputation_qc/"$filename"_qc7_cases.vmiss | awk '{print $2,$5}' | tail -n +2 | LC_ALL=C sort -g -k 1,1b > preimputation_qc/"$filename"_qc7_cases.vmiss_sorted
 cat preimputation_qc/"$filename"_qc7_cntls.vmiss | awk '{print $2,$5}' | tail -n +2 | LC_ALL=C sort -g -k 1,1b > preimputation_qc/"$filename"_qc7_cntls.vmiss_sorted

## Join the lists and create list of included and excluded variants 
 LC_ALL=C join preimputation_qc/"$filename"_qc7_cases.vmiss_sorted preimputation_qc/"$filename"_qc7_cntls.vmiss_sorted | awk '{print $1,$2,$3,$2-$3}' | awk '{if($4>0.02) print $1; else if ($4<-0.02) print $1}' > preimputation_qc/"$filename"_qc7.excluded_snplist
 LC_ALL=C join preimputation_qc/"$filename"_qc7_cases.vmiss_sorted preimputation_qc/"$filename"_qc7_cntls.vmiss_sorted | awk '{print $1,$2,$3,$2-$3}' | awk '{if($4<=0.02 && $4>=-0.02) print $1}' > preimputation_qc/"$filename"_qc7.included_snplist
 LC_ALL=C join preimputation_qc/"$filename"_qc7_cases.vmiss_sorted preimputation_qc/"$filename"_qc7_cntls.vmiss_sorted | awk '{print $1,$2,$3,$2-$3}' | awk '{if($4>0.02) print $0; else if ($4<-0.02) print $0}' | cat <(echo "ID CASE_MISS CNTL_MISS DIFF_CASE_CNTL") - > preimputation_qc/"$filename"_qc7.excluded_info

## Add number of variants after filter to log 
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of variants after 0.02 SNP call rate difference filter:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc7.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

## Add excluded variants to growing list of excluded variants
 cat preimputation_qc/"$filename"_qc6.excluded_snplist_combined preimputation_qc/"$filename"_qc7.excluded_snplist > preimputation_qc/"$filename"_qc7.excluded_snplist_combined

# STEP 8: Exclude invariant SNPs
## Get allele counts for SNPs
 plink --bfile starting_data/arrays --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc7.excluded_snplist_combined --freq counts --out preimputation_qc/"$filename"_qc8

## Create a list of included and excluded variants
 cat preimputation_qc/"$filename"_qc8.frq.counts | tail -n +2 | awk '{if ($5==0) print $2; else if ($6==0) print $2}' > preimputation_qc/"$filename"_qc8.excluded_snplist
 cat preimputation_qc/"$filename"_qc8.frq.counts | tail -n +2 | awk '{if ($5!=0 && $6!=0) print $2}' > preimputation_qc/"$filename"_qc8.included_snplist
 cat preimputation_qc/"$filename"_qc8.frq.counts | tail -n +2 | awk '{if ($5==0) print $0; else if ($6==0) print $0}' > preimputation_qc/"$filename"_qc8.excluded_info

## Add number of variants after filter to log
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of invariant SNPs:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc8.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

## Add excluded variants to growing list of excluded variants
 cat preimputation_qc/"$filename"_qc7.excluded_snplist_combined preimputation_qc/"$filename"_qc8.excluded_snplist > preimputation_qc/"$filename"_qc8.excluded_snplist_combined

# STEP 9: Exclude variants with HWE for controls <= -6
## Filter to EUR control subjects (largest ancestry group as v7 release) and calculate HWE
 gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/eur_subjects.txt .
 LC_ALL=C join <(LC_ALL=C sort eur_subjects | awk '{print $1":"$2}') <(LC_ALL=C sort preimputation_qc/"$filename"_qc7_cntls | awk '{print $1":"$2}') | awk 'BEGIN {FS=":"}{print $1, $2}' > preimputation_qc/"$filename"_qc9_cntls
 plink2 --bfile starting_data/arrays --keep preimputation_qc/"$filename"_qc9_cntls --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc8.excluded_snplist_combined --hardy --out preimputation_qc/"$filename"_qc9_cntls

## Create a list of excluded/included autosomal variants
 cat preimputation_qc/"$filename"_qc9_cntls.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -6) print $2}' > preimputation_qc/"$filename"_qc9_cntls.autosome.excluded_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 >= -6) print $2}' > preimputation_qc/"$filename"_qc9_cntls.autosome.included_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -6) print $0}' > preimputation_qc/"$filename"_qc9_cntls.autosome.info

## Create a list of excluded/included X chromosome variants
 cat preimputation_qc/"$filename"_qc9_cntls.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -6) print $2}' > preimputation_qc/"$filename"_qc9_cntls.x.excluded_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 >= -6) print $2}' > preimputation_qc/"$filename"_qc9_cntls.x.included_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -6) print $0}' > preimputation_qc/"$filename"_qc9_cntls.x.info

## Concatenate the separate lists to make included/excluded variants
 cat preimputation_qc/"$filename"_qc9_cntls.autosome.excluded_snplist preimputation_qc/"$filename"_qc9_cntls.x.excluded_snplist > preimputation_qc/"$filename"_qc9.excluded_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.autosome.included_snplist preimputation_qc/"$filename"_qc9_cntls.x.included_snplist > preimputation_qc/"$filename"_qc9.included_snplist
 cat preimputation_qc/"$filename"_qc9_cntls.autosome.info preimputation_qc/"$filename"_qc9_cntls.x.info > preimputation_qc/"$filename"_qc9.excluded_info

## Add number of variants after filter to log
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of variants with HWE for controls >= -6:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc9.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

## Add to growing list of excluded variants
 cat preimputation_qc/"$filename"_qc8.excluded_snplist_combined preimputation_qc/"$filename"_qc9.excluded_snplist > preimputation_qc/"$filename"_qc9.excluded_snplist_combined

# STEP 10: Exclude variants with HWE for cases <= -20
## Filter to EUR case subjects and calculate HWE <= -20
 LC_ALL=C join <(LC_ALL=C sort eur_subjects.txt | awk '{print $1":"$2}') <(LC_ALL=C sort preimputation_qc/"$filename"_qc7_cases | awk '{print $1":"$2}') | awk 'BEGIN {FS=":"}{print $1, $2}' > preimputation_qc/"$filename"_qc10_cases
 plink2 --bfile starting_data/arrays --keep preimputation_qc/"$filename"_qc10_cases --remove preimputation_qc/"$filename"_qc4.excluded_samples_combined --exclude preimputation_qc/"$filename"_qc9.excluded_snplist_combined --hardy --out preimputation_qc/"$filename"_qc10_cases

## Create a list of excluded/included autosomal variants
 cat preimputation_qc/"$filename"_qc10_cases.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -20) print $2}' > preimputation_qc/"$filename"_qc10_cases.autosome.excluded_snplist
 cat preimputation_qc/"$filename"_qc10_cases.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 >= -20) print $2}' > preimputation_qc/"$filename"_qc10_cases.autosome.included_snplist
 cat preimputation_qc/"$filename"_qc10_cases.hardy | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -20) print $0}' > preimputation_qc/"$filename"_qc10_cases.autosome.info

## Create a list of excluded/included X chromosome variants
 cat preimputation_qc/"$filename"_qc10_cases.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -20) print $2}' > preimputation_qc/"$filename"_qc10_cases.x.excluded_snplist
 cat preimputation_qc/"$filename"_qc10_cases.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 >= -20) print $2}' > preimputation_qc/"$filename"_qc10_cases.x.included_snplist
 cat preimputation_qc/"$filename"_qc10_cases.hardy.x | awk '{print $0,log($10)/log(10)}' | awk '{if ($11 < -20) print $0}' > preimputation_qc/"$filename"_qc10_cases.x.info 

## Concatenate the separate lists to make included/excluded variants
 cat preimputation_qc/"$filename"_qc10_cases.autosome.excluded_snplist preimputation_qc/"$filename"_qc10_cases.x.excluded_snplist > preimputation_qc/"$filename"_qc10.excluded_snplist
 cat preimputation_qc/"$filename"_qc10_cases.autosome.included_snplist preimputation_qc/"$filename"_qc10_cases.x.included_snplist > preimputation_qc/"$filename"_qc10.included_snplist
 cat preimputation_qc/"$filename"_qc10_cases.autosome.info preimputation_qc/"$filename"_qc10_cases.x.info > preimputation_qc/"$filename"_qc10.excluded_info

## Add number of variants after filter to log
 cat preimputation_qc/"$filename"_qc.log <(echo "Number of variants with HWE for controls >= -20:") > preimputation_qc/"$filename"_qc_temp.log
 cat preimputation_qc/"$filename"_qc_temp.log <(cat preimputation_qc/"$filename"_qc10.included_snplist | wc -l) > preimputation_qc/"$filename"_qc.log

## Add to growing list of excluded variants
 cat preimputation_qc/"$filename"_qc9.excluded_snplist_combined preimputation_qc/"$filename"_qc10.excluded_snplist > preimputation_qc/"$filename"_qc10.excluded_snplist_combined

# Save the combined variant and sample exclusion lists to the workspace bucket
 cat preimputation_qc/"$filename"_qc10.excluded_snplist_combined > "$filename".preimputation_excluded_variants
 gsutil -u $GOOGLE_PROJECT cp "$filename".preimputation_excluded_variants $WORKSPACE_BUCKET/data

 cat preimputation_qc/"$filename"_qc4.excluded_samples_combined > "$filename".preimputation_excluded_samples
 gsutil -u $GOOGLE_PROJECT cp "$filename".preimputation_excluded_samples $WORKSPACE_BUCKET/data

