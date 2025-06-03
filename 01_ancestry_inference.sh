### ANCESTRY INFERENCE

### Environment parameters
#Set up the environment with 8 CPUs, 7.2 GB RAM and standard persistent disk (300 GB)

######################## DO THIS IF STARTING WITH NEW ENVIRONMENT ############################

# Set environmental variables on start

#Set working directory 
 directory=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/full_imputation

#Imputation project name
 filename=full_imputation

# Set up folders and create directories

## Make project directory (should match filename) 
mkdir $directory

## Go to project directory 
cd $directory

## Make software folder (where any reference panels and programs should be stored)
mkdir ../software

#########################################

## Upload EIGENSOFT package included in the AllOfUs folder into the software folder
### Compiling this program as downloaded from the Alkes Price lab website is a challenge since the OpenBLAS library cannot be root installed.
### To ensure proper compilation, the program will only compile CONVERTF instead of the full EIGENSOFT installation since CONVERTF is all you need for the pipeline.

## Upload SNPweights with the modified infer_ancestry script.

## Upload the "global_ancestry_kh" folder and unzip the global_ancestry_kh folder
## This will have the modified SNPWeights program that uses vectorization instead of for loops.
 
## Upload the AllOfUs annotation file (this maps AllOfUs IDs to rsIDs; cannot download for GDA from Illumina website since the data has been lifted)
### Annotation was done using the UCSC file All_20180418_hg38.vcf.

## After uploading all the programs/files/folders listed below, you are ready to start.
##      global_ancestry folder 
##      SNPweights folder
##      AllOfUs annotation file
##      EIGENSOFT package

cd ../software 

tar -xvzf EIG-7.2.1_KH.tar.gz
cd EIG-7.2.1/src
rm convertf
make

tar -xvzf SNPweights2.1_KH.tar.gz
 
## Create the ancestry folder in the project directory
mkdir ancestry
cd ancestry

# Set environmental variables

## Location of global ancestry folder
 global_ancestry_folder=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/software/global_ancestry_kh

## SNP keep list for ancestry calls (this was made by joining the SNP list and the conversion list to get a list in the original AllOfUs IDs)
 snp_keep_list=$global_ancestry_folder/AllOfUs_SNPs.keep

## CONVERTF location
 eigensoft_loc=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/software/EIG-7.2.1/src/convertf

## SNPweights path
 snpweights_path=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/software/SNPweights2.1/inferancestry_kh.py

## Location of panel
 snpweights_snplist=hgdp_kgp_merged_gsaqced_v3_jointsample_v4_k6.snplist
 snpweightfile_path=hgdp_kgp_merged_gsaqced_v3_jointsample_v4_k6.snpweightrefpanel
 snpweight_clustercenters=hgdp_kgp_merged_gsaqced_v3_jointsample_v4_k6_forsnpweights.snpweightrefpanel_clustercenters.csv

## rs-ID update file
 illumina_snplist=AllOfUs_SNPs_switched.annotated.complete
 rsidfile="$global_ancestry_folder"/"$illumina_snplist"_nodot_first

# Download the arrays BIM, BED, FAM files
 gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v7/microarray/plink/arrays.* $directory/ancestry

# Split the .fam file into 50,000 lines to run through the ancestry pipeline
 split -l50000 -d <(cat arrays.fam | awk '{print $1, $2}') full_imputation_keep_

# Iterate through each chunk of 50,000 people (note: this will delete temporary files due to low storage)
 for i in {00..06}; do

# Subset to the chunk of 50,000 people and SNP list for ancestry calls (though this is done later, it is more efficient and better for storage)
 plink2 --bfile arrays --keep full_imputation_keep_"$i" --extract $snp_keep_list --make-bed --out ancestry_"$i"

# Update the phenotype information (convertf won't consider anyone with missing phenotype information)
 cat ancestry_"$i".fam | awk '{print $1, $2, $3, $4, $5, 1}' > ancestry_"$i".fam.temp
 cat ancestry_"$i".fam.temp > ancestry_"$i".fam

# Update the rs-IDs for the subsetted PLINK
 Rscript --vanilla $global_ancestry_folder/scripts/make_rsid_update_file_v1_may5_2016.R $global_ancestry_folder/$illumina_snplist

# Run ancestry script (estimated run time for 50,000 subjects ~ 1 hr 50 min)
 chmod u+rwx $global_ancestry_folder/scripts/call_ancestry_v3_2023.sh
 bash $global_ancestry_folder/scripts/call_ancestry_v3_2023.sh ancestry_"$i" $rsidfile $directory/ancestry $global_ancestry_folder/$snpweights_snplist $snpweights_path $global_ancestry_folder/$snpweightfile_path $eigensoft_loc

# Combine all files into the main folder
 cat temporary_files/ancestry_"$i"_anc_*.predpc_oneweek > ancestry_"$i".predpc_oneweek

# Remove the temporary files
 rm ancestry_"$i".bed
 rm ancestry_"$i".bim
 rm ancestry_"$i".fam
 rm ancestry_"$i".fam.temp
 rm ancestry_"$i".log
 rm -r temporary_files

 done
 
# Concatenate all 50,000 runs into one file
 cat ancestry_*.predpc_oneweek > "$filename".predpc_oneweek

# Back this file up to project bucket
 gsutil -u $GOOGLE_PROJECT cp "$filename".predpc_oneweek $WORKSPACE_BUCKET/data

# Run the Rscript plot code to generate the plots and files
Rscript --vanilla "$global_ancestry_folder"/scripts/ancestry_plots_v4_mar1_2017.R "$filename".predpc_oneweek "$global_ancestry_folder"/"$snpweight_clustercenters" GSA