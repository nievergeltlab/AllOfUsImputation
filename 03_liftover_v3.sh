### RSID UPDATE AND LIFTOVER

### Environment parameters
#Set up the environment with 8 CPUs, 7.2 GB RAM and standard persistent disk (500 GB works for largest chromosomes)

# Set environmental variables on start
## Set working directory 
 directory=/home/jupyter/workspaces/genomicarchitectureofptsdtinnitusandtbi/full_imputation
 cd $directory

## Create starting data directory
 mkdir starting_data

## Imputation project name
 filename=full_imputation

## Download the liftover information
 scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/infinium-global-diversity-array-8-v1-0_D1_b153_rsids.txt.snplist.hg38.mapped .
 scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/AllOfUs.chr_update .
 scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/AllOfUs.pos_update .
 
##Get checkflip and checkpos scripts, PLINK
 scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/software/* .
 
 #Allele frequencies in reference EUR pop
 scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/HRC.r1-1.EGA.GRCh37.allchr.impute.plink.EUR.frq2 .

# Download excluded variants/samples list from pre-imputation QC
 gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/"$filename".preimputation_excluded_samples .
 gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/"$filename".preimputation_excluded_variants .

#Download your list of European ancestry subjects (will need for allele frequency comparisons)
 gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/europeans.subjects .
  #Set the path to this file. 
   eur_subjects=europeans.subjects


# Filter to non-excluded subjects/variants
 plink2 --bfile starting_data/genotype_array_data_here --remove "$filename".preimputation_excluded_samples --exclude "$filename".preimputation_excluded_variants --make-bed --out "$filename"_preimpqc

# Remove starting data folder (free up storage ; or start with more space)
 rm -r starting_data

# Update variant names to rs IDs
## Make working directory
 mkdir updated


## Update to rs IDs
 

 #Update markers to rs-ids
  plink2 --bfile "$filename" --update-name infinium-global-diversity-array-8-v1-0_D1_b153_rsids.txt.snplist.hg38.mapped --make-bed --out "$filename"_rsupd

 #Identify all SNPs that do not have rs-ids
  grep -v rs "$filename"_rsupd.bim | awk '{print $2}' > non_updating_snps.snplist
  head non_updating_snps.snplist #This file should only contain SNPs without rs-ids
 
 #Update positions to hg19l; exclude SNPs without rs-ids
  plink2 --bfile updated/"$filename"_rsupd --exclude non_updating_snps.snplist --update-chr AllOfUs.chr_update --update-map AllOfUs.pos_update  --make-bed --out "$filename".lifted
 
 
#Check position overlap with reference panel. --ploc is the exact path to plink - keep in quotations
 perl checkpos6 --dbcol 1,2,3,4,5 --dbsnp HRC.r1-1.EGA.GRCh37.allchr.impute.plink.EUR.frq2 --ploc "plink/path/here"  "$filename".lifted

 #This will output a revised genotype file where only markers with reference matching positions are retained. 
 #IN the folder there is a .report file that summarizes everything. Check this out before proceeding.
 

#Flip non aligning alleles. 
 #Use the revised genotype file as the input. append .bim at the end of this filename
 #Notice that this uses your European subjects to calculate allele frequencies
 perl checkflip4.manc "$filename".lifted.ch --dbcol 0,3,4,5  --fth 0.15 --sfh 0.2 --info HRC.r1-1.EGA.GRCh37.allchr.impute.plink.EUR.frq2  --ploc "plink/path/here"  --keep $eur_subjects

#Save the output file, it is ready for imputation






#Test run (validated on tiffany
 #perl /scratch2/HRC/software/checkpos6 --dbcol 1,2,3,4,5 --dbsnp /scratch2/HRC/HRC.r1-1.EGA.GRCh37.allchr.impute.plink.EUR.frq2.AOU --ploc "/home/genetics/plink" pts_adni_mix_am-qc1.bim
 awk '{print $1,$2}' pts_adni_mix_am-qc1.ch.fam > pts_adni_mix_am-qc1.ch.fam.subs
 
  #perl /scratch2/HRC/software/checkflip4.manc pts_adni_mix_am-qc1.ch.bim --dbcol 0,3,4,5  --fth 0.15 --sfh 0.2 --info /scratch2/HRC/HRC.r1-1.EGA.GRCh37.allchr.impute.plink.EUR.frq2.AOU  --ploc  "/home/genetics/plink" --keep /home/adam/genotest/pts_adni_mix_am-qc1.ch.fam.subs

