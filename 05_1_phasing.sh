#!/bin/bash

# threads set to 32

projdir=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/full_imputation

#Install software
softwaredir=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/software/

chr=9
#Create/set working directory within project (one for each chromosome)
wdir="$projdir"/chr"$chr"_imputations
mkdir $wdir
refdir="$projdir"/reference
#Download QCed genotype data - should already be lifted over to hg19 in the previousstep
filename=full_imputation.lifted.ch.fl #name of the lifted over file from 03_liftover script

# USE THIS FILE FOR X/23
#filename=full_imputation.lifted.x.ch.fl

#gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/"$filename".lifted  $wdir/.

##Set environment variables

#Project name (i.e. bed file to be imputed)
filename=full_imputation

#Location of shapeit
shapeit5_loc="$softwaredir"/phase_common_static
#Location of minimac
minimac_loc="$softwaredir"/minimac4-4.1.6-Linux-x86_64/bin/minimac4
#Genetic map file (shapeit5 and minimac 4 use different genetic map files..
genetic_map_loc_shapeit="$refdir"/chr"$chr".b38.gmap.gz
genetic_map_loc_impute="$refdir"/genetic_map_chr"$chr"_combined_b38.txt
#SHAPEIT reference panel
reference_panel="$refdir"/HRC.r1-1.EGA.GRCh37.chr"$chr".impute.hg38Snochr.bcf.bgz

# # USE THIS FILE FOR X/23
#reference_panel="$wdir"/HRC.r1-1.EGA.GRCh37.chr"$chr".impute.hg38Snochr.xto23.bcf.bgz

#Set imputation window size
window=500000

#Make arrays for the chunk data to be looped over
chrom_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $1}'))
start_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $2}'))
stop_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $3}'))
num_chunks=(${#chrom_arr[@]})

mkdir phasing

# Phase each chunk for each set of subjects
for j in {00..06}; do
        # Store the start time into a file (then you can calculate the job length for a single set)
        date > phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}".start
        # Iterate through each chunk
        for i in $(seq $num_chunks); do
                # Set the start and stop values for the chunk and a region value
                current_start=${start_arr[$i-1]}
                current_stop=${stop_arr[$i-1]}
                region="${chrom_arr[$i-1]}":$((10#$current_start*1000000-$window))-$((10#$current_stop*1000000+$window))
                # Output the PLINK file as a bcf file
                plink --bfile chunking/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}" --recode vcf-iid bgz --out phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}"
                # Add the allele frequency, allele count, and total allele number tags to the bcf file
                bcftools +fill-tags phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".vcf.gz -Oz -o phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".annotated.vcf.gz -- -t AF,AC,AN
                # Index the bcf file
                tabix -f phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".annotated.vcf.gz
                # Run the phasing job (note: defaults are used)
                "$shapeit5_loc" --input phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".annotated.vcf.gz --reference "$reference_panel" --thread 32 --region "$region" --output phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.bcf --map "$genetic_map_loc_shapeit" --log phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.log
        done
		date > phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}".stop
done