#!/bin/bash
projdir=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/full_imputation

#Install software
softwaredir=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/software/
chr=1 # 2..22, X
#scp -P 242 jacqui@tiffany.ucsd.edu:/scratch2/HRC/HRC.r1-1.EGA.GRCh38.chr"$chr".header.impute.msav .
#Create/set working directory within project (one for each chromosome)
wdir="$projdir"/chr"$chr"_imputations
# File names for female/male X imputations
#wdir="$projdir"/chr"$chr"_imputations_f
#wdir="$projdir"/chr"$chr"_imputations_m

refdir="$projdir"/reference
#Download QCed genotype data - should already be lifted over to hg19 in the previousstep
filename=full_imputation.lifted.ch.fl #name of the lifted over file from 03_liftover script
# X/23 Imputations
#filename=full_imputation.lifted.x.ch.fl #name of the lifted over file from 03_liftover script

#gsutil -u $GOOGLE_PROJECT cp $WORKSPACE_BUCKET/data/"$filename".lifted  $wdir/.

##Set environment variables

#Project name (i.e. bed file to be imputed)
filename=full_imputation

#Location of shapeit
#shapeit5_loc="$softwaredir"/phase_common_static
#Location of minimac
minimac_loc="$softwaredir"/minimac4-4.1.6-Linux-x86_64/bin/minimac4
#MINIMAC reference panel
impute_ref_panel="$refdir"/HRC.r1-1.EGA.GRCh38.chr"$chr".header.impute.msav

#Unimputed over genotype data path
genotype_data="$filename".lifted
# X/23 Imputations
# genotype_data="$filename".lifted.x
#Set imputation window size
window=500000

#Make arrays for the chunk data to be looped over
chrom_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $1}'))
start_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $2}'))
stop_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $3}'))
num_chunks=(${#chrom_arr[@]})


## IMPUTE THE DATA
##### Note: you should reconfigure this environment to 32 CPUs and 28.8 GB of RAM;
##### The 00..02 loop was run first and 03..06 an hour later to optimize the thread usage

# Impute each set of subjects at a time
for j in {00..02}; do
#for j in {03..06}; do
	echo "$j"
	# Store the start time into a file (then you can calculate the job length for a single set)
	date > imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}".start
	# Iterate through each chunk
	for i in $(seq $num_chunks); do
		echo "$i"
		# Set the start and stop values for the chunk
		current_start=${start_arr[$i-1]}
		current_stop=${stop_arr[$i-1]}
		# Run the imputation job (note: ChunkLengthMb overrides the default length of 20 Mb for imputation, since some of the Snellius chunks are longer)
		echo "$impute_ref_panel"
		echo "${chrom_arr[$i-1]}"
		echo "$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.bcf
		echo $((10#$current_start*1000000))
		echo $((10#$current_stop*1000000))
		echo imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".imputed
		
		# USE THIS FOR CHR 1-22
		"$minimac_loc" --refHaps "$impute_ref_panel" --rsid --format HDS --chr "${chrom_arr[$i-1]}" --cpus 32 --ChunkLengthMb 32 --haps phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.bcf  --start $((10#$current_start*1000000)) --end $((10#$current_stop*1000000)) --window 500000  --prefix imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".imputed > imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".dose.vcf.log
		
		
		# USE THIS FOR CHR 23/X
		# Rename 23 to X for phased files
		#bcftools annotate --rename-chrs 23tox.txt phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.bcf --threads 32 -o phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.x.bcf
		# Use phased files with X annotation
		#bcftools index phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.x.bcf
		#"$minimac_loc" --refHaps "$impute_ref_panel" --rsid --format HDS --chr X --cpus 32 --ChunkLengthMb 32 --haps phasing/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".phased.x.bcf  --start $((10#$current_start*1000000)) --end $((10#$current_stop*1000000)) --window 500000  --prefix imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".imputed > imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}".dose.vcf.log
		
	done
	date > imputation/"$filename"."$j".chr"${chrom_arr[$i-1]}".stop
done