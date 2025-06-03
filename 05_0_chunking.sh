#!/bin/bash
#Make arrays for the chunk data to be looped over
chr=18 #23 for X
projdir=/home/jupyter/workspaces/duplicateofgenomicarchitectureofptsdtinnitusandtbi/full_imputation
refdir="$projdir"/reference

awk '$1=='$chr'{print $0}' "$refdir"/infosum_pos.chunks_10 > "$refdir"/infosum_pos.chunks_10.chr"$chr"


chrom_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $1}'))
start_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $2}'))
stop_arr=($(cat "$refdir"/infosum_pos.chunks_10.chr"$chr" | awk '{print $3}'))
num_chunks=(${#chrom_arr[@]})
filename=full_imputation
window=500000
wdir="$projdir"/chr"$chr"_imputations

cd $wdir
mkdir chunking

echo $num_chunks
for i in $(seq $num_chunks); do
    # Iterate through each set of subjects. The number of sets is 06 because that is the number of subject pieces, assuming thatN=50000 get imputed at a time
    for j in {00..06}; do
    # Set the chromosome start and stop positions
    current_start=${start_arr[$i-1]}
    current_stop=${stop_arr[$i-1]}
    # Subset the PLINK file to the set and chunk (note that you will need to multiply the start and stop positions to convert from bases to megabases

        plink2 --bfile "$projdir"/lifted/"$filename".lifted.ch.fl --chr ${chrom_arr[$i-1]} --from-bp $((10#$current_start*1000000-$window)) --to-bp $((10#$current_stop*1000000+$window)) --make-bed --keep "$projdir"/sample_subsets/subset_"$j" --out chunking/"$filename"."$j".chr"${chrom_arr[$i-1]}"."${start_arr[$i-1]}"_"${stop_arr[$i-1]}"

		# for X/23
		#plink2 --bfile "$projdir"/lifted/full_imputation.lifted.x.ch.fl --chr 23 --from-bp $((10#$current_start*1000000-$window)) --to-bp $((10#$current_stop*1000000+$window)) --make-bed --keep "$projdir"/sample_subsets_f/subset_"$j"_f --out chunking/"$filename"."$j".chr23."${start_arr[$i-1]}"_"${stop_arr[$i-1]}"
    done
 done