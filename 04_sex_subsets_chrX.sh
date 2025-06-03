# The X chromsome is imputed on females and males separately
# Here we split the subject chunks by sex for phasing and imputation


# get sex information created from ipynb
# fid,iid,sex
gsutil -m cp gs://fc-secure-e67cc73f-de75-4314-b5bb-69d18bf10508/data/all_subjects_sex.csv .
# 413458 total records

# get set of female
awk -F "," '$3==2 {print $2}' all_subjects_sex.csv > array_female.tsv
# get set of male
awk -F "," '$3==1 {print $2}' all_subjects_sex.csv > array_male.tsv

# create subset of male and subset of females
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_00 | less
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_00 > sample_subsets_f/subset_00_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_01 > sample_subsets_f/subset_01_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_02 > sample_subsets_f/subset_02_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_03 > sample_subsets_f/subset_03_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_04 > sample_subsets_f/subset_04_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_05 > sample_subsets_f/subset_05_f
awk 'NR==FNR{a[$1]; next} $2 in a' array_female.tsv sample_subsets/subset_06 > sample_subsets_f/subset_06_f
wc -l subset_00_f subset_01_f subset_02_f subset_03_f subset_04_f subset_05_f subset_06_f
# 29850 subset_00_f
# 29832 subset_01_f
# 29681 subset_02_f
# 29750 subset_03_f
# 29115 subset_04_f
# 29312 subset_05_f
#  6086 subset_06_f
#183626 total

awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_00 > sample_subsets_f/subset_00_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_01 > sample_subsets_m/subset_01_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_02 > sample_subsets_m/subset_02_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_03 > sample_subsets_m/subset_03_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_04 > sample_subsets_m/subset_04_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_05 > sample_subsets_m/subset_05_m
awk 'NR==FNR{a[$1]; next} $2 in a' array_male.tsv sample_subsets/subset_06 > sample_subsets_m/subset_06_m
wc -l subset_00_m subset_01_m subset_02_m subset_03_m subset_04_m subset_05_m subset_06_m
# 19087 subset_00_m
# 19044 subset_01_m
# 19302 subset_02_m
# 19322 subset_03_m
# 19826 subset_04_m
# 19634 subset_05_m
#  3819 subset_06_m
#120034 total

# touch 23tox.txt
# file contains: 23	X
# update full_imputation.lifted.x.ch.fl.bim to use X instead of 23
plink2 --bfile full_imputation.lifted.x.ch.fl --sort-vars --update-chr 23tox.txt --make-bed --out full_imputation.lifted.x.ch.fl.23tox

# update chr label for bcf.bgz
zcat HRC.r1-1.EGA.GRCh37.chrX.impute.hg38Snochr.bcf.bgz | sed 's/X:/23:/g' > HRC.r1-1.EGA.GRCh37.chrX.impute.hg38Snochr.xto23id.bcf
nohup bcftools annotate --rename-chrs xto23.txt HRC.r1-1.EGA.GRCh37.chrX.impute.hg38Snochr.bcf.bgz -Oz -o HRC.r1-1.EGA.GRCh37.chr23.impute.hg38Snochr.xto23.bcf.bgz &



