#!/bin/bash
module load anaconda/2

#Input arp files

arg_arr=("$@") 
echo "${arg_arr[@]}"
fsc2vcf=fsc2plink_vcf_missing_genos_chr.py

# Convert with python to vcfs

python $fsc2vcf "${arg_arr[@]}"

gzip $1

gzip *vcf
#gzip *arp
mkdir pop2 pop1_ref
mkdir pop2/unphased pop2/phased


mv *Pop2* pop2/phased
mv *unphased* pop2/unphased
mv *Pop1* pop1_ref

# Create chr targ chunk files for impute.py
ls pop2/unphased |egrep -o "$2"_[[:digit:]]*maf_[[:digit:]]*missing_ind_2_[[:digit:]]*_2_[[:digit:]]+_chr|sort -u > pop2.inds
sed -i s/$/\!_unphased.vcf.gz/ pop2.inds
split -l 10 pop2.inds chr_targs_pop2_

ls pop1_ref |egrep -o test_[[:digit:]]*maf_ind_1_[[:digit:]]*_1_[[:digit:]]+_chr|sort -u > pop1.inds
sed -i s/$/\!_unphased.vcf.gz/ pop1.inds
split -l 10 pop1.inds chr_targs_pop1_

mkdir output

for j in $(ls|grep chr_targs_pop2); do sh loop.sh 1 C $j; done



for j in {1..49}; do
	ls -1v|grep "$2"*maf_[[:digit:]]*missing_ind_2_"$j"_|sort -V > ind_tmp.list

	for i in $(cat ind_tmp.list); do
		zcat $i|grep -v '^#' >> "$2"_ind_2_"$j"_2_"$j"_allchrs_imputed.vcf
		gzip "$2"_ind_2_"$j"_2_"$j"_allchrs_imputed.vcf
	done

done


#for i in $(ls pop1_ref); do vcftools --gzvcf pop1_ref/$i --freq --out pop1_ref/$i ;done
