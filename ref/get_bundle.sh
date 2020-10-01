#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12

#############################
# get Broad hg38 BUNDLE data
#############################

LIST=$(cat <<'END_HEREDOC'
1000G_omni2.5.hg38.vcf.gz
1000G_omni2.5.hg38.vcf.gz.tbi
1000G_phase1.snps.high_confidence.hg38.vcf.gz
1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
Homo_sapiens_assembly38.dict
Homo_sapiens_assembly38.fasta.64.alt
Homo_sapiens_assembly38.fasta.fai
Homo_sapiens_assembly38.fasta.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
dbsnp_138.hg38.vcf.gz
dbsnp_138.hg38.vcf.gz.tbi
dbsnp_144.hg38.vcf.gz
dbsnp_144.hg38.vcf.gz.tbi
dbsnp_146.hg38.vcf.gz
dbsnp_146.hg38.vcf.gz.tbi
hapmap_3.3.hg38.vcf.gz
hapmap_3.3.hg38.vcf.gz.tbi
hapmap_3.3_grch38_pop_stratified_af.vcf.gz
hapmap_3.3_grch38_pop_stratified_af.vcf.gz.tbi
wgs_calling_regions.hg38.interval_list
beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
beta/Homo_sapiens_assembly38.known_indels.vcf.gz
beta/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
END_HEREDOC
)

# create folder and get data
#mkdir -p reference

# Broad bundle repo for hg38
bundleurl=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38

for file in $LIST; do
echo "# downloading $url"
wget -np --ftp-user=gsapubftp-anonymous ${bundleurl}/$file
#wget -P reference -np --ftp-user=gsapubftp-anonymous ${bundleurl}/$file
done

# parallel version for speed and retry on fail
# wgetopt="--retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 --continue"
# echo $LIST | xargs -I {} -n 1 -P ${nthr} sh -c "wget -P reference --ftp-user=gsapubftp-anonymous ${wgetopt} -np ${bundleurl}/{}"

# decompress reference fasta
gunzip -k reference/Homo_sapiens_assembly38.fasta.gz

# create PICARD dict from assembly
java -jar $PICARD/picard.jar CreateSequenceDictionary \
	R=reference/Homo_sapiens_assembly38.fasta \
	O=reference/Homo_sapiens_assembly38.dict2

# create chr22 BED interval file
gawk 'BEGIN{FS="\t"; OFS="\t"}{if ($2 ~/chr22/) {split($2,chr,":"); split($3,len,":"); print chr[2],0,len[2]}}' \
	reference/Homo_sapiens_assembly38.dict \
	> reference/chr22.bed

# add extra chr22 files from our GIT repo
wget -P reference -np https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020/raw/master/data/addedrefs.tgz &&\
tar -xzvf reference/addedrefs.tgz

# touch all tbi files to prevent date-tag issues
touch reference/*.tbi
