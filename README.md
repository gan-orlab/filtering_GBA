# filtering_GBA
Special quality control is necessary for GBA because of sequencing issues caused by GBAP1. This pipeline should be run on both the forced alignment GBA VCF (excluding the pseudogene) and the regular VCF, to compare.

## Start by annotating your VCF and selecting possible pathogenic variants.

````
#!/bin/bash

## RUN this in your local directory with your VCFs. Scripts are hard called. 

## Run after the VCF has been prepared by bioinformaticians (both forced alignment VCF and regular GBA):

###### STOP & READ: NEW_NAME variable must start with a capital letter. 
## It's part of a Dan Spiegelman script I can't change. 

# If extracting from normal VCF: --chr 1 --from-bp 155204139 --to-bp 155214753

VCF=$1 # name of the noGBAP1 VCF from Dan/Eric. 
NEW_NAME=$2 # what you want your VCF to be called after annotation. 
# it will be affixed with ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf

N=$3 # how many samples in this batch

## ANNOTATE

bash ~/runs/krohn/krohn/Home/annovar/annovar_MINIMAL.sh $VCF $NEW_NAME 
# this script annotates with GENE, gnomad_genome, and clinvar

# extract likely pathogenic positions:
grep nonsynonymous ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > nonsynonymous_variants.ls
grep frameshift ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > frameshift_variants.ls
grep stopgain ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > stopgain_variants.ls
grep splicing ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > splicing_variants.ls

cat frameshift_variants.ls nonsynonymous_variants.ls stopgain_variants.ls splicing_variants.ls > extract_possibly-pathogenic_variants.ls
````

## Next, parse the VCF to extract relevant information.
The most important piece of information here is the depth of coverage per allele (AD). 
Unfortunately, some sites don't have this information, and only have total depth of coverage (DP). 
For those sites, we extract the DP and the genotype call (GT). We don't have a way to filter these by
percentage of alternate allele call, so these sites should be treated skeptically. Ideally we can validate them with
another VCF from MIPs or with omnix GWAS data. 

````
# extract the coverage from the VCF. 
# it's essentially just parsing the information in the VCF. 
# AD = depth of coverage per allele
# DP = total depth of coverage
# GT = genotype call 
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="AD") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > AD_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="DP") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > DP_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="GT") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > GT_all_noGBAP1.vcf.tab
````

