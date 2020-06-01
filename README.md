# filtering_GBA
Special quality control is necessary for GBA because of sequencing issues caused by GBAP1. This pipeline should be run on both the forced alignment GBA VCF (excluding the pseudogene) and the regular VCF, to compare.

````
#!/bin/bash

## RUN this in your local directory with your VCFs. Scripts are hard called. 

## Run after the VCF has been prepared by bioinformaticians (both forced alignment VCF and regular GBA):

###### STOP & READ: NEW_NAME variable must start with a capital letter. Part of a Dan Spiegelman script I can't change. 

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
