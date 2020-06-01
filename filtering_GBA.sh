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

# extract the coverage from the VCF. 
# it's essentially just parsing the information in the VCF. 
# AD = depth of coverage per allele
# DP = total depth of coverage
# GT = genotype call 
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="AD") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > AD_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="DP") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > DP_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="GT") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > GT_all_noGBAP1.vcf.tab

# extract your list of pathogenic variants from the coverage file:
grep -f extract_possibly-pathogenic_variants.ls AD_all_noGBAP1.vcf.tab > to-tr_AD_patho_noGBAP1.vcf.tab

# see which variants don't have AD. Need another approach. 
cut -f2 to-tr_AD_patho_noGBAP1.vcf.tab | sort | uniq > vars_that_made_it.ls
sort extract_possibly-pathogenic_variants.ls | uniq > all_patho_vars.ls
comm -13 vars_that_made_it.ls all_patho_vars.ls > to-DP.ls

grep -f to-DP.ls DP_all_noGBAP1.vcf.tab > DP_select_noGBAP1.vcf.tab
grep -f to-DP.ls GT_all_noGBAP1.vcf.tab | cut -f4 > paste_this
paste DP_select_noGBAP1.vcf.tab paste_this > DP_GT_noGBAP1.vcf.tab

rm vars_that_made_it.ls all_patho_vars.ls paste_this DP_select_noGBAP1.vcf.tab to-DP.ls

# extract only variant sites (remove wild type or no call)
awk '{if ($5 != "0/0" && $4 != "./.") print $1="SPLIT",$2,$3,$4,$5}' DP_GT_noGBAP1.vcf.tab > carriers_DP_select_noGBAP1.vcf.tab

# separate the AD calls into two columns:
tr ',' '\t' < to-tr_AD_patho_noGBAP1.vcf.tab > temp_AD_patho_noGBAP1.vcf.tab

# extract tri-allelic positions that seem real
awk '{if ($6>0 && $5==0) print $1="TRI",$2,$3,$4,$6; else if ($6>0 && $5>0) print $1="TRI",$2,$3,$4,$5=0; else if ($6<0) print $1,$2,$3,$4,$5}' temp_AD_patho_noGBAP1.vcf.tab > AD_patho_noGBAP1.vcf.tab 

rm to-tr_AD_patho_noGBAP1.vcf.tab
rm temp*


# if A2 has 0 coverage, there are no carriers of this variant (all coverage on wild type):
awk '{if ($5 > 0 && $4>=0) print $0}' AD_patho_noGBAP1.vcf.tab > carriers_AD_patho_noGBAP1.vcf.tab

# make a list of all positions you're investigating: 
cut -f2 -d ' ' carriers_AD_patho_noGBAP1.vcf.tab | sort | uniq > carrier_variants_to-investigate_AD.ls
cut -f2 -d ' ' carriers_DP_select_noGBAP1.vcf.tab | sort | uniq > carrier_variants_to-investigate_DP.ls

## ADD VARIANT NAMES

# create input_annovar table (specifications for this on annovar website):
grep -f carrier_variants_to-investigate_AD.ls ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f1,2,4,5 > almost_input_annovar.txt

# for AD: if triallelic, take the second allele
awk '{print $1,$2,$2,$3,$4}' almost_input_annovar.txt > to-tr_almost_input_annovar.txt
tr ',' ' ' < to-tr_almost_input_annovar.txt > almost_input_annovar.txt
awk '{if ($6>0) print $1,$2,$3,$4,$6; else print $1,$2,$3,$4,$5}' almost_input_annovar.txt > input_annovar1.txt

grep -f carrier_variants_to-investigate_DP.ls ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f1,2,4,5 > almost_input_annovar.txt

# for DP: if triallelic, we're taking the first allele. 
# THIS IS AN ASSUMPTION. Check manually that the end (indicated by "SPLIT" in chr)
awk '{print $1,$2,$2,$3,$4}' almost_input_annovar.txt > to-tr_almost_input_annovar.txt
tr ',' ' ' < to-tr_almost_input_annovar.txt > almost_input_annovar.txt
awk '{print $1,$2,$3,$4,$5}' almost_input_annovar.txt > input_annovar2.txt

cat input_annovar1.txt input_annovar2.txt > input_annovar.txt


# run table_annovar.pl: 
bash ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/table_annovar.sh annotations_patho

# parse to get relevant information in clean table:
cut -f10 annotations_patho.hg19_multianno.txt > parse_proChange.tab
sed -i 's/':'/\t/g' parse_proChange.tab
cut -f11,13 parse_proChange.tab > parse_again.tab
sed -i 's/','/\t/g' parse_again.tab
awk '{print $1,$2}' parse_again.tab > paste_proChange.tab
awk '{print $1,$2,$9}' annotations_patho.hg19_multianno.txt > to_paste_anno.txt
paste to_paste_anno.txt paste_proChange.tab > likely_pathogenic_annotations.tab
rm *paste*
rm *parse*

tr ' ' '\t' < likely_pathogenic_annotations.tab > likely_pathogenic_annotations.tab2
mv likely_pathogenic_annotations.tab2 likely_pathogenic_annotations.tab

sed -i 's/'p.'//g' likely_pathogenic_annotations.tab
awk '{if ($3 !=".") print $0}' likely_pathogenic_annotations.tab > temp_anno.tab

tr ' ' '\t' < temp_anno.tab > likely_pathogenic_annotations.tab
rm temp_anno.tab

# convert GBA notations to manuscript notations:
cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/convert_gba_notations.py .
./convert_gba_notations.py

# combine AD and DP files 
cat carriers_AD_patho_noGBAP1.vcf.tab carriers_DP_select_noGBAP1.vcf.tab > temp_carriers_all_to-R.tab

# feeding in correct number of columns, weeding out any triallelic mistakes: 
awk '{print $1,$2,$3,$4,$5}' temp_carriers_all_to-R.tab > carriers_all_to-R.tab

# annotate the coverage files
cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/merge_AD_annotations.R .
R < merge_AD_annotations.R --no-save

rm temp*
rm carriers_all_to-R.tab

# delete exons 10 & 11 because we call these with Sanger sequencing. 
sed -i '/'exon10'/d' carriers_all_patho_noGBAP1_annotated.tab
sed -i '/'exon11'/d' carriers_all_patho_noGBAP1_annotated.tab
### NOTE: this step gets rid of variants that weren't annotated. Check them out manually. 
sed -i '/'NA'/d' carriers_all_patho_noGBAP1_annotated.tab



## count carriers:

cut -f5 carriers_all_patho_noGBAP1_annotated.tab > variant_count.ls # cut -f2 if you want the position instead

awk '{ 
         for ( i=1; i<=NF; i++ ) # loop over all fields/columns
            dict[$i]++;      	 # count occurrence in an array using the field value as index/key
     }
 END {                           # after processing all data
         for (key in dict)       # iterate over all array keys
             if(dict[key]>0)     # if the key occurred 
                 print key " : " dict[key]    # print counter and key
     }' variant_count.ls > temp.dict

rm variant_count.ls

tr ' ' '\t' < temp.dict > variant_counts.dict

# if a variant is called in more than 40% of samples, remove it.
# rationale: since this is a forced alignment VCF, there is a lot of garbage... some variants are very overcalled. 
# these need to be removed. 
awk -v N=$N '{if ($3>N*0.20) print $1}' variant_counts.dict > too_many_calls.ls
# it is good practice to look at too_many_calls.ls and see if these variants are actually common, or just poorly called. 

cat too_many_calls.ls  | while read line
do
   sed -i '/'$line'/d' carriers_all_patho_noGBAP1_annotated.tab # remove these variants from your file. 
done 

# split the ones without differential AD
grep "SPLIT" carriers_all_patho_noGBAP1_annotated.tab > carriers_DP_patho_noGBAP1_annotated.tab
sed -i '/'SPLIT'/d' carriers_all_patho_noGBAP1_annotated.tab

mv carriers_all_patho_noGBAP1_annotated.tab carriers_AD_patho_noGBAP1_annotated.tab

## some coded investigation of AD

## total coverage at this position, for this sample:
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8+$9}' carriers_AD_patho_noGBAP1_annotated.tab > temp_total.tab

## percentage of alternate allele call (we usually keep > 25%):
awk '{if ($9 > 0 && $10 > 0) {$(NF+1)=$9/$10; print;} \
else {$(NF+1)="NA"; print;}}' \
temp_total.tab > temp_all_info.tab

# make calls (BASED ON 25% ALT ALLELE FILTER & >15X DP):
# $10 == total coverage
# $11 == percentage of alternate allele calls
awk '{$(NF+1)=$5; \
if ($11 == "NA") {$(NF+1)="remove"; print;}
else if ($10 < 15) {$(NF+1)="low_coverage"; print;} \
else if ($11 < 0.05) {$(NF+1)="remove"; print;} \
else if ($11 < 0.10) {$(NF+1)="poor_quality"; print;} \
else if ($11 < 0.25) {$(NF+1)="het_low_confidence"; print;} \
else if ($11 >= 0.90) {$(NF+1)="hom"; print;} \
else if ($8 == 0 && $10 > 15) {$(NF+1)="hom"; print;} \
else if ($11 >= 0.25) {$(NF+1)="het"; print;} \
else {$(NF+1)="error"; print}}' temp_all_info.tab > temp_labeled_all_info.tab

# make separate table for all "removed" variants (for future reference):
grep remove temp_labeled_all_info.tab > all_removed1.tab
sed -i '/'remove'/d' temp_labeled_all_info.tab
 
 # make calls for DP (based on VCF calls):
 awk '{$(NF+1)=$5; \
if ($8>15 && $9=="0/1") {$(NF+1)="het_no-DP"; print;} \
else if ($8 > 15 && $9 == "1/1") {$(NF+1)="hom_no-DP"; print;} \
else if ($8 < 15 && $9 == "0/1") {$(NF+1)="het_LC_no-DP"; print;} \
else if ($8 < 15 && $9 == "1/1") {$(NF+1)="hom_LC_no-DP"; print;} \
else {$(NF+1)="uncertain_no-DP"; print;}}' carriers_DP_patho_noGBAP1_annotated.tab \
 > temp_labeled_DP_info.tab
 
grep uncertain temp_labeled_DP_info.tab > all_removed2.tab
cat all_removed1.tab all_removed2.tab > all_removed.tab

sed -i '/'uncertain'/d' temp_labeled_DP_info.tab

rm all_removed1.tab all_removed2.tab
 
# join calls with the variant notation: 
sed -i 's/\s/':'/12' temp_labeled_all_info.tab
sed -i 's/\s/':'/10' temp_labeled_DP_info.tab
 
awk '{print $1,$2,$3,$4,$5,$6,$7,$10,$12}' temp_labeled_all_info.tab > temp_join1.tab
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10}' temp_labeled_DP_info.tab > temp_join2.tab

cat temp_join1.tab temp_join2.tab > temp_all_info_filtered_AD-DP.tab

# add header
sed -i '1i\CHR POS EXON FUNCTION NOTATION1 NOTATION2 SAMPLE A1_COV A2_COV TOTAL_COV PERCENT_ALT REC_CALL' temp_labeled_all_info.tab
sed -i '1i\CHR POS EXON FUNCTION NOTATION1 NOTATION2 SAMPLE TOTAL_COV REC_CALL' temp_all_info_filtered_AD-DP.tab

# make it tab delimited because I like it
tr ' ' '\t' < temp_labeled_all_info.tab > all_info_filtered_AD_noGBAP1.tab
tr ' ' '\t' < temp_all_info_filtered_AD-DP.tab > temp_all_info_filtered_AD-DP.tab

# remove synonymous variants.
sed -e "/\<synonymous\>/d" < temp_all_info_filtered_AD-DP.tab > all_info_filtered_AD-DP.tab

cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/final_output_perSample.py .
./final_output_perSample.py

## FILE TO SAVE : gba_pipeline_finalOutput.csv

rm temp*

# add header to original annotated coverage file
sed -i '1i\CHR\tPOS\tEXON\tFUNCTION\tNOTATION1\tNOTATION2\tSAMPLE\tA1_COV\tA2_COV' carriers_AD_patho_noGBAP1_annotated.tab

## Clean up 
rm carriers_AD_patho_noGBAP1.vcf.tab
rm AD_patho_noGBAP1.vcf.tab

mkdir COVERAGE_FILES
mv *AD* COVERAGE_FILES

mkdir ANNOTATIONS
mv *annotations* ANNOTATIONS

mkdir VARIANT_LISTS
mv *.ls VARIANT_LISTS

rm input_annovar.txt
rm *.out

rm *.py
rm *.R


rm -r ${NEW_NAME}_ANNOVAR_analysis_files_*

