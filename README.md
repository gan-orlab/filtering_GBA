# filtering_GBA
Special quality control is necessary for GBA because of sequencing issues caused by GBAP1. This pipeline should be run on both the forced alignment GBA VCF (excluding the pseudogene) and the regular VCF, to compare.

## Enter your variables. 
If you're starting from a forced alignment VCF, it will already include only GBA.  
If not, extract GBA in hg19:  
* --chr 1 --from-bp 155204139 --to-bp 155214753  
  
**VCF**: name of your VCF, including *only* GBA.  
**NEW_NAME**: whatever name you want affixed to your annotated VCF. **note**: it must start with a capital letter.  
**N**: number of samples in this batch.
````
VCF=$1 
NEW_NAME=$2 
N=$3 
````

## Annotate your VCF and selecting possible pathogenic variants.

````

bash ~/runs/krohn/krohn/Home/annovar/annovar_MINIMAL.sh $VCF $NEW_NAME 
# this script annotates with GENE, gnomad_genome, and clinvar

# extract likely pathogenic positions:
grep nonsynonymous ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > nonsynonymous_variants.ls
grep frameshift ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > frameshift_variants.ls
grep stopgain ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > stopgain_variants.ls
grep splicing ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f2 > splicing_variants.ls

cat frameshift_variants.ls nonsynonymous_variants.ls stopgain_variants.ls splicing_variants.ls > extract_possibly-pathogenic_variants.ls
````

## Parse the VCF to extract relevant information.
The most important piece of information here is the depth of coverage per allele (**AD**).   
Unfortunately, some sites don't have this information, and only have total depth of coverage (**DP**).   
For those sites, we extract the DP and the genotype call (**GT**). We don't have a way to filter these by
percentage of alternate allele call, so these sites should be treated skeptically. Ideally we can validate them with
another VCF from MIPs or with omnix GWAS data. 

````
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="AD") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > AD_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="DP") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > DP_all_noGBAP1.vcf.tab
awk -F ' ' '/^#CHROM/ {split($0,samples);next;} /^#/ {next;} {dpcol=-1;n=split($9,fmt,/\:/);for(i=1;i<=n;++i) if(fmt[i]=="GT") { dpcol=i;break;} if(dpcol==-1) next; for(i=10;i<=NF;++i) {split($i,a,/\:/); printf("%s\t%s\t%s\t%s\n",$1,$2,samples[i],a[dpcol]);}}' ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf > GT_all_noGBAP1.vcf.tab
````

## Prepare your coverage files.  
**First: those with DP only.**
* extract your list of pathogenic variants from the coverage file:
````
grep -f extract_possibly-pathogenic_variants.ls AD_all_noGBAP1.vcf.tab > to-tr_AD_patho_noGBAP1.vcf.tab
````
* see which variants don't have AD (we will use another approach):
````
cut -f2 to-tr_AD_patho_noGBAP1.vcf.tab | sort | uniq > vars_that_made_it.ls
sort extract_possibly-pathogenic_variants.ls | uniq > all_patho_vars.ls
comm -13 vars_that_made_it.ls all_patho_vars.ls > to-DP.ls
````
* select those with only DP and make a separate file:
````
grep -f to-DP.ls DP_all_noGBAP1.vcf.tab > DP_select_noGBAP1.vcf.tab
grep -f to-DP.ls GT_all_noGBAP1.vcf.tab | cut -f4 > paste_this
paste DP_select_noGBAP1.vcf.tab paste_this > DP_GT_noGBAP1.vcf.tab

rm vars_that_made_it.ls all_patho_vars.ls paste_this DP_select_noGBAP1.vcf.tab to-DP.ls
````

* extract only variant sites (remove wild type or no call)
````
awk '{if ($5 != "0/0" && $4 != "./.") print $1="SPLIT",$2,$3,$4,$5}' DP_GT_noGBAP1.vcf.tab > carriers_DP_select_noGBAP1.vcf.tab
````
**Now preparing the AD coverage files.**  

* separate the AD calls into two columns:
````
tr ',' '\t' < to-tr_AD_patho_noGBAP1.vcf.tab > temp_AD_patho_noGBAP1.vcf.tab
````
* extract tri-allelic positions that seem real:  
(any will calls for *all 3* alleles are removed. All with more than 3 alleles are removed.)
````
awk '{if ($6>0 && $5==0) print $1="TRI",$2,$3,$4,$6; else if ($6>0 && $5>0) print $1="TRI",$2,$3,$4,$5=0; else if ($6<0) print $1,$2,$3,$4,$5}' temp_AD_patho_noGBAP1.vcf.tab > AD_patho_noGBAP1.vcf.tab 

rm to-tr_AD_patho_noGBAP1.vcf.tab
rm temp*
````
* keep only carriers:
````
awk '{if ($5 > 0 && $4>=0) print $0}' AD_patho_noGBAP1.vcf.tab > carriers_AD_patho_noGBAP1.vcf.tab
````

* make a list of all positions you're investigating: 
````
cut -f2 -d ' ' carriers_AD_patho_noGBAP1.vcf.tab | sort | uniq > carrier_variants_to-investigate_AD.ls
cut -f2 -d ' ' carriers_DP_select_noGBAP1.vcf.tab | sort | uniq > carrier_variants_to-investigate_DP.ls
````
