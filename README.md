# filtering_GBA
Special quality control is necessary for GBA because of sequencing issues caused by GBAP1. This pipeline should be run on both the forced alignment GBA VCF (excluding the pseudogene) and the regular VCF, to compare.  

**If on Beluga: you can take the bash script and run it as is _in the same directory as your VCF_.**  
bash filtering_GBA.sh $VCF $NEW_NAME $N

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
(any with calls for *all 3* alleles are removed. All with more than 3 alleles are removed.)
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
## Annotate the coverage files.  

We'll use [table_annovar.pl](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/startup/) to create a parse-able table of annotations, then a lot of parsing to match them to our covarage files.  
  
**Creating input_annovar table (specifications for this on annovar website)**  

*For those with AD data:*  
* take chr, base pair, A1, A2:
````
grep -f carrier_variants_to-investigate_AD.ls ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f1,2,4,5 > almost_input_annovar.txt
````
* if triallelic, take the second allele:
````
awk '{print $1,$2,$2,$3,$4}' almost_input_annovar.txt > to-tr_almost_input_annovar.txt
tr ',' ' ' < to-tr_almost_input_annovar.txt > almost_input_annovar.txt
awk '{if ($6>0) print $1,$2,$3,$4,$6; else print $1,$2,$3,$4,$5}' almost_input_annovar.txt > input_annovar1.txt
````
*For those with only DP:*
* take chr, base pair, A1, A2:
````
grep -f carrier_variants_to-investigate_DP.ls ${NEW_NAME}_ANNOVAR_output.gatk.b37.RESULTS.vcf | cut -f1,2,4,5 > almost_input_annovar.txt
````
* if triallelic, take the first allele:  
**Note:** This is an assumption. Since we don't have differential coverage (AD), we don't know which allele was called. The majority of the time, it is the "first" A1, which is why we make this assumption, but all of these should be validated in other data or taken with lots of precaution. *These variants will be indicated in the final file with "no_AD"*
````
awk '{print $1,$2,$2,$3,$4}' almost_input_annovar.txt > to-tr_almost_input_annovar.txt
tr ',' ' ' < to-tr_almost_input_annovar.txt > almost_input_annovar.txt
awk '{print $1,$2,$3,$4,$5}' almost_input_annovar.txt > input_annovar2.txt
````

* Combine the files and run table_annovar.pl:
````
cat input_annovar1.txt input_annovar2.txt > input_annovar.txt

cat input_annovar1.txt input_annovar2.txt > input_annovar.txt
bash ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/table_annovar.sh annotations_patho
````

* parse to get relevant information in clean table:
````
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
````
* remove "p." from protein change and exclude any which weren't annotated:
````
* remove "p." from protein change and exclude any which weren't annotated:
sed -i 's/'p.'//g' likely_pathogenic_annotations.tab
awk '{if ($3 !=".") print $0}' likely_pathogenic_annotations.tab > temp_anno.tab

tr ' ' '\t' < temp_anno.tab > likely_pathogenic_annotations.tab
rm temp_anno.tab
````

* convert GBA notations to manuscript notations:
````
cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/convert_gba_notations.py .
./convert_gba_notations.py
````
**The notation conversion python script for reference:**
```python
import pandas as pd

def convert(list):

    # Converting integer list to string list
    s = [str(i) for i in list]

    # Join list items using join()
    res = "".join(s)

    return(res)


def convert_gba_notation(variants):

    set_of_variants = set(variants)
    new_dict = {}

    for index, element in enumerate(set_of_variants):
        letters = []
        number = []
        garbage = []
        for i in element:
            if not i.isdecimal():
                letters.append(i)
            elif i.isdecimal():
                number.append(i)
            else: garbage.append(i)
        number = int(convert(number))
        if number > 39:
            new_number = number-39
        else:
            new_number = number-40

        new = [letters[0], str(new_number), letters[1]]
        new = convert(new)

        new_dict[element] = new


    df = pd.DataFrame(new_dict.items(), columns = ["PROTEIN_CHANGE", "NEW"])
    return df


df = pd.read_csv("likely_pathogenic_annotations.tab", delimiter="\t")

df.columns = ['Chr', 'Start', 'ExonicFunc', 'Exon', 'Protein_Change']

list_of_pros = df['Protein_Change'].values.tolist()

conversion = convert_gba_notation(list_of_pros)

conversion.to_csv(r'notation_conversion.tab', index = False, sep="\t")
````

* combine AD and DP files:
````
cat carriers_AD_patho_noGBAP1.vcf.tab carriers_DP_select_noGBAP1.vcf.tab > temp_carriers_all_to-R.tab
awk '{print $1,$2,$3,$4,$5}' temp_carriers_all_to-R.tab > carriers_all_to-R.tab
````

* annotate the coverage files:
````
cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/merge_AD_annotations.R .
R < merge_AD_annotations.R --no-save

rm temp*
rm carriers_all_to-R.tab
````

**The annotation R script for reference:**
```R
require(data.table)

data1 <- fread("likely_pathogenic_annotations.tab")
data2 <- fread("carriers_all_to-R.tab")
data3 <- fread("notation_conversion.tab")

colnames(data1) <- c("CHR", "POS", "FUNCTION", "EXON", "PROTEIN_CHANGE")
colnames(data2) <- c("chrom", "POS", "S_NUM", "COV_A1", "COV_A2")
colnames(data3) <- c("PRO_CHANGE", "NEW")


merge1 <- merge(data1, data3, by.x="PROTEIN_CHANGE", by.y="PRO_CHANGE", all.x=TRUE)
mdata <- merge1[,c("CHR", "POS", "FUNCTION", "EXON", "PROTEIN_CHANGE", "NEW")]
merged <- merge(mdata, data2, by="POS",all.y=TRUE)
dat <- merged[,c("chrom", "POS", "EXON", "FUNCTION", "NEW", "PROTEIN_CHANGE", "S_NUM", "COV_A1", "COV_A2")]

write.table(dat, file="carriers_all_patho_noGBAP1_annotated.tab", col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\t")

q()
```
## Quality Control
Now we have an annotated coverage file for possibly pathogenic variants. Gotta clean it up.  
  
* remove exons 10 & 11 because we call these with Sanger sequencing:
```unix
sed -i '/'exon10'/d' carriers_all_patho_noGBAP1_annotated.tab
sed -i '/'exon11'/d' carriers_all_patho_noGBAP1_annotated.tab
```
* exclude any which weren't annotated:
```unix
sed -i '/'NA'/d' carriers_all_patho_noGBAP1_annotated.tab
````

* count carriers of each variant & create a dictionary:
````
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
````
* remove variants which were called in more than 20% of samples:  
*This is a particular problem with the forced alignment VCF. Variants which are overcalled are generally poor quality.*
````
awk -v N=$N '{if ($3>N*0.20) print $1}' variant_counts.dict > too_many_calls.ls

cat too_many_calls.ls  | while read line
do
   sed -i '/'$line'/d' carriers_all_patho_noGBAP1_annotated.tab # remove these variants from your file. 
done 
````
* split the ones without differential AD:
````
grep "SPLIT" carriers_all_patho_noGBAP1_annotated.tab > carriers_DP_patho_noGBAP1_annotated.tab
sed -i '/'SPLIT'/d' carriers_all_patho_noGBAP1_annotated.tab

mv carriers_all_patho_noGBAP1_annotated.tab carriers_AD_patho_noGBAP1_annotated.tab
````

## Calling Variants  
Call | Criteria
-----|---------
Remove | No alternate allele call; alternate allele <5% total calls; no call (DP data)
Low coverage (LC) | <15 total coverage
Poor quality | 5% < alt allele < 10% total calls
Het: low confidence | 10% < alt allele < 25%
Het | alt allele >= 25% total calls
Hom | alt allele >= 90% total calls 

**AD Samples**
* calculate total DP:
````
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8+$9}' carriers_AD_patho_noGBAP1_annotated.tab > temp_total.tab
````

* calculate percentage of calls for the alternate allele:
````
awk '{if ($9 > 0 && $10 > 0) {$(NF+1)=$9/$10; print;} \
else {$(NF+1)="NA"; print;}}' \
temp_total.tab > temp_all_info.tab
````

* make calls ($10 == total coverage, $11==percentage of alternate allele calls)

````
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
````

* make separate table for all "removed" variants (for future reference) & remove them from main file:
````
grep remove temp_labeled_all_info.tab > all_removed1.tab
sed -i '/'remove'/d' temp_labeled_all_info.tab
````
**DP Samples**

* make calls for DP (based on VCF calls):
````
 awk '{$(NF+1)=$5; \
if ($8>15 && $9=="0/1") {$(NF+1)="het_no-DP"; print;} \
else if ($8 > 15 && $9 == "1/1") {$(NF+1)="hom_no-DP"; print;} \
else if ($8 < 15 && $9 == "0/1") {$(NF+1)="het_LC_no-DP"; print;} \
else if ($8 < 15 && $9 == "1/1") {$(NF+1)="hom_LC_no-DP"; print;} \
else {$(NF+1)="uncertain_no-DP"; print;}}' carriers_DP_patho_noGBAP1_annotated.tab \
 > temp_labeled_DP_info.tab
 ````
* make separate table for all removed variants & remove them from main file:
````
grep uncertain temp_labeled_DP_info.tab > all_removed2.tab
cat all_removed1.tab all_removed2.tab > all_removed.tab

sed -i '/'uncertain'/d' temp_labeled_DP_info.tab

rm all_removed1.tab all_removed2.tab
````
 
* join calls with the variant notation: 
````
sed -i 's/\s/':'/12' temp_labeled_all_info.tab
sed -i 's/\s/':'/10' temp_labeled_DP_info.tab
````

## Finalize files
* combine AD&DP files:
````
awk '{print $1,$2,$3,$4,$5,$6,$7,$10,$12}' temp_labeled_all_info.tab > temp_join1.tab
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10}' temp_labeled_DP_info.tab > temp_join2.tab

cat temp_join1.tab temp_join2.tab > temp_all_info_filtered_AD-DP.tab
````

* add headers:
````
sed -i '1i\CHR POS EXON FUNCTION NOTATION1 NOTATION2 SAMPLE A1_COV A2_COV TOTAL_COV PERCENT_ALT REC_CALL' temp_labeled_all_info.tab
sed -i '1i\CHR POS EXON FUNCTION NOTATION1 NOTATION2 SAMPLE TOTAL_COV REC_CALL' temp_all_info_filtered_AD-DP.tab
````

* make it tab delimited because I like it:
````
tr ' ' '\t' < temp_labeled_all_info.tab > all_info_filtered_AD_noGBAP1.tab
tr ' ' '\t' < temp_all_info_filtered_AD-DP.tab > temp_all_info_filtered_AD-DP.tab
````
* remove synonymous variants (somehow they creep in, will resolve later):
````
sed -e "/\<synonymous\>/d" < temp_all_info_filtered_AD-DP.tab > all_info_filtered_AD-DP.tab
````

## Make your CSV with calls per sample!
````
cp ~/runs/go_lab/mips/unfiltered/forced_alingment_gba/CALLING/final_output_perSample.py .
./final_output_perSample.py
````
**Python script for reference:**
```python
import pandas as pd

pdata = pd.read_csv("all_info_filtered_AD-DP.tab", delimiter="\t")

df = pdata

list_of_names = df['SAMPLE'].values.tolist()
list_of_calls = df['REC_CALL'].values.tolist()


d = {}

for name, pro in zip(list_of_names, list_of_calls):
    if name not in d:
        d[name] = pro
    elif type(d[name]) == list:
        d[name].append(pro)
    else:
        d[name] = [d[name], pro]


df_dict = pd.DataFrame(d.items(), columns = ["Sample", "GBA_variants"])

print(df_dict.head())

df_dict.to_csv(r'gba_pipeline_finalOutput.csv', index = False)
````

**FILE TO SAVE** : gba_pipeline_finalOutput.csv

* add header to original annotated coverage file:
````
sed -i '1i\CHR\tPOS\tEXON\tFUNCTION\tNOTATION1\tNOTATION2\tSAMPLE\tA1_COV\tA2_COV' carriers_AD_patho_noGBAP1_annotated.tab
````

**Clean up:** 
````
rm temp*
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
````
