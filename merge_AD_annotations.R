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
