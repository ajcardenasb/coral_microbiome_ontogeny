library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(data.table)
library(dplyr)

### ANCON analyses only without larvae, each stage vs all other samples
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
map$Egg=ifelse(map$stage == "Egg", "yes", "no")
map$Seawater=ifelse(map$stage == "Seawater", "yes", "no")
map$Sperm=ifelse(map$stage == "Sperm", "yes", "no")
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)
#phy.gen=aggregate_taxa(phy, "Genus")
phy.noL=subset_samples(phy, !stage %like% "^Larvae" )

res1=ancombc2(phy.noL, fix_formula ="Egg",p_adj_method = "fdr", group = NULL, tax_level = "Genus")
res1_df=data.frame(res1[["res"]])[,c(1,3,5,7,9,11,13)]
colnames(res1_df)=c(	"Taxon","LFG",	"LFC_se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "Others", "Egg")
res1_sig$Family=tax$Family[match(res1_sig$Taxon,tax$Genus)]
res1_sig$Order=tax$Order[match(res1_sig$Taxon,tax$Genus)]
res1_sig$Class=tax$Class[match(res1_sig$Taxon,tax$Genus)]
res1_sig$Phylum=tax$Phylum[match(res1_sig$Taxon,tax$Genus)]
message("Number of total DA genera: ", nrow(res1_sig), "\nNumber of DA genera enriched in Egg: ", nrow(subset(res1_sig, Diff_more_abundant == "Egg" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res1_sig, Diff_more_abundant == "Others")))

res2=ancombc2(phy.noL, fix_formula ="Sperm",p_adj_method = "fdr", group = NULL, tax_level = "Genus")
res2_df=data.frame(res2[["res"]])[,c(1,3,5,7,9,11,13)]
colnames(res2_df)=c(	"Taxon","LFG",	"LFC_se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "Others", "Sperm")
res2_sig$Family=tax$Family[match(res2_sig$Taxon,tax$Genus)]
res2_sig$Order=tax$Order[match(res2_sig$Taxon,tax$Genus)]
res2_sig$Class=tax$Class[match(res2_sig$Taxon,tax$Genus)]
res2_sig$Phylum=tax$Phylum[match(res2_sig$Taxon,tax$Genus)]
message("Number of total DA genera: ", nrow(res2_sig), "\nNumber of DA genera enriched in Sperm: ", nrow(subset(res2_sig, Diff_more_abundant == "Sperm" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res2_sig, Diff_more_abundant == "Others")))

res3=ancombc2(phy.noL, fix_formula ="Seawater",p_adj_method = "fdr", group = NULL, tax_level = "Genus")
res3_df=data.frame(res3[["res"]])[,c(1,3,5,7,9,11,13)]
colnames(res3_df)=c(	"Taxon","LFG",	"LFC_se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "Others", "Seawater")
res3_sig$Family=tax$Family[match(res3_sig$Taxon,tax$Genus)]
res3_sig$Order=tax$Order[match(res3_sig$Taxon,tax$Genus)]
res3_sig$Class=tax$Class[match(res3_sig$Taxon,tax$Genus)]
res3_sig$Phylum=tax$Phylum[match(res3_sig$Taxon,tax$Genus)]
message("Number of total DA genera: ", nrow(res3_sig), "\nNumber of DA genera enriched in Seawater: ", nrow(subset(res3_sig, Diff_more_abundant == "Seawater" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res3_sig, Diff_more_abundant == "Others")))

all=rbind(res1_sig, res2_sig, res3_sig) %>% filter(!Diff_more_abundant == "Others")

write.table(all, "outputs/Larvae_ANCOM_Results3_May2024.txt",  quote = F, sep = "\t")
