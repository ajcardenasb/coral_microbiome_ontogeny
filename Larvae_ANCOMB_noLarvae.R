library(ANCOMBC)
library(phyloseq)
library(microbiome)

### ANCON analyses only without larvae, each stage vs all other samples
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
map$Egg=ifelse(map$stage == "Egg", "yes", "no")
map$Seawater=ifelse(map$stage == "Seawater", "yes", "no")
map$Sperm=ifelse(map$stage == "Sperm", "yes", "no")

asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:23]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,24:29]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)
phy.gen=aggregate_taxa(phy, "Genus")
phy.noL=subset_samples(phy.gen, !stage %like% "^Larvae" )

res1=ancombc(phy.noL,formula="Egg",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Egg",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res1_df=data.frame(res1[["res"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "Others", "Egg")
res1_sig$Genus=rownames(res1_sig)
res1_sig$Family=tax$Family[match(rownames(res1_sig),tax$Genus)]
res1_sig$Comparison="Egg vs all"
message("Number of total DA genera: ", nrow(res1_sig), "\nNumber of DA genera enriched in Egg: ", nrow(subset(res1_sig, Diff_more_abundant == "Egg" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res1_sig, Diff_more_abundant == "Others")))

res2=ancombc(phy.noL,formula="Sperm",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Sperm",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res2_df=data.frame(res2[["res"]])
colnames(res2_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "Others", "Sperm")
res2_sig$Genus=rownames(res2_sig)
res2_sig$Family=tax$Family[match(rownames(res2_sig),tax$Genus)]
res2_sig$Comparison="Sperm vs all"
message("Number of total DA genera: ", nrow(res2_sig), "\nNumber of DA genera enriched in Seprm: ", nrow(subset(res2_sig, Diff_more_abundant == "Sperm" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res2_sig, Diff_more_abundant == "Others")))

res3=ancombc(phy.noL,formula="Seawater",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Seawater",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res3_df=data.frame(res3[["res"]])
colnames(res3_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "Others", "Seawater")
res3_sig$Genus=rownames(res3_sig)
res3_sig$Family=tax$Family[match(rownames(res3_sig),tax$Genus)]
res3_sig$Comparison="Seawater vs all"
message("Number of total DA genera: ", nrow(res3_sig), "\nNumber of DA genera enriched in Seawater: ", nrow(subset(res3_sig, Diff_more_abundant == "Seawater" )),  "\nNumber of DA genera enriched in Others: ", nrow(subset(res3_sig, Diff_more_abundant == "Others")))

all=rbind(res1_sig, res2_sig, res3_sig)
write.table(all, "outputs/Larvae_ANCOM_Results3_noLarvae.txt",  quote = F, sep = "\t")


