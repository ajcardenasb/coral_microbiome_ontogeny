library(phyloseq)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
map$stage2=gsub("0-9", "", map$stage)
map$stage=factor(map$stage, levels = c("Egg","Sperm","Larvae15d" ,"Larvae17d", "Larvae28d","Seawater"))

asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:23]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,24:29]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)

P6=c("#FF928B","#B9C0DA", "#645617", "#B79D2A","#E1CE7A", "#2589BD")

# #compositional approach
# phy.t=microbiome::transform(phy, transform = "clr", target = "OTU", shift = 0, scale = 1)
# ord = ordinate(phy.t, method = "RDA", distance = "euclidean")
# plot_ordination(phy.t,ord, color = "stage") + 
#   geom_point(size = 4, alpha = 1)  + 
#   scale_colour_manual(values=P6) + ggtitle("") +  
#   theme_bw() + theme( legend.position = 'bottom')

## classic approach
phy.t=microbiome::transform(phy, transform = "compositional", target = "OTU", shift = 0, scale = 1)
ord = ordinate(phy.t, method = "NMDS", distance = "bray") # 0.07205362 
pdf("./outputs/Larvae_ordination.pdf", width=5,height=5, pointsize = 10)
plot_ordination(phy.t,ord, color = "stage") + 
  geom_point(size = 4, alpha = 1)  + 
  scale_colour_manual(values=P6) + ggtitle("") +  
  theme_bw() + theme( legend.position = 'bottom')+
  annotate(geom="text", x=-2.5, y=-3.2, label= "Stress: 0.072")
dev.off()


#####################################################
######## Stats on community composition ############
####################################################
otu.n=data.frame(t(otu_table(phy.t)))
otu.n$Stage=map$stage[match(rownames(otu.n), rownames(map))]

#betadisper
distance=vegdist(otu.n[,1:1294], method = "bray")
beta_phenotype=betadisper(distance, otu.n$Stage)
permutest(beta_phenotype, permutations = 999, pairwise = TRUE) 

##overal model
adonis=adonis(otu.n[,1:1294]~ otu.n$Stage , method = "bray" )
adonis_df=as.data.frame(adonis[["aov.tab"]])
write.table(adonis_df, "outputs/Larvae_overall_adonis.txt", sep = "\t", row.names = T, quote = F)

#all comparisons
all_pairWadonis_df=pairwise.adonis(otu.n[,1:1294], otu.n$Stage,  sim.method = "bray", p.adjust.m = "fdr", perm = 999) 
write.table(all_pairWadonis_df, "outputs/Larvae_pairwiseAdonis.txt", sep = "\t", row.names = F, quote = F)
s
