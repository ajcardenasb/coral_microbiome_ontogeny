library(phyloseq)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Ontogeny/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
map$stage2=gsub("0-9", "", map$stage)
map$stage2=gsub("Larvae15d", "Larvae 15 days", map$stage2)
map$stage2=gsub("Larvae28d",  "Larvae 28 days", map$stage2)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

map$stage2=factor(map$stage2, levels = c("Egg","Sperm","Larvae 15 days" , "Larvae 28 days","Seawater"))


otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)

#P6=c( "#f25f5c", "#f6ae2d", "#b0ff5c", "#4a8f00" , "#247ba0")
P6=c( "#FB6F92", "#29BF12", "#FBC174", "#F48C06", "#247ba0")
#P6=c("#2EC4B6","#F48C06", "#FB6F92", "#29BF12", "#FFD000") #Valentine's version


#compositional approach
phy.t=microbiome::transform(phy, transform = "clr", target = "OTU", shift = 0, scale = 1)
ord = ordinate(phy.t, method = "RDA", distance = "euclidean")
pdf("./outputs/Larvae_ordination.pdf", width=4,height=4, pointsize = 10)
plot_ordination(phy.t,ord, color = "stage2") +
  geom_point(size = 3, alpha = 1)  +
  scale_colour_manual(values=P6) + ggtitle("") +
  theme_bw() + theme( legend.position = c(0.80,0.25), legend.title=element_blank()) 
dev.off()


## classic approach
# phy.t=microbiome::transform(phy, transform = "compositional", target = "OTU", shift = 0, scale = 1)
# ord = ordinate(phy.t, method = "NMDS", distance = "bray") # 0.047
# #pdf("./outputs/Larvae_ordination.pdf", width=4,height=4, pointsize = 10)
# plot_ordination(phy.t,ord, color = "stage2") +
#   geom_point(size = 4, alpha = 1)  +
#   scale_colour_manual(values=P6) + ggtitle("") +
#   theme_bw() + theme( legend.position = 'bottom')+
#   annotate(geom="text", x=-2.0, y=-3.2, label= "Stress: 0.047")# +
#dev.off()

#### dissimilarities #####

# names(asv)
# asv_clr=apply(asv,2,clr)
# dist=vegdist(t(asv_clr), method="euclidean", upper = T)
# dist_df=as.matrix(dist)
# 
# 
# long_dist = reshape2::melt(dist_df, value.name = c("Distance"))#### Melting Data into Long Format
# long_dist$Stage1=map$stage2[match(long_dist$Var1, rownames(map))]
# long_dist$Stage2=map$stage2[match(long_dist$Var2, rownames(map))]
# long_dist$Comparison=paste(long_dist$Stage1, " vs ",long_dist$Stage2, sep = "")
# 
# long_s=subset(long_dist, Stage1 %in% c("Larvae 14 days" , "Larvae 28 days") & 
#                                              Stage2 %in% c("Egg", "Sperm", "Seawater"))
# long_s$Stage2=factor(long_s$Stage2, levels = c("Egg", "Sperm", "Seawater"))
# 
# #Effect of season on depth per species
# P3=c( "#FB6F92", "#29BF12", "#247ba0")
# #pdf("Outputs/dissimilarities.pdf", height = 3, width = 4, pointsize = 12)
# disi_plot=ggplot(long_s, aes(x=Stage2, y=Distance, fill= Stage2)) + 
#   geom_boxplot()  +
#   labs(y = "Dissimilarity to\nlarval microbiomes", x= "") + theme_minimal() +
#   theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   scale_fill_manual(values=P3)  + theme( legend.position = 'none') #facet_grid(~NucleicAcid1) + 
# dev.off()
# 
# pdf("outputs/CoralDevelopment_beta_diversity.pdf", height = 3.5, width = 5, pointsize = 12)
# ordi_plot + disi_plot + plot_annotation(tag_levels = "A")
# dev.off()

#####################################################
######## Stats on community composition ############
####################################################
otu.n=data.frame(t(otu_table(phy.t)))
otu.n$Stage=map$stage[match(rownames(otu.n), rownames(map))]

#betadisper
distance=vegdist(otu.n[,1:1490], method = "euclidean")
beta_phenotype=betadisper(distance, otu.n$Stage)
beta_disp_df=permutest(beta_phenotype, permutations = 999, pairwise = TRUE) 

##overal model
adonis=adonis2(otu.n[,1:1490]~ otu.n$Stage , method = "euclidean" )
write.table(adonis, "outputs/Larvae_overall_adonis.txt", sep = "\t", row.names = T, quote = F)

#all comparisons
all_pairWadonis_df=pairwise.adonis(otu.n[,1:1490], otu.n$Stage,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 9999) 
write.table(all_pairWadonis_df, "outputs/Larvae_pairwiseAdonis.txt", sep = "\t", row.names = F, quote = F)
