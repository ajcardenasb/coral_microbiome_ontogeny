library(VennDiagram)
library(gridExtra)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:23]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,24:29]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

## aggregate to genus
fam.wid.agg=aggregate(asv, by = list(tax[, 6]), FUN =  sum)
all.l=melt(fam.wid.agg, id.vars=c("Group.1"), variable.name = "Genus", value.name = "Abundance")
colnames(all.l)=c("Genus","Sample","Abundance")
present=subset(all.l, !Abundance == 0) ## can be any cutoff!!
present$Stage=map$stage[match(present$Sample, rownames(map))]

list_egg=subset(present, Stage == "Egg")$Genus
list_sper=subset(present, Stage == "Sperm")$Genus
list_lar15=subset(present, Stage == "Larvae15d")$Genus
list_lar17=subset(present, Stage == "Larvae17d")$Genus
list_lar28=subset(present, Stage == "Larvae28d")$Genus
list_lar=subset(present, Stage %like% "^Larvae")$Genus
list_sea=subset(present, Stage == "Seawater")$Genus

Ven=venn.diagram(list(list_egg, list_sper, list_lar, list_sea ), NULL, 
                category.names = c("Egg","Sperm", "Larvae", "Seawater"), 
                fill = c("#FF928B","#B9C0DA","#B79D2A", "#2589BD") )

Ven_la=venn.diagram(list(list_lar15, list_lar17,list_lar28), NULL, 
                 category.names = c("Larvae15d", "Larvae17d", "Larvae28d"), 
                 fill = c("#645617", "#B79D2A","#E1CE7A") )

pdf("./outputs/Larvae_VennDiagram.pdf", width=5,height=3, pointsize = 10)
grid.arrange( grobTree(Ven),grobTree(Ven_la), ncol = 2, nrow =1)
dev.off()
