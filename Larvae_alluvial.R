library(ggalluvial)
library(reshape2)
library(ggplot2)
library(scales)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:23]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,24:29]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

## identify top 20 families 
asv.r=asv.r=apply(asv,2,function(x){x/sum(x)})
fam.wid.agg=aggregate(asv.r, by = list(tax[, 6]), FUN =  sum)
larvae15_samples=rownames(subset(map, stage == "Larvae15d"))
topGenusLarvae15d=fam.wid.agg[order(rowSums(fam.wid.agg[, names(fam.wid.agg) %in% larvae15_samples]),decreasing = TRUE),][1:10,1]
fam.wid.agg$Group.1=ifelse(fam.wid.agg$Group.1 %in% topGenusLarvae15d, as.character(fam.wid.agg$Group.1), "zOthers")
fa.gg2=aggregate(fam.wid.agg[, 2:ncol(fam.wid.agg)], by = list(fam.wid.agg[, 1]), FUN =  sum)
rownames(fa.gg2)=fa.gg2$Group.1
fa.gg2=fa.gg2[,-1]
colnames(fa.gg2)=map$stage[match(colnames(fa.gg2),rownames(map))]
means=as.data.frame(t(apply(fa.gg2, 1,function(x) tapply(x, colnames(fa.gg2), mean))))
means$Taxa=rownames(means)
means.l=melt(means, id.vars=c("Taxa"), variable.name = "Stage", value.name = "Abundance")

final=means.l %>% group_by(Stage, Taxa) %>% dplyr::summarise(Abundance=sum(Abundance)) %>% filter(Stage %in% c("Egg","Sperm","Larvae15d" ,"Larvae28d" ,"Seawater"))
final$Stage=factor(final$Stage, levels = c("Egg","Sperm","Seawater","Larvae15d" ,"Larvae28d"))


# sort taxa by source
anc=read.table( "outputs/Larvae_ANCOM_Results3_noLarvae.txt", sep = "\t", header = T)
anc_rele=subset(anc, Genus %in% topGenusLarvae15d & !Diff_more_abundant == "Others" )[,c(7,8)]
as.vector(anc_rele$Genus)

final$Taxa=factor(final$Taxa, levels = c("Pseudophaeobacter","Ruegeria", "Algicola", "Spongiibacter","Alcanivorax","Oleibacter","Vibrio","Alteromonas", "Marinobacter", "Peredibacter","zOthers"))
source=subset(final, Stage %in% c("Egg","Sperm","Seawater"))
larvae=subset(final, !Stage %in% c("Egg","Sperm","Seawater"))


## Plot
P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#961658", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")
sou_plot=ggplot() + geom_bar(aes(y = Abundance, x = Stage, fill = Taxa), 
                    data = source, stat="identity", position = "fill")  +
  labs( y= "", x="") + 
  scale_fill_manual("Genus",values=P10)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'none', ) + 
  guides(fill=guide_legend(nrow=2)) + scale_y_continuous( expand = c(0, 0)) + facet_wrap(~Stage, ncol=1, scales = "free")

larv_plot=ggplot() + geom_bar(aes(y = Abundance, x = Stage, fill = Taxa), 
                    data = larvae, stat="identity", position = "fill")  +
  labs( y= "", x="") + 
  scale_fill_manual("Genus",values=P10)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'bottom') + 
  guides(fill=guide_legend(nrow=2)) + scale_y_continuous( expand = c(0, 0)) + facet_wrap(~Stage, ncol=2, scales = "free_x")

pdf("outputs/combined_barplots.pdf", height = 7, width = 7)
sou_plot+larv_plot
dev.off()
