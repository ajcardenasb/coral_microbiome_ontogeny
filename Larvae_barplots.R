library(reshape2)
library(ggplot2)
library(scales)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:23]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,24:29]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)

## identify top 20 families 
names(asv)
fam.wid.agg=aggregate(asv, by = list(tax[, 6]), FUN =  sum)
topFamilies=fam.wid.agg[order(rowSums(fam.wid.agg[, 2:ncol(fam.wid.agg)]),decreasing = TRUE),][1:20,1]
fam.wid.agg$Group.1=ifelse(fam.wid.agg$Group.1 %in% topFamilies, as.character(fam.wid.agg$Group.1), "zOthers")
fa.gg2=aggregate(fam.wid.agg[, 2:ncol(fam.wid.agg)], by = list(fam.wid.agg[, 1]), FUN =  sum)
all.l=melt(fa.gg2, id.vars=c("Group.1"), variable.name = "Taxa", value.name = "Abundance")
colnames(all.l)=c("Taxa","Sample","Abundance")

## Add sample information
all.l$Stage=map$stage[match(all.l$Sample, rownames(map))]

final=all.l %>% group_by(Stage, Sample,Taxa) %>% dplyr::summarise(Abundance=sum(Abundance))
final$Stage=factor(final$Stage, levels = c("Egg","Sperm","Larvae15d" ,"Larvae17d", "Larvae28d","Seawater"))

## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
#P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")

#pdf("./outputs/Larvae_barplots_families.pdf",  width = 10.5, height =4, pointsize = 12) 
pdf("./outputs/Larvae_barplots_genus.pdf",  width = 10.5, height =4, pointsize = 12) 

ggplot() + geom_bar(aes(y = Abundance, x = Sample, fill = Taxa), 
                   data = final, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="") + 
  scale_fill_manual("Genus",values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'bottom', axis.text.x=element_text(angle=90,hjust=1),) + 
  guides(fill=guide_legend(nrow=3)) + facet_wrap(~Stage, scales = "free", ncol = 6) + scale_y_continuous( expand = c(0, 0)) 

dev.off() 
