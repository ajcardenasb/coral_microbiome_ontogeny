library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Ontogeny/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)
map$stage=gsub("Larvae15d", "Larvae 15 days", map$stage)
map$stage=gsub("Larvae28d",  "Larvae 28 days", map$stage)

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
final$Stage=factor(final$Stage, levels = c("Egg","Sperm","Larvae 15 days" , "Larvae 28 days","Seawater"))

## Plot
old_palette=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
P21 = createPalette(21,  c("#B33C86", "#197BBD", "#C46DF7"))

#pdf("./outputs/Larvae_barplots_families.pdf",  width = 10.5, height =4, pointsize = 12) 
pdf("./outputs/Larvae_barplots_genus.pdf",  width = 8.5, height =4, pointsize = 12) 

ggplot() + geom_bar(aes(y = Abundance, x = Sample, fill = Taxa), 
                   data = final, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="") + 
   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'bottom', axis.text.x=element_text(angle=90,hjust=1),) + 
  guides(fill=guide_legend(nrow=3)) + facet_wrap(~Stage, scales = "free", ncol = 6) + scale_y_continuous( expand = c(0, 0)) +
  scale_fill_manual( values=old_palette)

dev.off() 



#exporting relative abundances at genus and family level
gen_lev=aggregate(asv, by = list(tax[, 6]), FUN =  sum)
rownames(gen_lev)=gen_lev$Group.1
gen_lev=gen_lev[,-1]
gen.rel=sweep(gen_lev,2,colSums(gen_lev),"/")
write.table(gen.rel, "outputs/Genus_relAbun.txt", sep = "\t", quote = F, row.names = T)

fam_lev=aggregate(asv, by = list(tax[, 5]), FUN =  sum)
rownames(fam_lev)=fam_lev$Group.1
fam_lev=fam_lev[,-1]
fam.rel=sweep(fam_lev,2,colSums(fam_lev),"/")
write.table(fam.rel, "outputs/Family_relAbun.txt", sep = "\t", quote = F, row.names = T)


#Top genus per stage
colnames(gen.rel)=map$stage[match(colnames(gen.rel), rownames(map))]
colnames(gen.rel)=gsub(".[0-9]", "", colnames(gen.rel))
ave_gen=as.data.frame(t(apply(gen.rel, 1, function(x) tapply(x, colnames(gen.rel), mean))))
se_gen=as.data.frame(t(apply(gen.rel, 1, function(x) tapply(x, colnames(gen.rel), sd))))


rownames(arrange(ave_gen, Egg))[1:20]

ave_gen %>% arrange(desc(Egg))
ave_gen[order(Egg),]


top_egg=gen.rel[order(gen.rel$Egg,decreasing = TRUE),][1:20,1]
top_sperm=gen_agg[order(rowSums(gen_agg[,colnames(gen_agg) %in% rownames(subset(map, stage == "Sperm"))]),decreasing = TRUE),][1:20,1]
