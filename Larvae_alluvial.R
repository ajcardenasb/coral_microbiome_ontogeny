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


                            library(reshape2)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Ontogeny/paper1/")

sperm_list=read.table( "outputs/sperm_larvae_shared_ASVs", sep = "\t", header = T)
egg_list=read.table( "outputs/egg_larvae_shared_ASVs", sep = "\t", header = T)
all_list=read.table( "outputs/all_shared_ASVs", sep = "\t", header = T)

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)
map$stage=gsub("Larvae15d", "Larvae 15 days", map$stage)
map$stage=gsub("Larvae28d",  "Larvae 28 days", map$stage)

asv.rel=sweep(asv,2,colSums(asv),"/")

colnames(asv.rel)=map$stage[match(colnames(asv.rel),rownames(map))]
asv_rel_mean=as.data.frame(t(apply(asv.rel, 1, function(x) tapply(x, colnames(asv.rel), mean))))


###### alluvial barplots  ####
contri=asv_rel_mean
contri$source=ifelse(rownames(contri) %in% sperm_list$x, "Sperm",
                     ifelse(rownames(contri) %in% egg_list$x, "Egg", 
                            ifelse(rownames(contri) %in% egg_list$x, "Egg",
                                   ifelse(rownames(contri) %in% all_list$x, "Egg+Sperm", "Unknown"))))

contri_gg=aggregate(contri[,-6], by = list(contri[, 6]), FUN =  sum)
#contri_gg=contri %>% group_by(Genus, Group) %>% tally() %>% arrange(-n)
contri_l=reshape2::melt(contri_gg, id.vars=c("Group.1"), variable.name = "Group", value.name = "Abundance")
colnames(contri_l)=c("Source","Stage","Abundance")
contri_l$Stage=factor(contri_l$Stage, levels = c("Egg", "Sperm", "Larvae 15 days", "Larvae 28 days", "Seawater" ))
contri_l$Source=factor(contri_l$Source, levels = c("Egg", "Sperm", "Egg+Sperm", "Unknown" ))

## Add sample information
#all.l$Stage=map$stage[match(all.l$Sample, rownames(map))]
#final=all.l %>% group_by(Stage, Sample,Taxa) %>% dplyr::summarise(Abundance=sum(Abundance))


P3=c("#FB6F92", "#29BF12", "#A1915B","#ababab")
sou_plot=ggplot() + geom_bar(aes(y = Abundance, x = Stage, fill = Source),  
                    data = subset(contri_l, Stage %in% c("Egg", "Sperm")), stat="identity", position = "fill")  +
  labs( y= "", x="") + 
  scale_fill_manual("Source",values=P3)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'right', ) + 
  guides(fill=guide_legend(ncol=1)) + scale_y_continuous( expand = c(0, 0)) + facet_wrap(~Stage, ncol=1, scales = "free")

larv_plot=ggplot() + geom_bar(aes(y = Abundance, x = Stage, fill = Source), 
                              data = subset(contri_l, Stage %in% c("Larvae 15 days", "Larvae 28 days", "Seawater")), stat="identity", position = "fill")  +
  labs( y= "", x="") + 
  scale_fill_manual("Source",values=P3)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'right') + 
  guides(fill=guide_legend(ncol=1)) + scale_y_continuous( expand = c(0, 0)) + facet_wrap(~Stage, ncol=3, scales = "free_x")

pdf("outputs/combined_barplots_V2.pdf", height = 4.2, width = 7)
sou_plot+larv_plot + plot_layout(guides = "collect") #& 
dev.off()

####  Log Fold 2 change ########

lfc=subset(asv_rel_mean, rownames(asv_rel_mean) %in% c(sperm_list$x,egg_list$x, all_list$x))
lfc$Source=ifelse(rownames(lfc) %in% egg_list$x, "Egg", ifelse(rownames(lfc) %in% sperm_list$x, "Sperm", "Egg+Sperm"))
lfc$LFC=ifelse(lfc$Source ==  "Egg", log2(lfc$`Larvae 28 days`/lfc$Egg), 
               ifelse(lfc$Source ==  "Sperm", log2(lfc$`Larvae 28 days`/lfc$Sperm), 
                      ave(log2(lfc$`Larvae 28 days`/lfc$Egg),log2(lfc$`Larvae 28 days`/lfc$Sperm) )))

#contri_l$Stage=factor(contri_l$Stage, levels = c("Egg", "Sperm", "Larvae 15 days", "Larvae 28 days", "Seawater" ))
#contri_l$Source=factor(contri_l$Source, levels = c("Egg", "Sperm", "Egg+Sperm", "Unknown" ))


lfc$abs_LFC=abs(lfc$LFC)
lfc$ASV=rownames(lfc)
lfc$Genus=tax$Genus[match(lfc$ASV, rownames(tax))] 
lfc$Genus=ifelse(is.na(lfc$Genus), paste("Unclassified", tax$Family, sep = " "),lfc$Genus ) #lfc$Genus == "NA"
lfc$taxa=paste(lfc$ASV,lfc$Genus)
lfc$enrichment=ifelse(lfc$LFC < 0, "Gametes", "Larvae")

lfc_2 = lfc %>% filter(!Genus %in% c("Cutibacterium", "Corynebacterium", "Lactobacillus", "Streptococcus", "Cutibacterium", "Cloacibacterium", "Gemella", "Staphylococcus", "Neisseria", "Haemophilus")) %>%
  group_by(Source, enrichment) %>% slice_max(order_by = abs_LFC, n = 5)



egg=ggplot(subset(lfc_2, Source =="Egg"), aes(x=reorder(taxa,LFC),y=LFC, label = taxa) ) + 
  geom_text(aes(y = LFC), vjust=0, size =1.8, hjust=0) + geom_bar(stat = "identity", aes(fill = Source)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#FB6F92")) + 
  ylim(-25,25) + labs(y="Log Fold Change", x="Egg",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

spe=ggplot(subset(lfc_2, Source =="Sperm"), aes(x=reorder(taxa,LFC),y=LFC, label = taxa) ) + 
  geom_text(aes(y = LFC), vjust=0, size =1.8, hjust=0) + geom_bar(stat = "identity", aes(fill = Source)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#29BF12")) + 
  ylim(-25,25)  + labs(y="Log Fold Change", x="Sperm",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

sea=ggplot(subset(lfc_2, Source =="Egg+Sperm"), aes(x=reorder(taxa,LFC),y=LFC, label = taxa) ) + 
  geom_text(aes(y = LFC), vjust=0, size =1.8, hjust=0) + geom_bar(stat = "identity", aes(fill = Source)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#A1915B")) + 
  ylim(-25,25) + labs(y="Log Fold Change", x="Egg+Sperm",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

pdf("outputs/Larvae_Log2FC.pdf", width = 3, height = 6, pointsize = 6)
egg/spe/sea
dev.off()


### exploration

sperm_list=read.table( "outputs/sperm_larvae_shared_ASVs", sep = "\t", header = T)
egg_list=read.table( "outputs/egg_larvae_shared_ASVs", sep = "\t", header = T)
all_list=read.table( "outputs/all_shared_ASVs", sep = "\t", header = T)

shared=subset(tax, rownames(tax) %in% c(sperm_list$x,sperm_list$x, all_list$x ))
shared$Source=ifelse(rownames(shared) %in% egg_list$x, "Egg", ifelse(rownames(shared) %in% sperm_list$x, "Sperm", "Egg+Sperm"))


