library(UpSetR)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Ontogeny//paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)
map$stage=gsub("Larvae15d", "Larvae 15 days", map$stage)
map$stage=gsub("Larvae28d",  "Larvae 28 days", map$stage)


###### Analysis at the ASV level including seawater ############

binary=asv %>% mutate_if(is.numeric, ~1 * (. > 0))
colnames(binary)=map$stage[match(colnames(binary),rownames(map))]
binary_sums=as.data.frame(t(apply(binary, 1, function(x) tapply(x, colnames(binary), sum))))
present_s3=binary_sums %>% mutate_if(is.numeric, ~1 * (. > 0))

pdf("outputs/larvae_Upset.pdf", height = 3, width = 4)
upset(present_s3,  sets = c( "Egg","Sperm" ,"Larvae 15 days","Larvae 28 days","Seawater"), 
      order.by="freq",  point.size=3, sets.bar.color=c("#FB6F92", "#29BF12", "#FBC174", "#F48C06", "#247ba0"),
      empty.intersections = "on")
dev.off()

#export table
out_shared=present_s3
out_shared$o_Egg=ifelse(out_shared$Egg == 1, "Egg", "")
out_shared$o_Sperm=ifelse(out_shared$Sperm == 1, "Sperm", "")
out_shared$o_15_larvae=ifelse(out_shared$`Larvae 15 days` == 1, "15-day-old larvae", "")
out_shared$o_28_larvae=ifelse(out_shared$`Larvae 28 days` == 1, "28-day-old larvae", "")
out_shared$Shared=paste(out_shared$o_Egg, out_shared$o_Sperm,out_shared$o_15_larvae, out_shared$o_28_larvae, sep = "," )
out_shared$Shared=gsub("^,*|,*$", "",out_shared$Shared )

out_shared2=subset(out_shared, !Shared == "")
out_shared3=merge(out_shared2,tax,by="row.names",all.x=TRUE)

write.table(out_shared3[,-(2:10)], "outputs/shared_ASVs", row.names = F, quote = F, sep = "\t")

#### exploration

test=subset(out_shared3, Egg == 0 & Sperm == 1 & `Larvae 15 days` == 0 & `Larvae 28 days` == 0) %>% 
  group_by(Genus, Shared) %>% tally()


### barplots

sperm_list=rownames(subset(present_s3, Egg == 0 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 1 & Seawater=0))
egg_list=rownames(subset(present_s3, Egg == 1 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 0 & Seawater=0))

bar_upset= data.frame(ASV=c(sperm_list,egg_list))
bar_upset$Group=ifelse(bar_upset$ASV %in% sperm_list, "Sperm", "Egg")

bar_upset$Genus=tax$Genus[match(bar_upset$ASV, rownames(tax))]
bar_upset$Family=tax$Family[match(bar_upset$ASV, rownames(tax))]
bar_upset$Genus=ifelse(is.na(bar_upset$Genus), "Unclassified", bar_upset$Genus)
final=bar_upset %>% group_by(Genus, Group) %>% tally() %>% arrange(-n)
unique(final$Genus)
top_genera=unique(final$Genus)[c(2:11,13:18, 20,21,23,30,33,38,39)]
final$Genus=ifelse(final$Genus %in% top_genera, as.character(final$Genus), "zOthers")

final_v2=final%>% group_by(Genus, Group) %>% dplyr::summarise(n=sum(n))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711",  "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455","#ff7f0e", "#ffbb78" ,"#e377c2", "#f7b6d2","#C0C0C0", "#858585")
P10=c("#2077b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df89", "#17becf", "#9edae5", "#e377c2", "#C0C0C0", "#858585")

pdf("outputs/larvae_Upset_Comp_barplots_V2.pdf", height = 3.5, width = 3.5)
ggplot() + geom_bar(aes(y = n, x = Group, fill = Genus), 
                    data = final_v2, stat="identity", position = "stack")  +
  labs( y= "Number of taxa shared", x="") + theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'right') + guides(fill=guide_legend(ncol=1))  + scale_y_continuous( expand = c(0, 0)) +
  scale_fill_manual( values=P21)
dev.off()


###### Analysis at the ASV level excluding seawater ############

binary=asv %>% mutate_if(is.numeric, ~1 * (. > 0))
colnames(binary)=map$stage[match(colnames(binary),rownames(map))]
binary_sums=as.data.frame(t(apply(binary, 1, function(x) tapply(x, colnames(binary), sum))))
present_s3=binary_sums[,-4] %>% mutate_if(is.numeric, ~1 * (. > 0))

pdf("outputs/larvae_Upset_V2.pdf", height = 3, width = 4)
upset(present_s3,  sets = c( "Egg","Sperm" ,"Larvae 15 days","Larvae 28 days"), 
      order.by="freq",  point.size=3, sets.bar.color=c("#FB6F92", "#29BF12", "#FBC174", "#F48C06"),
      empty.intersections = "on")
dev.off()


### barplots

sperm_list=rownames(subset(present_s3, Egg == 0 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 1))
egg_list=rownames(subset(present_s3, Egg == 1 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 0))
all_list=rownames(subset(present_s3, Egg == 1 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 1))


bar_upset= data.frame(ASV=c(sperm_list,egg_list, all_list))
bar_upset$Group=ifelse(bar_upset$ASV %in% sperm_list, "Sperm", ifelse(bar_upset$ASV %in% egg_list,  "Egg", "All"))

bar_upset$Genus=tax$Genus[match(bar_upset$ASV, rownames(tax))]
bar_upset$Family=tax$Family[match(bar_upset$ASV, rownames(tax))]
bar_upset$Genus=ifelse(is.na(bar_upset$Genus), "Unclassified", bar_upset$Genus)
final=bar_upset %>% group_by(Genus, Group) %>% tally() %>% arrange(-n)
unique(final$Genus)
top_genera=unique(final$Genus)[c(1:7, 9,11:17, 19, 22,27:28, 41, 44, 45, 50, 52)]
final$Genus=ifelse(final$Genus == "Unclassified", "zUnclassified", 
                   ifelse(final$Genus %in% top_genera, as.character(final$Genus), "zzOthers"))

final_v2=final%>% group_by(Genus, Group) %>% dplyr::summarise(n=sum(n))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711",  "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455","#ff7f0e", "#ffbb78" ,"#e377c2", "#f7b6d2","#C0C0C0", "#858585")
P10=c("#2077b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df89", "#17becf", "#9edae5", "#e377c2", "#f7b6d2", "#DDDD77","#C0C0C0")

pdf("outputs/larvae_Upset_Comp_barplots_v2.pdf", height = 2, width = 3.5)
pdf("outputs/larvae_Upset_Comp_barplots_v3.pdf", height = 4, width = 4)
ggplot() + geom_bar(aes(y = n, x = Group, fill = Genus), 
                    data = final_v2, stat="identity", position = "stack")  +
  labs( y= "Number of taxa shared", x="") + theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'right') + guides(fill=guide_legend(ncol=1))  + scale_y_continuous( expand = c(0, 0)) +
  scale_fill_manual( values=P21)
dev.off()

sperm_list=rownames(subset(present_s3, Egg == 0 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 1))
egg_list=rownames(subset(present_s3, Egg == 1 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 0))
all_list=rownames(subset(present_s3, Egg == 1 & `Larvae 15 days` == 1 & `Larvae 28 days` == 1 & Sperm == 1))

write.table(sperm_list, "outputs/sperm_larvae_shared_ASVs", quote = F, row.names = F, sep = "\t")
write.table(egg_list, "outputs/egg_larvae_shared_ASVs", quote = F, row.names = F, sep = "\t")
write.table(all_list, "outputs/all_shared_ASVs", quote = F, row.names = F, sep = "\t")
#comparisons

explo_gen=bar_upset %>% group_by(Genus, Group) %>% tally() %>% arrange(-n)
explo_fam=bar_upset %>% group_by(Family, Group) %>% tally() %>% arrange(-n)





