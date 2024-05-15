library(reshape2)
library(ggplot2)
library(patchwork)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

anc=read.table( "outputs/Larvae_ANCOM_Results3_May2024.txt", sep = "\t", header = T)

anc$abs_LFC=abs(anc$LFG)
comp_l =anc %>%
  group_by(Diff_more_abundant) %>% arrange(Diff_more_abundant,abs_LFC)


#P6=c("#FF928B","#B9C0DA", "#645617", "#B79D2A","#E1CE7A", "#2589BD")

egg=ggplot(subset(comp_l, Diff_more_abundant =="Egg"), aes(x=reorder(Taxon,abs_LFC),y=abs_LFC, label = Taxon) ) + 
  geom_text(aes(y = abs_LFC), vjust=0, size =2, hjust=0) + geom_bar(stat = "identity", aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#2EC4B6")) + 
  ylim(0,10) + labs(y="Log Fold Change", x="Egg",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

spe=ggplot(subset(comp_l, Diff_more_abundant =="Sperm"), aes(x=reorder(Taxon,abs_LFC),y=abs_LFC, label = Taxon) ) + 
  geom_text(aes(y = abs_LFC), vjust=0, size =2, hjust=0) + geom_bar(stat = "identity", aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#F48C06")) + 
  ylim(0,10) + labs(y="Log Fold Change", x="Sperm",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

sea=ggplot(subset(comp_l, Diff_more_abundant =="Seawater"), aes(x=reorder(Taxon,abs_LFC),y=abs_LFC, label = Taxon) ) + 
  geom_text(aes(y = abs_LFC), vjust=0, size =2, hjust=0) + geom_bar(stat = "identity", aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#FFD000")) + 
  ylim(0,10) + labs(y="Log Fold Change", x="Seawater",title = "") + theme_bw() + 
  theme(legend.position = 'none') 

pdf("outputs/Larvae_biomarkers_updated_May2024.pdf", width = 3, height = 8, pointsize = 6)
egg/spe/sea
dev.off()

pdf("outputs/Larvae_biomarkers_updated_May2024_V2.pdf", width = 4, height = 6, pointsize = 6)
ggplot(comp_l, aes(x=reorder(Taxon,abs_LFC),y=abs_LFC, label = Taxon) ) + 
  geom_text(aes(y = abs_LFC), vjust=0, size =3, hjust=0) + geom_bar(stat = "identity", aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#F48C06","#2EC4B6", "#FFD000")) + 
  ylim(0,10) + labs(y="Log Fold Change", x="",title = "") + theme_bw() + 
  theme(legend.position = 'none')  + facet_grid(Diff_more_abundant~., scales = "free", space = "free")
dev.off()
