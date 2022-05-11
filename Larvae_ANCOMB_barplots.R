library(reshape2)
library(ggplot2)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

anc=read.table( "outputs/Larvae_ANCOM_Results3_noLarvae.txt", sep = "\t", header = T)
ind=read.table( "outputs/Larvae_indic_species.txt", sep = " ", header = T)
#comp=subset(anc, (Comparison == "Egg vs Larvae" & Diff_more_abundant == "Egg") | 
 #             (Comparison == "Larvae vs Sperm" & Diff_more_abundant == "Sperm") |
  #            (Comparison == "Larvae vs Seawater" & Diff_more_abundant == "Seawater"))
comp=subset(anc, !Diff_more_abundant == "Others")

comp$abs_W=abs(comp$W)
comp_l =comp %>%
  group_by(Comparison) %>%
  slice_max(order_by = abs_W, n = 20) %>%
  arrange(Comparison,abs_W)
comp_l$label = comp_l$Genus

#P6=c("#FF928B","#B9C0DA", "#645617", "#B79D2A","#E1CE7A", "#2589BD")

egg=ggplot(subset(comp_l, Diff_more_abundant =="Egg"), aes(x=reorder(label,abs_W),y=abs_W, label = label)) + 
  geom_text(aes(y = abs_W), vjust=0, size =3, hjust=0) + geom_bar(stat = "identity", 
                                                                  aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#FF928B")) + 
  ylim(0,25) + labs(y="Effect size", x="Differentially enriched genera",title = "") + theme_bw() + 
  theme(legend.position = 'bottom') 

spe=ggplot(subset(comp_l, Diff_more_abundant =="Sperm"), aes(x=reorder(label,abs_W),y=abs_W, label = label)) + 
  geom_text(aes(y = abs_W), vjust=0, size =3, hjust=0) + geom_bar(stat = "identity", 
                                                                  aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#B9C0DA")) + 
  ylim(0,25) + labs(y="Effect size", x="Differentially enriched genera",title = "") + theme_bw() + 
  theme(legend.position = 'bottom') 

sea=ggplot(subset(comp_l, Diff_more_abundant =="Seawater"), aes(x=reorder(label,abs_W),y=abs_W, label = label)) + 
  geom_text(aes(y = abs_W), vjust=0, size =3, hjust=0) + geom_bar(stat = "identity", 
                                                                  aes(fill = Diff_more_abundant)) + 
  coord_flip()  + scale_x_discrete(breaks = NULL) + 
  theme_bw() + scale_fill_manual("",values =c("#2589BD")) + 
  ylim(0,25) + labs(y="Effect size", x="Differentially enriched genera", title = "") + theme_bw() + 
  theme(legend.position = 'bottom') 

pdf("outputs/Larvae_biomarkers.pdf", width = 7, height = 4, pointsize = 6)
egg+spe+sea
dev.off()



