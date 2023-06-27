library(UpSetR)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)
map$stage=gsub("Larvae15d", "Larvae T0", map$stage)
map$stage=gsub("Larvae28d",  "Larvae T1", map$stage)

## aggregate to genus
fam.wid.agg=aggregate(asv, by = list(tax[, 6]), FUN =  sum)
rownames(fam.wid.agg)=fam.wid.agg$Group.1
fam.wid.agg=fam.wid.agg[,-1]


#### genus rule 3 ####
# determine presence in at least 3 samples per group
binary=fam.wid.agg %>% mutate(across(everything(), ~replace(., . >0 , 1)))
colnames(binary)=map$stage[match(colnames(binary),rownames(map))]
binary_sums=as.data.frame(t(apply(binary, 1, function(x) tapply(x, colnames(binary), sum))))
present_s1=binary_sums %>% mutate(across(everything(), ~replace(., . <3 , 0)))
present_s2=present_s1 %>% mutate(across(everything(), ~replace(., . >=3 , 1))) %>% filter_all( any_vars(. != 0))

upset(present_s2,  sets = c( "Egg","Sperm" ,"Larvae T0","Larvae T1","Seawater"), 
             order.by="freq",  point.size=4, sets.bar.color=col=c("#2EC4B6","#F48C06", "#FB6F92", "#29BF12", "#FFD000"))

pdf("outputs/upset_indic_temp_pver.pdf", height = 4, width = 4)
upset1 
dev.off()


##### genos no 3 rule ######

binary=fam.wid.agg %>% mutate(across(everything(), ~replace(., . >0 , 1)))
colnames(binary)=map$stage[match(colnames(binary),rownames(map))]
binary_sums=as.data.frame(t(apply(binary, 1, function(x) tapply(x, colnames(binary), sum))))
present_s2=binary_sums %>% mutate(across(everything(), ~replace(., . >=1 , 1))) %>% filter_all( any_vars(. != 0))

upset(present_s2,  sets = c( "Egg","Sperm" ,"Larvae T0","Larvae T1","Seawater"), 
      order.by="freq",  point.size=4, sets.bar.color=c("#2EC4B6","#F48C06", "#FB6F92", "#29BF12", "#FFD000"))


##### ASV ######

binary=asv %>% mutate(across(everything(), ~replace(., . >0 , 1)))
colnames(binary)=map$stage[match(colnames(binary),rownames(map))]
binary_sums=as.data.frame(t(apply(binary, 1, function(x) tapply(x, colnames(binary), sum))))
present_s1=binary_sums %>% mutate(across(everything(), ~replace(., . <3 , 0)))
present_s2=present_s1 %>% mutate(across(everything(), ~replace(., . >=3 , 1))) %>% filter_all( any_vars(. != 0))

upset(present_s2,  sets = c( "Egg","Sperm" ,"Larvae T0","Larvae T1","Seawater"), 
      order.by="freq",  point.size=4, sets.bar.color=col)

