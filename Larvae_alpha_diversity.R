library(vegan)
library(ggplot2)
library(GUniFrac)
library(patchwork)
library(emmeans)
library(lme4)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Ontogeny/paper1/")

map=read.table("inputs/Larvae_metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,1:17]
tax=read.table("outputs/Larvae_ASVs_noContanoOut.txt", header = T, sep = "\t", row.names = 1)[,19:24]
tax$Family=gsub("Alcanivoracaceae1", "Alcanivoracaceae", tax$Family)
map$stage=gsub("Larvae15d", "Larvae 15 days", map$stage)
map$stage=gsub("Larvae28d",  "Larvae 28 days", map$stage)
P6=c( "#FB6F92", "#29BF12", "#FBC174", "#F48C06", "#247ba0")
###rarefying
cnts=t(asv)
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, min(rowSums(cnts)))$otu.tab.rff
rowSums(asv.rar)

# test no rarefaction
#asv.rar=cnts

############################################################
##################### Alpha-diversity ######################
############################################################

alpha=as.data.frame(t(estimateR(asv.rar)))
alpha$Shannon=diversity(asv.rar, index = "shannon")
alpha$Stage=map$stage[match(rownames(alpha), rownames(map))]

#write.table(alpha, "outputs/alpha_table.txt", sep = "\t", row.names = T, quote = F)

stats= alpha  %>% group_by(Stage) %>% 
  dplyr::summarise(Shannon=mean(Shannon), Shannon_se= sd(Shannon), 
                   Chao1= mean(S.chao1), Chao1_se= sd(S.chao1) )

alpha  %>% group_by(Stage) %>% 
  dplyr::summarise( Shannon_se= sd(Shannon))
##################################################
##################### Stats ######################
##################################################
#Shannon
sha_model=lm(Shannon~ Stage, data = alpha)
shannon_model=as.data.frame(anova(sha_model)) 
summary(sha_model) 
anova(sha_model)

plot(sha_model)
hist(sha_model$residuals, main = "Residual Histogram")
shapiro.test(residuals(sha_model))

#pairwise 
sha_pairs = emmeans(sha_model, pairwise ~ Stage, weights = "proportional", adjust="fdr")
sha_pairs_df=as.data.frame(sha_pairs[["contrasts"]])
#write.table(sha_pairs_df, "outputs/anova_shannon_pairwise.txt", sep = "\t", row.names = T, quote = F)


#Chao1
cha_model=lm(S.chao1~ Stage, data = alpha)
chao_model=as.data.frame(anova(cha_model)) 
summary(cha_model) 

#pairwise 
cha_pairs = emmeans(cha_model, pairwise ~ Stage, weights = "proportional", adjust="none")
rbind(cha_pairs$contrasts, adjust="fdr")
cha_pairs_df=as.data.frame(cha_pairs[["contrasts"]])
#write.table(cha_pairs_df, "outputs/anova_chao1_pairwise.txt", sep = "\t", row.names = T, quote = F)


## QC models

res_sha <- resid(sha_model)
plot(fitted(sha_model), res_sha)
qqnorm(res_sha)
qqline(res_sha) 
plot(density(res_sha))

res_cha <- resid(cha_model)
plot(fitted(cha_model), res_cha)
qqnorm(res_cha)
qqline(res_cha) 
plot(density(res_cha))

# boxplots
alpha$Stage=factor(alpha$Stage, levels = c("Egg","Sperm","Larvae 15 days" ,"Larvae 28 days","Seawater"))

shan=ggplot(alpha, aes(x=Stage, y=Shannon, fill=Stage)) + 
  stat_boxplot(geom = "errorbar")  + 
  geom_boxplot(alpha = 1) +  
  scale_fill_manual(values=P6)  + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none") + 
  labs( y= "Shannon diversity", x="")  + 
  annotate(geom="text", x=1, y=0.5, label= "A") + 
  annotate(geom="text", x=2, y=0.5, label= "A") + 
  annotate(geom="text", x=3, y=0.5, label= "A") + 
  annotate(geom="text", x=4, y=0.5, label= "B") + 
  annotate(geom="text", x=5, y=0.5, label= "C") 

cha1=ggplot(alpha, aes(x=Stage, y=S.chao1, fill=Stage)) + 
  stat_boxplot(geom = "errorbar")  + 
  geom_boxplot(alpha = 1) +  
  scale_fill_manual(values=P6)  + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none") + 
  labs( y= "Chao1 estimated richness", x="")  +
  annotate(geom="text", x=1, y=0.5, label= "A") + 
  annotate(geom="text", x=2, y=0.5, label= "B") + 
  annotate(geom="text", x=3, y=0.5, label= "C") + 
  annotate(geom="text", x=4, y=0.5, label= "A") + 
  annotate(geom="text", x=5, y=0.5, label= "A") 

pdf("./outputs/Larvae_AlphaDiversiy.pdf", width=7,height=4, pointsize = 12)
shan+cha1+ plot_annotation(tag_levels = 'A')
dev.off()
