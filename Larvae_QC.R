setwd("~/Documents/Bioinformatics_scripts/R_scripts/Valentine_et_al/paper1/")
library(stringr)
library(dplyr)


##################################################################
#####Identifying and removing contaminant ASVs normalized data#####
###################################################################

asv=read.table("inputs/Larvae_ASV_table.txt", header = TRUE, row.names = 1, sep = "\t")[,c(1:27,29:30)]#excluding negative control 5 as instructed by Valentine
tax=read.table("inputs/Larvae_ASV_table.txt", header = TRUE, row.names = 1, sep = "\t")[,32:37]
map=read.table("inputs/Larvae_metadata.txt", header = T, row.names = 1, sep = "\t")

# all samples have more than 10000 reads
hist(colSums(asv),  breaks = 50, labels = F)
axis(side=1, at=seq(0,150000, 10000))

asv.o=asv[, colSums(asv) > 10000]
message(ncol(asv.o)," samples with > 10000 reads were retained out of ", ncol(asv), " total samples")

#Identify and removing contaminant ASVs raw data
asv.r=as.data.frame(sweep(asv.o,2,colSums(asv.o),`/`))
asv.r$Sum=rowSums(asv.r[,1:ncol(asv.r)])
names(asv.r)
asv.r$NegSum=rowSums(asv.r[,24:29])
asv.r$contaFactor=(asv.r$NegSum/asv.r$Sum)*100
rownames(asv.r)=rownames(asv)
Conta=subset(asv.r, asv.r$contaFactor > 10)
Conta$Family=tax$Family[match(rownames(Conta), rownames(asv))]
message("Number of total ASVs: ", nrow(asv)) #111
message("Number of identified contaminant ASVs removed from the analysis: ", 
        length(rownames(Conta)), "\n", Conta$Family[1],"\n", 
        Conta$Family[2],"\n", Conta$Family[3],"\n", 
        Conta$Family[4],"\n", Conta$Family[5])

#remove any chloroplast or mitochobdria
unwant_tax=tax %>% filter_all(any_vars(str_detect(., 'Mitochondria|Chloroplast')))
colnames(asv.o)
asv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(Conta) & !rownames(asv.o) %in% rownames(unwant_tax))[,1:23]
colnames(asv.noConta)
dim(asv.noConta)

#remove ASVs with only zeros
asv.noRare=asv.noConta[rowSums(asv.noConta[])>0,]
dim(asv.noRare)
#rename sample name and export
colnames(asv.noRare)=gsub("_.*", "", colnames(asv.noRare))

# Export QC filtered ASV tables
asv.noConta.f=merge(asv.noRare, tax, by="row.names")

write.table(asv.noConta.f, "./outputs/Larvae_ASVs_noContanoOut.txt",  quote = FALSE, row.names=F, sep = "\t") 
message("Number of ASVs used in the analysis: ", length(rownames(asv.noConta.f)))
