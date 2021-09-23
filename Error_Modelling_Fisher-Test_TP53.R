## mutations object is generated in the Db_Gen_from-vcf.sh script ##
## tabella object is generated in the Db_Gen_from-vcf.sh script ##


library(dplyr)
library(fastmatch)

`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}

args = commandArgs(trailingOnly=TRUE)

mutations=read.csv(args[1],sep="\t",header=TRUE)
#colnames(mutations)=c("ChR","POS","REF","ALT","RO","AO","TYPE","AA_CHANGE","Effect","Sample_Name")

mutations<-subset(mutations,mutations$AA_CHANGE!="p.P72R")
mutations<-subset(mutations,mutations$Effect!="silent")
mutations<-subset(mutations,mutations$TYPE=="snp")

mutations$RO=as.numeric(as.character(mutations$RO))
mutations$AO=as.numeric(as.character(mutations$AO))
mutations$VAF=((mutations$AO)/((mutations$RO) + (mutations$AO)))*100
mutations$REF=gsub("TRUE","T",mutations$REF)
mutations$ALT=gsub("TRUE","T",mutations$ALT)

tabella=read.csv(args[2],sep="",header=FALSE)

tabella_fisher=mutations
params=data.frame()

for (i in 1:nrow(mutations)){
  for (n in 1:nrow(tabella)){

	result1 =(interaction(paste(tabella[n,1],tabella[n,2],tabella[n,3],sep="_")) 
              %fin% 
                interaction(paste(mutations[i,2],mutations[i,3],mutations[i,4],sep="_")))
    if(result1 == TRUE){
      database_copia=tabella
      tabella_fisher[i,12]=paste(mean(as.numeric(database_copia[n,4:ncol(database_copia)])))
    } 
}}
for (i in 1:nrow(mutations)){
  for (n in 1:nrow(tabella)){
    result2 =(interaction(paste(tabella[n,1],tabella[n,2],tabella[n,3],sep="_")) 
              %fin% 
                interaction(paste(mutations[i,2],mutations[i,3],mutations[i,5],sep="_")))

    if(result2 == TRUE){
      database_copia=tabella
      tabella_fisher[i,13]=paste(mean(as.numeric(database_copia[n,4:ncol(database_copia)])))
    } 
    }}

for ( i in 1:nrow(mutations)){if(as.numeric(tabella_fisher[i,15])>0.01){tabella_fisher[i,16]=paste("NS")}else{tabella_fisher[i,16]=paste("S")}}
colnames(tabella_fisher)=c("Sample_Name","ChR","POS","REF","ALT","RO","AO","TYPE","AA_CHANGE","Effect","VAF","REF-Mean","ALT-Mean","Fisher_Exact_Test")
write.table(tabella_fisher,"Fisher_Exact_Test.tab",sep="\t",row.names=FALSE)
 

