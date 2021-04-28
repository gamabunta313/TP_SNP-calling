## mutations object is generated in the Db_Gen_from-vcf.sh script ##
## database object is generated in the Db_Gen_from-vcf.sh script ##

library(Johnson)
library(fastmatch)

`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

args = commandArgs(trailingOnly=TRUE)

mutations=read.csv(args[1],sep="\t",header=FALSE)

mutations<-subset(mutations,mutations$AA_CHANGE!="p.P72R")
mutations<-subset(mutations,mutations$Effect!="silent")
mutations<-subset(mutations,mutations$TYPE=="snp")

mutations$RO=as.numeric(as.character(mutations$RO))
mutations$AO=as.numeric(as.character(mutations$AO))
mutations$VAF=((mutations$AO)/((mutations$RO) + (mutations$AO)))*100
mutations$REF=gsub("TRUE","T",mutations$REF)
mutations$ALT=gsub("TRUE","T",mutations$ALT)

#tabella_risultati=data.frame()

database=read.csv(args[2],sep="",header=FALSE)

tabella_risultati=mutations

for (i in 1:nrow(mutations)){
  for (n in 1:nrow(database)){
    stats=vector()
    spline_values=vector()
    medians=vector()
    SD=vector()  
    error=vector()
    result =(interaction(paste(database[n,1],database[n,2],database[n,3],sep="_")) 
             %fin% 
               interaction(paste(mutations[i,2],mutations[i,3],mutations[i,5],sep="_")))
    if(result == TRUE) 
    {row=paste(n) 
    error[i]=paste(i)    
    list_names=list()
    mut_pos=as.numeric(mutations[i,3])
    database_copia=database
    #print(i)
    counter=0
    x<-as.numeric(database_copia[n,4:(ncol(database_copia))]) ; x <- x[!is.na(x)] 
    theta<-as.numeric((sum(x == 0)/length(x))); N<-as.numeric(length(x))
    tabella_risultati[i,12]<-paste(theta)
    w<-x[x!=0]
    w2<-c(w,as.numeric(mutations[i,11]))

    w_transf<-RE.Johnson(w2)$transformed
    tabella_risultati[i,13]<-paste(w_transf[length(w_transf)])
	
	if ( theta > 0.9 ) 
	{w3=w}
	else
    {w3<-w_transf[-length(w_transf)]}
	
    tabella_risultati[i,14]<-paste(mean(w3)); tabella_risultati[i,15]<-paste(sd(w3))
    tabella_risultati[i,16]=paste(length(w3))
    tabella_risultati[i,17]=paste(as.numeric(tabella_risultati[i,14]) + (as.numeric(tabella_risultati[i,15]) * 2.75))
	  tabella_risultati[i,18]=if(tabella_risultati[i,13] > tabella_risultati[i,17]){paste("S")}else{paste("NS")}
    
    }
  }
}
colnames(tabella_risultati)=c("Sample_Name","ChR","POS","REF","ALT","RO","AO","TYPE","AA_CHANGE","Effect","VAF","theta","ALT_trans","Mean_trans","SD_trans","Non-0-observations","Distribution_Cutoff","Significance")

write.table(tabella_risultati,"/path/to/Distributions.tab",sep="\t",row.names=FALSE)

