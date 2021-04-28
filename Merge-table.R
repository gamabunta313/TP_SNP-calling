
files=list.files(getwd(),"tab")

read=function(files){
df=read.csv(files,"\t",header=TRUE,stringsAsFactor=FALSE)
df$Source=paste(files)
df$Source=gsub("_IARC-Annotated.tab","",df$Source)
df=df[,c(10,1:9)]
df}

m=do.call("rbind",lapply(files,read))
colnames(m)=c("Source","ChR","POS","REF","ALT","RO","AO","TYPE","AA_CHANGE","Effect")
write.table(m,"Good_Mutations.tab",sep="\t",row.names=FALSE)
