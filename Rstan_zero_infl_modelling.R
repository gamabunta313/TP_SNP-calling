
library(fitdistrplus)
library(dplyr)
library(logspline)
library(tidyverse)
library(rstan)

args = commandArgs(trailingOnly=TRUE)

model_data <- list(
  n_obs = nrow(df), 
  y_shape = shape, 
  y_scale = scale, 
  non_zero = theta,   
  y=df$y
)

## is the Good_Mutations.tab generated from the Db_Gen_from-vcf.py script 
mutations=read.csv(args[1],sep="\t",header=FALSE,stringsAsFactors=FALSE)

mutations=mutations[which(mutations$`AA_CHANGE` != "p.P72R" & mutations$Effect != "silent" & mutations$TYPE == "snp"),]

# is the DB_VAF_Custom.tab generated from the Db_Gen_from-vcf.py script 
db=read.csv(args[2],sep="",header=FALSE,stringsAsFactors=TRUE)

# is the absolute path of the .stan model
model_name=args[3]

for (i in 1:nrow(mutations)){
    result=db[which( db$V1==mutations[i,2] & db$V2==mutations[i,3] & db$V3==mutations[i,5]),]
    x=as.numeric(result[,c(4:ncol(db))])
    theta<-(sum(x == 0)/length(x))
    n<-as.numeric(length(x))
    mutations[i,11]<-n
    mutations[i,12]<-paste(theta)
    if (theta<=0.3) {
        w<-x[x!=0]
          
        fit.weibull<-fitdist(w,"weibull")
          
        stats <- replicate(4000, {      
        r <- rweibull(n*(1-theta),fit.weibull$estimate["shape"],fit.weibull$estimate["scale"])
            
        estfit.weibull<-fitdist(r,"weibull")

        as.numeric(ks.test(r,"pweibull",estfit.weibull$estimate["shape"],estfit.weibull$estimate["scale"])$statistic
        )  })
        zi_weibull_fit <- logspline(stats)
        values<-replicate(4000, {1 - plogspline(ks.test(w,"pweibull",shape= fit.weibull$estimate["shape"]
                                                    , scale = fit.weibull$estimate["scale"])$statistic,zi_weibull_fit)})
        mutations[i,13]<-paste(median(as.numeric(values)))
    }
    else {
        x<-x[!is.na(x)]
        w<-x[x!=0]
        fit.weibull<-fitdist(w,"weibull")
        shape = as.numeric(fit.weibull$estimate["shape"])
        scale = as.numeric(fit.weibull$estimate["scale"])
        df<-as.data.frame(x)
            colnames(df)<-"y"
        model_data <- list(
            n_obs = nrow(df), 
            y_shape = shape, 
            y_scale = scale, 
            zeros = theta,   
            y=df$y
        )
        model <- stan(file=model_name,data = model_data,cores = getOption("mc.cores", 1L))
        parameters<-extract(model)
        vector_zeros<-as.numeric(parameters[[1]])
        vector_shapes<-as.numeric(parameters[[2]])
        vector_scales<-as.numeric(parameters[[3]])
        for (v in 1:4000){
            r <- as.numeric(as.matrix(tibble(y=rweibull(n*(1-theta),
            fit.weibull$estimate["shape"],
            fit.weibull$estimate["scale"])) 
            %>%
            bind_rows(tibble(y=rep(0, n*theta)))))
            r_ext <- as.numeric(as.matrix(tibble(y=rweibull(n*(1-vector_zeros[v]),vector_shapes[v],vector_scales[v])) 
            %>%
            bind_rows(tibble(y=rep(0, n*vector_zeros[v])))))
        
            stats[v]<-as.numeric(ks.test(r,r_ext)$statistic)
  }
      
      zi_weibull_fit <- logspline(stats)
      for (s in 1:4000){
      spline_values[s]<-(1 - plogspline(ks.test(x,
      as.numeric(as.matrix(tibble(y=rweibull(n*(1-vector_zeros[s]),
      vector_shapes[s],vector_scales[s]))
      %>%
      bind_rows(tibble(y=rep(0, n*vector_zeros[s]))))))$statistic,zi_weibull_fit))
      }
      
      mutations[i,13]<-paste(median(spline_values))
  }
        if (as.numeric(mutations[i,13]) >0.01) {
        mutations[i,14]<-paste("Good_Weibull_fit")}
        else {        
        mutations[i,14]<-paste("wrong_model")
        }

}

colnames(mutations)[11:14]<-c("n_samples","theta","p.value")
write.table(mutations,args[4],sep="\t",row.names=FALSE)

