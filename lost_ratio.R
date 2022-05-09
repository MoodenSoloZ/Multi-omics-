library("Rtsne")
library("umap")
library("ggplot2")
library("Rdimtools")
library("MASS")
library("kknn")
library("survival")
library("survminer")
library("NMI")


library(Spectrum)
library(SNFtool)

source("run_data_col.R")
source("survive_analyze.R")
source("specClust.R")
assignInNamespace('specClust',specClust,ns='kknn')
environment(specClust) <- asNamespace('kknn')

data<-simu_data
cluster_number<-3
final_result1<-c()
time_final<-c()
r1<-data_col(data,incomplete_data=F,incomplete_sample_name,remain_view=1,dim1=20,dim2=10,pca_scale=F)
cl1 <- specClust(r1$l_space,cluster_number, nn=3)
used_label<-true_label
for (ratio in seq(0.02,0.4,by=0.02)){
  print(paste("This time ratio is",ratio))

  N<-c()
  for (i in seq(1,cluster_number)){
    N<-c(N,round(summary(as.factor(used_label))[i]*ratio))
  }
  C=seq(1,cluster_number)
  temp_result<-c()
  time_result<-c()
  for (iter in seq(1,10)){
    seed<-sample(seq(1,100000),1)
    t1<-proc.time()[3]
    incomplete_sample_<-get_incomplete(colnames(data[[1]]),used_label,C,N,seed)
    r2<-data_col(data,incomplete_data=T,
                 incomplete_sample_$name,
                 remain_view=1,dim1=3,dim2=3,pca_scale=F,seed=seed)
    
    convert_label<-get_convert_labels(used_label,incomplete_sample_$index)
    cl2 <- specClust(r2$l_space, cluster_number, nn=20)
    t2<-proc.time()[3]
    temp_result<-c(temp_result,NMI(X<-data.frame(seq(1,length(convert_label)),convert_label),
                                   Y<-data.frame(seq(1,length(convert_label)),cl2$cluster))$value)
    time_result<-c(time_result,t2-t1)
  }
  final_result1<-c(final_result1,mean(temp_result))
  time_final<-c(time_final,mean(time_result))
}
print(final_result1)

