library("Rtsne")
library("umap")
library("ggplot2")
library("Rdimtools")
library("MASS")
library("kknn")
library("survival")
library("survminer")
library("NMI")

source("run_data_col.R")
source("survive_analyze.R")
source("specClust.R")
assignInNamespace('specClust',specClust,ns='kknn')
environment(specClust) <- asNamespace('kknn')

############################simulation data######################################################
##complete
true_label<-c()
for (i in seq(1,100)){
  if (i<=35){
    true_label<-c(true_label,1)
  }
  if (35<i && i<=65){
    true_label<-c(true_label,2)
  }
  if(i>65){
    true_label<-c(true_label,3)
  }
}

v1<-read.csv("Simulation_data/v1_2.csv",header = FALSE)
row.names(v1) <- seq(1,100)
v2<-read.csv("Simulation_data/v2_2.csv",header = FALSE)
row.names(v2) <- seq(1,100)
v3<-read.csv("Simulation_data/v3_2.csv",header = FALSE)
row.names(v3) <- seq(1,100)
data<-list(RNAseq=t(v1),miRNAseq=t(v2),protein=t(v3),clinical=NULL)
simul_result<-data_col(data,incomplete_data=F,
                       incomplete_sample_name,remain_view=1,dim1=3,dim2=3,pca_scale=F)
cl_s <- specClust(simul_result$l_space, 3, nn=20)
umap_r<-umap(simul_result$l_space,preserve.seed=FALSE)
plot(umap_r$layout,col=cl_s$cluster,pch=16,asp = 0.7,
     xlab = "UMAP_1",ylab = "UMAP_2",
     #main = "A UMAP visualization of the dataset"
)
print(NMI(X<-data.frame(seq(1,length(true_label)),true_label),
          Y<-data.frame(seq(1,length(cl_s$cluster)),cl_s$cluster)))

###incomplete
N=c(15,15)
C=c(1,3)
seed<-sample(seq(1,100000),1)

incomplete_sample_<-get_incomplete(colnames(data[[1]]),true_label,C,N,seed)
convert_label<-get_convert_labels(true_label,incomplete_sample_$index)
result_s<-data_col(data,incomplete_data=T,
                   incomplete_sample_$name,
                   remain_view=2,dim1=3,dim2=3,pca_scale=F,seed=seed)

cl_s <- specClust(result_s$l_space, 3, nn=8)
umap_r<-umap(result_s$l_space,preserve.seed=FALSE)

plot(umap_r$layout,col=cl_s$cluster)
print(seed)
print(NMI(X<-data.frame(seq(1,length(convert_label)),convert_label),
          Y<-data.frame(seq(1,length(cl_s$cluster)),cl_s$cluster)))



