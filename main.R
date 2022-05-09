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

################################### Completed data clustering result evaluation ##################################################

#import the dataset 
data<-brainlower_re

#Compute the latent representation by data_col under the completed data situation
result1<-data_col(data,incomplete_data=F,
                 incomplete_sample_name,remain_view=1,dim1=30,dim2=5,pca_scale=F)

#First clustering
cl1 <- specClust(result1$l_space, 5, nn=3)
umap<-umap(result1$l_space)
plot(umap$layout,col=(cl1$cluster),pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2",
     main = "A UMAP visualization of the dataset")
legend("topright",title = "Species",inset = 0.01,
       legend = unique(cl1$cluster),pch=16,
       col = unique(cl1$cluster))
#survive analyze
survive_result(data[[4]],cl1$cluster)


################################# Incomplete data clustering result evaluation#################################################

#First step (get proper cluster labels) 
result2<-data_col(data,incomplete_data=F,
                  incomplete_sample_name,remain_view=1,dim1=30,dim2=5,pca_scale=F)
cl2 <- specClust(result2$l_space,5, nn=3)
umap2<-umap(result2$l_space)
plot(umap2$layout,col=(cl2$cluster),pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2"
     )
legend("topright",title = "Cluster",inset = 0.01,
       legend = sort(unique(cl2$cluster)),pch=16,
       col = sort(unique(cl2$cluster)))



#Second step (convert the incomplete data to same latent space)


N=c(15)#the number of the selected incomplete samples
C=c(2)#the ID of the cluster containing the selected incomplete samples
#seed<-sample(seq(1,100000),1)
seed<-24652

#select the random sample and make them be incomplete
incomplete_sample_<-get_incomplete(colnames(data[[1]]),cl2$cluster,C,N,seed)

#Compute the latent representation by data_col under the incomplete data situation
result3<-data_col(data,incomplete_data=T,
                  incomplete_sample_$name,
                 remain_view=2,dim1=30,dim2=5,pca_scale=F,seed=seed)

convert_label<-get_convert_labels(cl2$cluster,incomplete_sample_$index)

cl3 <- specClust(result3$l_space, 4, nn=8)
umap3<-umap(result3$l_space)

plot(umap3$layout,col=(convert_label),pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2"
     )


legend(-6,6.5,title = "Clusters of complete data",
       legend = seq(1,length(unique(convert_label))),pch=16,
       col = seq(1,length(unique(convert_label))),horiz=T,xpd=T,box.lwd = 0)

legend(2,6.5,title = "Incomplete data",inset = 0.01,
       legend = c(5), pch=c(15),
       col = c("purple"),box.lwd = 0,xpd=T,horiz=T)

#here 425 is the total sample number 
points(x=as.matrix(umap3$layout[(425-N[1]+1):425,])[,1],
       y=as.matrix(umap3$layout[(425-N[1]+1):425,])[,2],col='purple',pch=15,cex=1.5)
print(seed)

#evaluate the clustering performance under the incomplete data situation
print(NMI(X<-data.frame(seq(1,length(convert_label)),convert_label),
    Y<-data.frame(seq(1,length(convert_label)),cl3$cluster)))









############original points distribution#################
plot(umap2$layout,col=(cl2$cluster),pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2",
     #main = "A UMAP visualization of the dataset"
     )
legend(-6.3,8,title = "Clusters",inset = 0.01,
       legend = c(1,2,3,4),pch=16,
       col = c(1,2,3,4),horiz=T,xpd=T,box.lwd = 0)
legend(1,8,title = "Selected samples",inset = 0.01,
       legend = c(2), pch=c(15),
       col = c("purple"),box.lwd = 0,xpd=T,horiz=T)


points(x=as.matrix(umap2$layout[incomplete_sample_$index[1:N[1]],])[,1],
       y=as.matrix(umap2$layout[incomplete_sample_$index[1:N[1]],])[,2],col='purple',pch=15,cex=1.5)

