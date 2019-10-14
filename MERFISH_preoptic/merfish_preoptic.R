
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd("..")

library(Rtsne)
library(ggplot2)
library(gridExtra)


fname <- "MERFISH_preoptic/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv"
if (!file.exists(fname)) { download.file("https://datadryad.org/bitstream/handle/10255/dryad.192695/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv?sequence=1", fname) }


data5k <- read.csv("MERFISH_preoptic/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv", nrows = 5000) #read data for 5000 cells
heads <- colnames(lines)
mermeta <- data5k[,1:10]

merfish<-data5k[,10:170]
row.names(merfish) <- mermeta$Cell_ID

 
 merfish <- merfish[,!colnames(merfish)=="Fos"] #the "Fos" column included NA values
 genes1 <-colnames(merfish)[1:80]  #80 genes in view 1
 genes2 <-colnames(merfish)[81:160]  #80 genes in view 2

#sampel view 1
set.seed(0)
ix1<-sample(seq_len(nrow(merfish)), 1500)
mer.data1 <- merfish[ix1,genes1]
celltype.1 <-as.character( mermeta$Cell_class [ix1] )
#sampel view 2
ix2<-sample(seq_len(nrow(merfish)), 2000)
mer.data2 <- merfish[ix2,genes2]
celltype.2 <- as.character(mermeta$Cell_class [ix2] )
  
# unique(celltype.1)
#[1] "Astrocyte"     "Inhibitory"    "Pericytes"     "Endothelial 3" "Microglia"     "Ambiguous"     "OD Immature 1" "Endothelial 2" "Endothelial 1"
#[10] "Excitatory"    "OD Mature 2"   "OD Mature 1"   "OD Mature 4"   "OD Mature 3"   "OD Immature 2"
# Reduce and organize the cell types so that the Figure colour codes will be perceivable

celltype.1[substr(celltype.1,start=1,stop = 13)=="Endothelial 1"]= "Endothelial" #  
celltype.1[substr(celltype.1,start=1,stop = 13)=="Endothelial 2"]= "Endothelial" #  
celltype.1[substr(celltype.1,start=1,stop = 13)=="Endothelial 3"]= "Endothelial" 

celltype.2[substr(celltype.2,start=1,stop = 13)=="Endothelial 1"]= "Endothelial" 
celltype.2[substr(celltype.2,start=1,stop = 13)=="Endothelial 2"]= "Endothelial" 
celltype.2[substr(celltype.2,start=1,stop = 13)=="Endothelial 3"]= "Endothelial"   

celltype.1[substr(celltype.1,start=1,stop = 11)=="OD Immature"]= "OD Immature" 
celltype.1[substr(celltype.1,start=1,stop = 9)=="OD Mature"]= "OD Mature" 

celltype.2[substr(celltype.2,start=1,stop = 11)=="OD Immature"]= "OD Immature" 
celltype.2[substr(celltype.2,start=1,stop = 9)=="OD Mature"]= "OD Mature" 


## prepare 5 cell type memberships as shared features between the views
xsh1<-matrix(0, nrow= length(celltype.1),ncol=5)
xsh1[celltype.1=="Inhibitory",1]=1
xsh1[substr(celltype.1,start=1,stop = 6)=="Endoth" ,2]=1
xsh1[substr(celltype.1,start=1,stop = 2)=="OD" ,3]=1
xsh1[celltype.1=="Astrocyte",4]=1
xsh1[celltype.1=="Excitatory",5]=1

xsh2<-matrix(0, nrow= length(celltype.2),ncol=5)
xsh2[celltype.2=="Inhibitory",1]=1
xsh2[substr(celltype.2,start=1,stop = 6)=="Endoth" ,2]=1
xsh2[substr(celltype.2,start=1,stop = 2)=="OD" ,3]=1
xsh2[celltype.2=="Astrocyte",4]=1
xsh2[celltype.2=="Excitatory",5]=1

######### prepare the arguments for satsne_annealing
n1= nrow(mer.data1)
n2=nrow(mer.data2)

perplex_in1=min(floor(n1/5),250) # initial perplexity for view 1
perplex_in2=min(floor(n2/5),250) # initial perplexity for view 2
perplex_fin=40 #(minimum perplexity for either of the views)
perplex_steps=0.8 #reduce perplexities by 0.8* 

#### make independent t-SNE plots
set.seed(0)
ts1<- Rtsne(mer.data1)
ts2<- Rtsne(mer.data2)

df1<-data.frame(x=ts1$Y[,1],y=ts1$Y[,2],group=as.factor(celltype.1))                     
p1<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1) + xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#
df2<-data.frame(x=ts2$Y[,1],y=ts2$Y[,2],group=as.factor(celltype.2))                     
p2<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1) + xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#
grid.arrange(p1,p2 , nrow = 2 , ncol = 1)  
######## run satsne_annealing

source('code/d2p.R')
source('code/satsne_p.R')
source('code/satsne_annealing.R')

tsneX <- satsne_annealing (X10=mer.data1,X20=mer.data2,X1shared=xsh1,X2shared=xsh2,labels1=celltype.1, labels2=celltype.2,
                           perplex_in1=perplex_in1,perplex_in2=perplex_in2,perplex_fin=perplex_fin,perplex_steps=perplex_steps,
                           no.initiations=1, do.plot=TRUE )

#### make the aligned t-SNE plots
df1<-data.frame(x=tsneX$Y1[,1],y=tsneX$Y1[,2],group=as.factor(celltype.1))                     
p1<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1)+ xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#+ggtitle("Multi")#
df2<-data.frame(x=tsneX$Y2[,1],y=tsneX$Y2[,2],group=as.factor(celltype.2))
p2<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group))+coord_fixed(ratio = 1)+ xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#
grid.arrange(p1,p2 , nrow = 2 , ncol = 1)
