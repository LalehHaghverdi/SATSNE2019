# Data simulation of the branching event not resolved by the feature set in view 2. 

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
setwd("..")

#dir.create("results", showWarning=FALSE)


library(Rtsne)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

load("Simulations/toydata.RData")
set.seed(0)

glob=matrix(rnorm(200*ncol(toydata),sd=30), ncol=ncol(toydata)) 
#seedi=1
#glob=matrix(rep(toydata[seedi,,drop=FALSE],each=500),nrow=500,ncol=ncol(toydata))+glob
glob=matrix(rep(c(0,0,100,100),each=200),nrow=200,ncol=ncol(toydata))+glob

timepoint.a<-as.numeric(cut(toydata[,1],breaks = 3)) 
timepoint.b<-sample(seq_len(3), nrow(glob),replace = TRUE)
timepoint= c(timepoint.a,timepoint.b)        ###from here  
  
toydata1=rbind(toydata,glob)
toydata2=rbind(toydata,glob)

set.seed(0)
ix1<-sample(seq_len(nrow(toydata1)),800)
ix2<-sample(seq_len(nrow(toydata2)),1200 ) 
Y1<-toydata1[ix1,c(1,2,3,4)]              
Y2<-toydata2[ix2,c(1,2)]        # Mask the branching features in view 2

Y1<-scale(Y1,center =TRUE)
Y2<-scale(Y2,center =TRUE)

timepoint1<- timepoint[ix1]
timepoint2<- timepoint[ix2] #as.numeric(cut(Y2[,1],breaks = 3))

ngenes1 <- 100
ncells1<-nrow(Y1)
set.seed(0)
proj1 <- matrix(rnorm(ngenes1*ncells1), nrow=ngenes1, ncol=ncol(Y1))
A1 <- Y1 %*% t(proj1) + matrix(rnorm(nrow(Y1)*(ngenes1) ),nrow=nrow(Y1),ncol=(ngenes1))

ngenes2 <- 200
ncells2<-nrow(Y2)
set.seed(10)
proj2 <- matrix(rnorm(ngenes2*ncells2), nrow=ngenes2, ncol=ncol(Y2))
A2 <- Y2 %*% t(proj2) + matrix(rnorm(nrow(Y2)*(ngenes2) ),nrow=nrow(Y2),ncol=(ngenes2))

ts1<- Rtsne(A1[,2:ncol(A1)],perplexity = 20)
ts2<- Rtsne(A2[,2:ncol(A2)],perplexity = 20)

A1<-cbind(timepoint1,A1)
A2<-cbind(timepoint2,A2)

################## blue branches

rbPal <- colorRampPalette(c('blue','yellow','red'))
datCol1 <- rbPal(3)[as.numeric(cut(Y1[,1],breaks = 3))]
datCol2 <- rbPal(3)[as.numeric(cut(Y2[,1],breaks = 3))]

df1<-data.frame(x=ts1$Y[,1],y=ts1$Y[,2],group=as.factor(timepoint1))                   
p10<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1)+ xlab("t-SNE1") + ylab("t-SNE2") +theme(legend.position="none")#+ggtitle("Multi")#+ scale_color_hue(l=40, c=35)
df2<-data.frame(x=ts2$Y[,1],y=ts2$Y[,2],group=as.factor(timepoint2)) 
p20<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group))+coord_fixed(ratio = 1) +xlab("t-SNE1") + ylab("t-SNE2") + theme(legend.position="none")#+ scale_color_hue(l=40, c=35)
grid.arrange(p10,p20 , nrow = 2 , ncol = 1)
# save as 700*700 png
#stop()

########## prepare the two views and params for satsne
X10<-A1[,2:100]
X20<-A2[,2:200]
X1shared=A1[,1,drop=FALSE]
X2shared=A2[,1,drop=FALSE]

nk1=250 #nns in data1
nk2=250 #nns in data2
L1=10#8
L2=10#8

n1<-nrow(X10)
n2<-nrow(X20)
perplex_in1=min(floor(n1/5),250) #130
perplex_in2=min(floor(n2/5),250)
perplex_fin=30 
perplex_steps=0.8
no.initiations=3

#### run SATSNE
source('code/d2p.R')
source('code/satsne_p.R')
source('code/satsne_annealing.R')

tsneX <- satsne_annealing (X10,X20,X1shared,X2shared,labels1=timepoint1,labels2=timepoint2,
                            nk1=nk1, nk2=nk2,L1=L1, L2=L2, no_dims=2,
                            perplex_in1=perplex_in1,perplex_in2=perplex_in2,perplex_fin=perplex_fin, perplex_steps=perplex_steps,
                            no.initiations=no.initiations,Y1_init=NULL,Y2_init=NULL, max_iter = 200, do.plot=TRUE )

df1<-data.frame(x=tsneX$Y1[,1],y=tsneX$Y1[,2],group=as.factor(timepoint1))                     
p1<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1)+ xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#+ggtitle("Multi")#
df2<-data.frame(x=tsneX$Y2[,1],y=tsneX$Y2[,2],group=as.factor(timepoint2))
p2<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group))+coord_fixed(ratio = 1)+ xlab("t-SNE1") + ylab("t-SNE2")  + theme(legend.position="none")#
grid.arrange(p1,p2 , nrow = 2 , ncol = 1)

