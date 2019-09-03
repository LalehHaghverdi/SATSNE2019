satsne_annealing <-function (X10,X20,X1shared=NULL,X2shared=NULL,nk1=200, nk2=200,L1=10, L2=10, no_dims=2,
                              perplex_in1=250,perplex_in2=250,perplex_fin=30, perplex_steps=0.8,no.initiations=1,
                              Y1_init=NULL,Y2_init=NULL,max_iter = 200 , do.plot=TRUE ) {
###  
n1<-nrow(X10)
n2<-nrow(X20)

#nk1=250 #nns in data1 #150
#nk2=250 #nns in data2
#L1=10#8
#L2=10#8

 
#no.initiations=1
#perplex_in1=min(floor(n1/5),250) #130
#perplex_in2=min(floor(n2/5),250)
#perplex_fin=40 
#perplex_steps=20

max_perplex_in = max(perplex_in1,perplex_in2)


perplex_rounds=floor((max_perplex_in - perplex_fin)/perplex_steps) #4

source('code/satsne_p.R')
source('code/d2p.R')
library(FNN)
library(ggplot2)
library(gridExtra)

firstview="v1"  #"RNA-seq"
secondview="v2"  #ATAC-seq"

#X1shared= NULL
#X2shared= NULL
X1sh<-X1shared
X2sh<-X2shared


Y1_init=NULL#cbind(X1shared,X1shared)#NULL
Y2_init=NULL#cbind(X2shared,X2shared)#NULL

#############
if  ( (ncol(X10) > 30) & ncol(X20) > 30 ){
X1<-prcomp(X10)$x[,1:20]
X2<-prcomp(X20)$x[,1:20]
} else {
  X1<- X10
  X2<-X20
}

perplex1=perplex_in1
perplex2=perplex_in2

P1<-d2p(as.matrix(dist(X1)),u=perplex1)$P
P2<-d2p(as.matrix(dist(X2)),u=perplex2)$P

set.seed(1) 
#Take rather long (>=500)iterations in the begining with the largest perplexity 
#tsneX<-satsne_p(P1,P2,X1shared=X1shared,X2shared=X2shared,Y1_init=Y1_init,Y2_init=Y2_init,
#                  nk1=nk1, nk2=nk2,L1=L1, L2=L2,no_dims =no_dims,max_iter = 500) #700

par(mfrow=c(1,2))
cost<-vector()

#i=0
while ( ((perplex1) > perplex_fin) & ((perplex2) > perplex_fin)  ) {
#i=i+1
  
   #try several initialisations at i=1
    mincost=10^8 # a very large number
   if ( (is.null(Y1_init) | is.null(Y2_init)) ) { # A) first (largest) perplexities
    for (j in 1:no.initiations)
    tsneXn<-satsne_p(P1,P2,X1shared=X1shared,X2shared=X2shared,Y1_init=Y1_init,Y2_init=Y2_init,
                  nk1=nk1, nk2=nk2,L1=L1, L2=L2,no_dims =no_dims,max_iter = 700)   #for the first perplexity use at least 700 iterations
    if ( min(tsneXn$Costs) < mincost ) {
    tsneX<-tsneXn
    mincost= min(tsneX$Costs) }else {rm(tsneXn)}
  } else { # B) after the first perplexities
      
      # perplex1=perplex1-perplex_steps
      # perplex2=perplex2-perplex_steps 
      perplex1=floor(perplex1*perplex_steps)
      perplex2=floor(perplex2*perplex_steps) 
    
      P2<-d2p(as.matrix(dist(X2)),u=perplex1)$P
      P1<-d2p(as.matrix(dist(X1)),u=perplex2)$P
     
      if ( (perplex1<=perplex_fin) | perplex2<=perplex_fin) { # use at least 700 iterations for the final perplexity.
        max_iter <- 1000
      }
      tsneX<-satsne_p(P1,P2,X1shared=X1sh,X2shared=X2sh,Y1_init=Y1_init,Y2_init=Y2_init,
                      nk1=nk1, nk2=nk2,L1=L1, L2=L2,no_dims =no_dims,max_iter = max_iter) #700
    }


Y1_init=tsneX$Y1
Y2_init=tsneX$Y2

cost<-c(cost,tsneX$Costs)
print(paste0("perplexity=",c(perplex1, perplex2)) )

if (do.plot==TRUE) {
library(ggplot2)
library(gridExtra)
  
if (exists("datCol1") & exists("datCol2") ) { 
par(mfrow=c(2,1))
plot(tsneX$Y1,col=datCol1,pch=20, main=perplex1,xlab="",ylab="")
plot(tsneX$Y2,col=datCol2,pch=20,main=perplex2,xlab="",ylab="")
 }else if (exists("celltype.1") & exists("celltype.1")) {
   
   df1<-data.frame(x=tsneX$Y1[,1],y=tsneX$Y1[,2],group=as.factor(celltype.1))                     
   p1<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1) +ggtitle(perplex1) + theme(legend.position="none")#+ scale_color_hue(l=40, c=35)
   df2<-data.frame(x=tsneX$Y2[,1],y=tsneX$Y2[,2],group=as.factor(celltype.2))
   p2<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group))+coord_fixed(ratio = 1)+ggtitle(perplex2) + theme(legend.position="none")#+ scale_color_hue(l=40, c=35)
   grid.arrange(p1,p2 , nrow = 2 , ncol = 1)  
 }
}

}
return(list("Y1"=tsneX$Y1, "Y2"=tsneX$Y2, "Costs"=cost))
}
#############
