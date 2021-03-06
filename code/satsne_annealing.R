satsne_annealing <-function (X1view, X2view, X1shared=NULL, X2shared=NULL, labels1=NULL, labels2=NULL,
                             is_distance= FALSE, nk1=200, nk2=200,L=10, no_dims=2,
                             perplex_in1=250, perplex_in2=250, perplex_fin=30, perplex_steps=0.8, no.initiations=1,
                             Y1_init=NULL, Y2_init=NULL, max_iter = 300, do.plot=TRUE ) {
###  
# X1view data view 1 [cells *features]
# X2view data view 2 [cells *features]
#X1shared shared feature values in view1 
#X2shared shared feature values in view2 

#nk1: nns in data1 #150
#nk2: nns in data2
#L: box size
 
#no.initiations=1  no. of random embedding initiations
#perplex_in1=min(floor(n1/5),250) # initial perplexity for view 1
#perplex_in2=min(floor(n2/5),250) # initial perplexity for view 2
#perplex_fin:  approx. final (minimum) perplexity for either of the views
#perplex_steps=0.8 reduce perplexities by 0.8* 

  start_time <- Sys.time()
  
  library(FNN)
  library(ggplot2)
  library(gridExtra)
  library(here)  
  
  source(here("projects/SATSNE2019/code/support_functions.R"))
  
  #uses labels1 and labels2 for colouring of the plots (if defined and do.plot=TRUE)  
  
  n1<-nrow(X1view)
  n2<-nrow(X2view)
  X1shared<- as.matrix(X1shared,drop=FALSE)
  X2shared<- as.matrix(X2shared,drop=FALSE)
  
  #############
  perplex1=perplex_in1
  perplex2=perplex_in2
  
  if (is_distance) {
    D1 <- as.matrix(X1view)
    D2 <- as.matrix(X2view)
  }else if  ( (ncol(X1view) > 30) & ncol(X2view) > 30 ){
    #X1<-prcomp(X1view)$x[,1:20] #use 20 first pcs
    #X2<-prcomp(X2view)$x[,1:20]
    X1<-X1view %*% svd(X1view, nu=20, nv=20)$v #use 20 first svd components
    X2<-X2view %*% svd(X2view, nu=20, nv=20)$v
    D1 <- as.matrix(dist(X1))
    D2 <- as.matrix(dist(X2))
  }else {
    X1<- X1view
    X2<- X2view
    D1 <- as.matrix(dist(X1))
    D2 <- as.matrix(dist(X2))
  }
  
  P1<-d2p(D1,u=perplex1)$P
  P2<-d2p(D2,u=perplex2)$P
  
  
  set.seed(1)
  par(mfrow=c(1,2))
  cost<-vector()
  #cost1 = rep(NA, max_iter)
  
  
  while ( ((perplex1) > perplex_fin) & ((perplex2) > perplex_fin)  ) {
    
    mincost=10^8 # start with a very large number
    if ( (is.null(Y1_init) | is.null(Y2_init)) ) { # A) first (largest) perplexities
      for (j in 1:no.initiations) {
        tsneXn<-satsne_p(P1,P2,X1shared=X1shared,X2shared=X2shared,Y1_init=Y1_init,Y2_init=Y2_init,
                         nk1=nk1, nk2=nk2,L=L,no_dims =no_dims,max_iter = 300,initiation.round=TRUE)   #for the first perplexity use at least 200 iterations
        if ( min(tsneXn$Costs) < mincost ) {
          tsneX<-tsneXn
          mincost= min(tsneX$Costs) }else {rm(tsneXn)}
      }
      print(mincost) ### to test, remove later
    } else { # B) next perplexities after the first one
      
      perplex1=floor(perplex1*perplex_steps)
      perplex2=floor(perplex2*perplex_steps) 
      
      P1<-d2p(D1,u=perplex1)$P
      P2<-d2p(D2,u=perplex2)$P
      
      if ( (perplex1<=perplex_fin) | perplex2<=perplex_fin) { # use at least 200 iterations for the final perplexity.
        max_iter <- 300 # 500
        no_dims=2
        Y1_init=Y1_init[,1:no_dims]
        Y2_init=Y2_init[,1:no_dims]
      }
      tsneX<-satsne_p(P1,P2,X1shared=X1shared,X2shared=X2shared,Y1_init=Y1_init,Y2_init=Y2_init,
                      nk1=nk1, nk2=nk2,L=L,no_dims =no_dims,max_iter = max_iter,initiation.round=FALSE) #700
    }
    
    Y1_init=tsneX$Y1
    Y2_init=tsneX$Y2
    
    cost<-c(cost,tsneX$Costs)
    print(paste0("perplexity1=",perplex1," perplexity2=", perplex2)) 
    
    if (do.plot==TRUE) {
      # plottings for each perplexity round    
      if (exists("labels1") & exists("labels2")) {
        
        df1<-data.frame(x=tsneX$Y1[,1],y=tsneX$Y1[,2],group=as.factor(labels1))                     
        p1<-ggplot(df1)+geom_point(aes(x,y,color=group,fill=group)) +coord_fixed(ratio = 1) +ggtitle(perplex1) + theme(legend.position="none")#+ scale_color_hue(l=40, c=35)
        df2<-data.frame(x=tsneX$Y2[,1],y=tsneX$Y2[,2],group=as.factor(labels2))
        p2<-ggplot(df2)+geom_point(aes(x,y,color=group,fill=group))+coord_fixed(ratio = 1)+ggtitle(perplex2) + theme(legend.position="none")#+ scale_color_hue(l=40, c=35)
        grid.arrange(p1,p2 , nrow = 1 , ncol = 2)  
      }
    }
    
  }
  
  end_time <- Sys.time()
  minutes= end_time-start_time
  
  return(list("Y1"=tsneX$Y1, "Y2"=tsneX$Y2, "Costs"=cost, "Runtime"=minutes))
}
#############
