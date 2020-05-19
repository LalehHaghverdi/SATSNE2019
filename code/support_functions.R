satsne_p <-function (P1,P2,X1shared=NULL,X2shared=NULL,nk1=20, nk2=20,L=10, no_dims=3,Y1_init=NULL,Y2_init=NULL,max_iter = 200,initiation.round=FALSE) {

# Adapted from tsne-p Matlab code by Laurens van der Maaten, 2010.
# Taking precalculated affinity matrices P1 and P2 (with fixed perplexity), this function iterates between coupled and uncoupled t-sne modes
# for a soft alignment of the embeddings. 
# outputs is a list of: the partly aligned embeddings "Y1" and "Y2", sum of the two views' t-sne costs ("Costs") and the execution time in minutes ("Runtime").   

#nk1=10 #nns of data2 in data1 
#nk2=10  #nns of data1 in data2 
#L #tsne box size

start_time <- Sys.time()
  
library(FNN)

# Initialize some variables
n1 = nrow(P1)                                     # number of instances (cells)
n2 = nrow(P2)                                     # number of instances (cells)

if (initiation.round==TRUE) {bindingforce=1}
else{bindingforce=0.05} #0.5 

momentum = 0.5;                                     #% initial momentum
#final_momentum = 0.8;                               #% value to which momentum is changed
mom_switch_iter = 250;                              #% iteration at which momentum is changed
stop_lying_iter = floor(max_iter/4);                              #% iteration at which lying about P-values is stopped
#max_iter = 1000;                                    #% maximum number of iterations
epsilon = 500;                                      #% initial learning rate
min_gain = .01;                                     #% minimum gain for delta-bar-delta


coupled_period=30  #(1st)
uncoupled_period=20 #(2nd)

one_whole_period =coupled_period + uncoupled_period 
# set max_iter such than the last iterration lands in an early uncoupled period stage 
max_iter= max_iter  - (max_iter %% one_whole_period ) + coupled_period + 10 *(!initiation.round) # 10


# Make sure P-vals are set properly
diag(P1) = 0                                 #% set diagonal to zero
P1 = 0.5 * (P1 + t(P1));                     #% symmetrize P-values
realmin=1e-18
P1 = (P1 / sum(P1))
P1[P1< realmin]= realmin;                   #% make sure P-values sum to one #maybe I don't need this

diag(P2) = 0                                 #% set diagonal to zero
P2 = 0.5 * (P2 + t(P2));                    #% symmetrize P-values
realmin=1e-18
P2 = (P2 / sum(P2))
P2[P2< realmin]= realmin;

const1 = sum(P1 * log(P1) );                     #% constant in KL divergence
const2 = sum(P2 * log(P2) );


if ( is.null(Y1_init) ) {
P1 = P1 * 4                                      # lie about the P-vals to find better local minima
#set.seed(125)
ydata1 <- .0001 * matrix(rnorm(n1*no_dims), nrow = n1, ncol = no_dims)
}else {    
ydata1 <- Y1_init
}
if ( is.null(Y2_init) ) {
  P2 = P2 * 4                                     # lie about the P-vals to find better local minima
#  set.seed(120)
  ydata2 <- .0001 * matrix(rnorm(n2*no_dims), nrow = n2, ncol = no_dims)
}else {     
ydata2 <- Y2_init
}

y_incs1  = matrix(0,nrow(ydata1),ncol(ydata1) );
gains1 = matrix(1,nrow(ydata1),ncol(ydata1) );

y_incs2  = matrix(0,nrow(ydata2),ncol(ydata2) );
gains2 = matrix(1,nrow(ydata2),ncol(ydata2) );


cost1=vector()
cost2=vector()
######Finding MNNs on Xshareds and organize
if (!is.null(X1shared) ) { #& (iter > (max_iter/2) ) ) {
  ydata11= X1shared 
  ydata22= X2shared 
} else {
  ydata11=ydata1 
  ydata22=ydata2 
  }

mat <- find.mutual.nn(ydata22,ydata11,nk1,nk2)$mnns  # mnns in Xshared dimensions

Z1 <- sapply(seq_len(n2),function(i) mat[(mat[,1]==i),2])  #nns of data2 in data1
Z2 <- sapply(seq_len(n1),function(i) mat[(mat[,2]==i),1])  #nns of data1 in data2


mat21<-do.call(rbind,lapply(Z1,
                            function(x)
                              if(length(x)==0) rep(0,nk1)
                            else if(length(x)<nk1) c(x,rep(0,nk1-length(x)))
                            else x)) # nns of data 2 in data1

mat12<-do.call(rbind,lapply(Z2,
                            function(x)
                              if(length(x)==0) rep(0,nk2)
                            else if(length(x)<nk2) c(x,rep(0,nk2-length(x)))
                            else x))

mat12[mat12==0]<-NaN
mat21[mat21==0]<-NaN  

#########begin the uncoupled/coupled iterations
for (iter in 1:max_iter) {

if ( ((iter %% one_whole_period ) < coupled_period) ) { #start with a coupled period 
 # adopt binding forces for the coupling 
 bfc1=rep(0,n1)
 bfc2=rep(0,n2)
 bfc1=rep(bindingforce,n1)   #to pull tsne2 towards tsne1
 bfc2=rep(bindingforce,n2)   #to pull tsne1 towards tsne2   
 #momentum=0.5
 #epsilon=500

 } else{ #in uncoupled periods 
  # turns off couplings  
 bfc1=rep(0,n1) 
 bfc2=rep(0,n2)
 #momentum=0.5
 #epsilon=500
 }

# Compute joint probability that point i and j are neighbors
euc.dist_y1  <- dist(ydata1, method = "euclidean")       # pairwise distances in low-dim
D_y1 <- as.matrix(euc.dist_y1)
N1 <- 1/(1 + (D_y1)^2)
diag(N1)<-0
Q1 <- N1 / sum(N1)
Q1[Q1< realmin]= realmin;

euc.dist_y2  <- dist(ydata2, method = "euclidean")       # pairwise distances in low-dim
D_y2 <- as.matrix(euc.dist_y2)
N2 <- 1/(1 + (D_y2)^2)
diag(N2)<-0
Q2 <- N2 / sum(N2)
Q2[Q2< realmin]= realmin;

# Compute the gradients (faster implementation)
R1 = (P1 - Q1) * N1;
y_grads1 = 4 * (diag(rowSums(R1)) - R1) %*% ydata1;  

R2 = (P2 - Q2) * N2;
y_grads2 = 4 * (diag(rowSums(R2)) - R2) %*% ydata2;  

######
# Update the solution
gains1 = (gains1 + .2) * (sign(y_grads1) != sign(y_incs1))          #% note that the y_grads are actually -y_grads
+ (gains1 * .8) * (sign(y_grads1) == sign(y_incs1));
gains1[gains1 < min_gain] = min_gain;

gains2 = (gains2 + .2) * (sign(y_grads2) != sign(y_incs2))          #% note that the y_grads are actually -y_grads
+ (gains2 * .8) * (sign(y_grads2) == sign(y_incs2));
gains2[gains2 < min_gain] = min_gain;

####
if (iter < (max_iter /2) && is.null(Y1_init) && is.null(Y2_init) ){     # A) in early iterrations couple randomly with one of the sharedX MNNs 
  match21= lapply(Z2, function(x)
                        if(length(x)==0) NaN
                        else sample(seq_len(length(x)),1) ) 
  
  match12= lapply(Z1, function(x)
                        if(length(x)==0) NaN
                        else sample(seq_len(length(x)),1) )
  match21 <-unlist(match21)
  match12 <-unlist(match12) 
  
} else{           # B) in late iterrations couple with the MNNs which is the closeset in embedding space 
  
yd1 <- ydata1[mat21,] 
yd2 <- do.call(rbind, replicate(nk1, ydata2, simplify=FALSE))
tsdist12 <- matrix( rowSums((yd2-yd1)^2) , n2, nk1)
##
yd2 <- ydata2[mat12,] 
yd1 <- do.call(rbind, replicate(nk2, ydata1, simplify=FALSE))
tsdist21 <- matrix( rowSums((yd1-yd2)^2) , n1, nk2)

  match21 <- apply(tsdist21,1,which.min)
  match12 <- apply(tsdist12,1,which.min)

match21 <- sapply( match21,
                  function(x) if(length(x)==0) NaN
                               else x) 
match12 <- sapply( match12,
                   function(x) if(length(x)==0) NaN
                   else x)
}
indsmat21 <- cbind(seq_len(n1),match21) #same dims as matching 21
matching21 <- mat12[indsmat21] 

indsmat12 <- cbind(seq_len(n2),match12)
matching12 <- mat21[indsmat12]

####
mom.dirct1 <- (-ydata1+ydata2[matching21,] )#1:no_dims])
mom.dirct1[is.na(mom.dirct1)]=0
y_incs1 = bfc1 * mom.dirct1 + momentum * y_incs1 - epsilon * (gains1 * y_grads1);
  ydata1 = ydata1 + y_incs1;

mom.dirct2 <- (-ydata2+ydata1[matching12,]) #1:no_dims])
mom.dirct2[is.na(mom.dirct2)]=0
  y_incs2 = bfc2 * mom.dirct2 + momentum * y_incs2 - epsilon * (gains2 * y_grads2);
  ydata2 = ydata2 + y_incs2;

# scale the embeddings in a fixed-size box of length L
ydata1<-L*apply(ydata1, 2, function(y) (y-min(y))/(max(y)-min(y)) ) 
ydata2<-L*apply(ydata2, 2, function(y) (y-min(y))/(max(y)-min(y)) )

# Update the momentum if necessary
#if (iter == mom_switch_iter){
#momentum = final_momentum;
#}
if (iter == stop_lying_iter && is.null(Y1_init) ) {
P1 = P1 / 4;
P2 = P2 / 4;
}

cost1[iter] = const1 - sum(P1 * log(Q1) );
cost2[iter] = const2 - sum(P2 * log(Q2) );

}
end_time <- Sys.time()
minutes= end_time-start_time
return(list("Y1"=ydata1,"Y2"=ydata2,"Costs"=cost1+cost2, "Runtime"=minutes))
}

#########
find.mutual.nn <- function(data1, data2, k1, k2)
  # Finds mutal neighbors between data1 and data2, using k1 and k2 nns respectively
{
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n.total <- n1 + n2
  
  #print(dim(data1))
  
  W21 <- FNN::get.knnx(data2, query=data1, k=k1)
  W12 <- FNN::get.knnx(data1, query=data2, k=k2)
  
  inds12<-cbind( as.vector(rep(seq_len(n1), each=k1)),  as.vector(t(W21$nn.index)) )
  inds21<-cbind(as.vector(t(W12$nn.index)) , as.vector(rep(seq_len(n2), each=k2)) )
  
  A<-rbind(inds12,inds21)
  keeps=duplicated(A)  ##duplicated rows of A are mutaul nns
  A<-A[keeps,]
  
  A1<-A[,1]
  A2<-A[,2]
  # Report cells that are MNNs.
  list(mnns=A)
}




# % Function that computes the Gaussian kernel values given a vector of
# % squared Euclidean distances, and the precision of the Gaussian kernel.
# % The function also computes the perplexity of the distribution.

Hbeta <- function(Di, betai) {
  #betai scalar (=2*sigma_i^2)
  #Di vector =D[i,] dist matrix
  Pi = exp(-Di * betai);
  #Pi[i]=0 # diagonal P elements
  
  sumPi = sum(Pi) -1; #minus 1 for diagonal P elements
  realmin=1e-18
  if (sumPi <realmin) {sumPi=realmin}
  
  Hi = log(sumPi) + betai * sum(Di * Pi) / sumPi;
  
  Pi = Pi / sumPi;
  (list("Hi"=Hi,"Pi"=Pi)) #Hi scalar, Pi vector
}

###########
d2p <- function (D, u=30, tol=1e-4) {
  #u is the perplexity
  # %d2p Identifies appropriate sigma's to get kk NNs up to some tolerance,
  # and turns a distance matrix d to a transition probability matrix P
  
  # % Identifies the required precision (= 1 / variance^2) to obtain a Gaussian
  # % kernel with a certain uncertainty for every datapoint. The desired
  # % uncertainty can be specified through the perplexity u (default = 30). The
  # % desired perplexity is obtained up to some tolerance that can be specified
  # % by tol (default = 1e-4).
  # % The function returns the final Gaussian kernel in P, as well as the
  # % employed precisions per instance in beta.
  # %
  # %
  # % (C) Laurens van der Maaten, 2008
  # % Maastricht University
  
  # Initialize some variables
  n = nrow(D);                     # number of instances
  P = matrix(0,n, n);                    #% empty probability matrix
  #beta = matrix(1,n, 1);                  #% empty precision vector
  beta = rep(1,n);
  logU = log(u);                      #% log of perplexity (= entropy)
  
  # Run over all datapoints
  for (i in 1:n) {
    
    # Set minimum and maximum values for precision
    betamin = -Inf;
    betamax = Inf;
    
    # Compute the Gaussian kernel and entropy for the current precision
    Hbetai= Hbeta(D[i, ], beta[i]);
    Hi =Hbetai$Hi
    thisP =Hbetai$Pi
    #thisP[i] =0 # diagonal P elements
    # Evaluate whether the perplexity is within tolerance
    Hdiff = Hi - logU;
    tries = 0;
    
    while ( (abs(Hdiff) > tol) && (tries < 50) ) {
      # If not, increase or decrease precision
      if (Hdiff > 0){
        betamin = beta[i]
        if (is.infinite(betamax) ) {
          beta[i] = beta[i] * 2;
        }else{
          beta[i] = (beta[i] + betamax) / 2;
        }
      }else{
        betamax = beta[i];
        if (is.infinite(betamin)) {
          beta[i] = beta[i] / 2;
        }else{
          beta[i] = (beta[i] + betamin) / 2;
        }
      }
      # Recompute the values
      Hbetai= Hbeta(D[i, ], beta[i]);
      Hi =Hbetai$Hi
      thisP =Hbetai$Pi
      #thisP[i] =0 # diagonal P elements
      Hdiff = Hi - logU;
      tries = tries + 1;
    } # while end
    
    #Set the final row of P
    P[i, ] = thisP;
  } #enf for i
  
  return(list("P"=P,"beta"=beta))
}

cosine.norm <- function(X)
  # Computes the cosine norm, with some protection from zero-length norms.
  # Cell rows and Gene columns
{
  cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
  X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
}

