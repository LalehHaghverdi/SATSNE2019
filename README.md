# SATSNE2019
R code for Soft Alignment of t-SNEs


If the shared features are gene expression values, we recommend cosine (L2) normalisation before passing it as an argument of satins_p (alignment with fixed perplexities) or satsne_annealing (alignment with annealing) functions.

Example: 

n1=500 #number of cells in view 1
n2=300 #number of cells in view 2
X10=matrix(rnorm(n1*5,sd=c(30, 20,20,10, 30)), ncol=5 ) 
X20=matrix(rnorm(n2*3,sd=c(30,10,5)), ncol=3 ) 
X1shared <- X10[,1, drop=FALSE]
X2shared <- X20[,1, drop=FALSE]
X1shared <- cosine.norm(X1shared)
X2shared <- cosine.norm(X2shared)

library(RColorBrewer)
rbPal <- colorRampPalette(c('blue','yellow','red'))
datCol1 <- rbPal(10)[as.numeric(cut(X1shared[,1],breaks = 10))]
datCol2 <- rbPal(10)[as.numeric(cut(X2shared[,1],breaks = 10))]

tsneX <- satsne_annealing (X10,X20,X1shared,X2shared,perplex_in1=60,perplex_in2=floor(n2/n1*60),perplex_fin=30)
plot(tsneX$Y1, col=datCol1)
plot(tsneX$Y2, col=datCol2)