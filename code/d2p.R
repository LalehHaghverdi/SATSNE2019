
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


d2p <- function (D, u=30, tol=1e-4) {
#u is the perplexity
# %d2p Identifies appropriate sigma's to get kk NNs up to some tolerance

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

