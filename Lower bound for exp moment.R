### Search for a lower bound for EexpWn

N = 10000; delta = 1; c = 0.99999;
n=1:N;
k = numeric(N); for (j in n){ m=1:j; k[j] = sum(m*exp(m) <= j); }
h = delta*sqrt(k)*qnorm(c)-delta^2*k/2;
LB = exp(h)*( 1-(pnorm( (h+delta^2*k/2)/(delta*sqrt(k))))^(n/k) );

plot(n,LB,type="l")


# Recursive computation of EU = EexpW for 1,...,n
EexpW = function(N,delta){
    EU = numeric(N);                         # EU = EexpW
    EexpSplus = 2*pnorm(delta*sqrt(1:N)/2);  
    EU[1] = EexpSplus[1];
    for (n in 1:(N-1)){ 
        EU[n+1] = sum(EU[n:1] * EexpSplus[1:n])/n;
        }
    return(EU);  }  # end of function EexpW

# Generic lower bound. We'll choose k, h later.
LBound = function(n,k,h,delta){   # lower bound for P(Wn > h)
  F = exp(h)*(1 - (pnorm(h,-k*delta^2/2,delta*sqrt(k)))^(n/k));
  return(F) }  # end of function LBound


N = 1000; delta = 1;
n = 1:N; h = (log(n)); k = 2*h/delta^2;
plot( n, EexpW(N,delta), type="l", lwd=2, col="red" )
lines( n, LBound(n,k,h,delta), lwd=2, col="blue" )
lines( n, n*delta^2/4/h^(3/2), col="green" )




####################### Earlier ################################

n=1000; h = sqrt(log(n)/1.2); 
delta = seq(0.01,1,0.01);
F = LBound(n,1,h,delta);
plot(delta, F, type="l", ylim=c(0,2*max(F)) )
for (k in 2:25){
  F = LBound(n,k,h,delta);
  lines(delta, F, col=k );
  text(delta[which.max(F)]-0.5, 1.1*max(F), paste("k =",k), col=k) }
  F = LBound(n,n,h,delta);
  lines(delta, F );
  text(delta[which.max(F)]-0.5, 1.1*max(F), paste("k =",n))

