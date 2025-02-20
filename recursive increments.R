### Recursive formula for the CUSUM exponential moments and their increments ###

# Normal distribution case. F = N(theta0, sigma), G = Normal(theta1, sigma).
# Let delta = (theta1-theta0)/sigma.
# Then Y = LLR = Normal( -delta^2/2, delta ) under F 
#  and Y = LLR = Normal( +delta^2/2, delta ) under G 

Lag = 1;      # 0 = EexpW = EU, 1 = D = 1st increments, 2 = D2 = double differences, etc.

setwd("C:\\Users\\baron\\Documents\\Research\\Malov\\Spitzer formula paper\\R")

N = 20; d = 0:N; 
 DELTA = c(0.1, 0.3, 0.5, 1.0, 2.0, 3.0); 
Ndelta = length(DELTA); 
EU = numeric(N+1);                       # First element is EU[0], then EU[1], etc.
EU[1] = 1;                               # EU[1] = EexpW[0] = Eexp0 = 1
M = pnorm(DELTA/2) - pnorm(-DELTA/2);

delta = max(DELTA);
if (Lag==0){ ymin=0; ymax = max(M)*N+1; ylabel = "Eexp(W_n)" }
if (Lag==1){ ymin=0; ymax = 5; ylabel = "Difference Delta_n" }
if (Lag==2){ ymin=-0.05; ymax = 0.01; ylabel = "2nd difference Delta2_n" }    
if (Lag==3){ ymin=-0.005; ymax = 0.01; ylabel = "3rd difference Delta3_n" }    
if (Lag>3){ ymax = 0.1; ylabel = "Higher difference DeltaR_n" }  

plot(c(0,N),c(0,0),type="l",ylim=c(ymin,ymax),xlab="n",
ylab=ylabel,main=paste(ylabel,"for ",expression(delta)," from",min(DELTA),"to",max(DELTA)))

# Recursive computation of EU = EexpW

Lag=2;

for (j in 1:Ndelta){
    delta = DELTA[j];
    EexpSplus = 2*pnorm(delta*sqrt(0:N)/2);  # First element is ESplus[1], ...
    for (n in 1:N){ 
        EU[n+1] = sum(EU[1:n] * EexpSplus[(n+1):2])/n;
        }
    if (j==1){ legend.text=paste("delta =",delta); } else {
       legend.text=c(legend.text,paste("delta =",delta)); }
    if (Lag==0){graph = EU; bound = M[j]*d + 1;}
    if (Lag==1){graph  = c(1,diff(EU)); bound = M[j]*rep(1,n+1);}
    if (Lag==2){graph  = c(1,1,diff(diff(EU))); bound = rep(1,n+1);}
    if (Lag==3){graph  = c(1,1,1,diff(diff(diff(EU)))); bound = rep(1,n+1);}
    lines(d,graph,col=j+2,lwd=2);
    lines(d,bound,col=j+2,lwd=2, lty=3);  # 3=dotted
  }

 legend(0, ymax, legend.text, col=3:Ndelta, lwd=2, cex=0.6, seg.len=9);

