# Trying to get help from our lower bound, which gives a lower bound for EexpWn

delta = 1;
n = 1:10000;
# h = log(n);
# LowerBound = exp(h)*(1-(pnorm((h+delta^2/2)/delta))^n);
# plot( n, LowerBound, type="l" )

# plot( n, qnorm(exp(-1/n))-1/2, type="l" )
# lines( n, (log(log(n)))^(1.5), col="blue" )

# plot( n, n*dnorm(qnorm((exp(-1/n)))), type="l" )

# x = seq(0.9,1,0.0001)
# plot( x, qnorm(x), type="l" )
# lines( x, sqrt(-2*log(1-x)), col="blue" )

c = 1; delta = .5;
# h = (1-exp(-c))*exp(delta*sqrt(-2*log(1-exp(-c/n)))-delta^2/2)
# plot(n,h,type="l")

EexpWLB = n^(delta*sqrt(2/log(n)))
plot(n,EexpWLB,type="l") 