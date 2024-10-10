


?arma.spec()

# Causality: An ARMA model is causal if its autoregressive polynomial has all roots outside the unit circle. This means that past values influence future values but not vice versa.
# Invertibility: An ARMA model is invertible if its moving average polynomial has all roots outside the unit circle. This indicates that current and past errors affect future errors in a stable way.


arma.spec(ar = c(1.5, -.75), ma = c(-.8,.4), col=4, lwd=2)

par(mfrow=c(2,1))
sois  = mvspec(soi, spans=c(7,7), taper=.1, col=4, lwd=2)

rm(list = ls())

x1 = 2*cos(2*pi*1:100*5/100)  + 3*sin(2*pi*1:100*5/100)
x2 = 4*cos(2*pi*1:100*10/100) + 5*sin(2*pi*1:100*10/100)
x3 = 6*cos(2*pi*1:100*40/100) + 7*sin(2*pi*1:100*40/100)
x  = x1 + x2 + x3

mvspec(x, col=4, lwd=2, type='o', pch=20)

library(devtools)
install_github("nickpoison/astsa/astsa_build")
library(astsa)
