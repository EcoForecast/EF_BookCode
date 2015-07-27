## Propagating, Analyzing, and Reducing Uncertainty

### Uncertainty vs Variability

par(mfrow=c(2,1))
X = rnorm(50)
xseq = seq(min(X),max(X),length=100)
plot(density(X),xlab="X",main="Variability",lwd=3)
points(X,rep(0,length(X)),pch="|")
mu = mean(X)
se = sd(X)/sqrt(length(X))
plot(xseq,dnorm(xseq,mu,se),type='l',xlab="E[X]",ylab="Density",main="Uncertainty",lwd=3)
abline(v=mu,col=2,lwd=3)
abline(v=mu+1.96*se*c(-1,1),col=2,lwd=2,lty=2)


### Sensitivity Analysis

### Group sampling
k = 128 ## number of parameters
ni = 4  ## number of important parameters
ng = 16  ## number of groups
ns = ceiling(k/ng)
nMC = 10000
correct = fp = fn = c2 = fp2 = rep(0.0,nMC) 
do.plot = FALSE
for(q in 1:nMC){

beta = runif(128,-0.5,0.5)
j = sample.int(k,ni)
beta[j] = sign(beta[j])*runif(ni,10,20)
#hist(beta)

X = rnorm(k)

### Group Sampling
stage1 = sample.int(k)
group1 = rep(1:ng,each=ns)[stage1]
Y1 = rep(NA,ng)
for(i in 1:ng){
  Xi = X
  Xi[group1 == i] = Xi[group1==i] + 1
  Y1[i] = sum(Xi*beta)
}
imp.g1 = abs(Y1-mean(Y1))>sd(Y1) ## important groups
imp.p1 = rep(0,k); imp.p1[group1 %in% which(imp.g1)] = 1
if(do.plot){
  par(mfrow=c(1,3))
  barplot(Y1-mean(Y1),names.arg = 1:ng,horiz=TRUE,cex.names=0.75,col=imp.g1)
}

stage2 = sample.int(k)
stage2 = c(which(imp.p1==1),stage2[-which(imp.p1==1)]) ##distrib imp. evenly
group2 = rep(1:ng,times=ns)[stage2]
Y2 = rep(NA,ng)
for(i in 1:ng){
  Xi = X
  Xi[group2 == i] = Xi[group2==i] + 1
  Y2[i] = sum(Xi*beta)
}
imp.g2 = abs(Y2-mean(Y2))>sd(Y2)
imp.p2 = rep(0,k); imp.p2[group2 %in% which(imp.g2)] = 1
if(do.plot){
  barplot(Y2-mean(Y2),names.arg = 1:ng,horiz=TRUE,cex.names=0.75,col=imp.g2)
}

stage3 = sample.int(k)
stage3 = c(which((imp.p1+imp.p2)==2),stage3[-which((imp.p1+imp.p2)==2)]) ##distrib imp. evenly
group3 = rep(1:ng,each=ns)[stage3]
Y3 = rep(NA,ng)
for(i in 1:ng){
  Xi = X
  Xi[group3 == i] = Xi[group3==i] + 1
  Y3[i] = sum(Xi*beta)
}
imp.g3 = abs(Y3-mean(Y3))>sd(Y3)
imp.p3 = rep(0,k); imp.p3[group3 %in% which(imp.g3)] = 1
if(do.plot){
barplot(Y3-mean(Y3),names.arg = 1:ng,horiz=TRUE,cex.names=0.75,col=imp.g3)
}
#par(mfrow=c(1,1))
#barplot(beta,col=(imp.p1+imp.p2+imp.p3)>2,horiz=TRUE)
if(do.plot){
table(imp.p1,abs(beta)>1)
table(imp.p1+imp.p2,abs(beta)>1)
table(imp.p1+imp.p2+imp.p3,abs(beta)>1)
}
correct[q] = sum((imp.p1+imp.p2+imp.p3) == 3 & abs(beta)>1)
c2[q] = sum((imp.p1+imp.p2+imp.p3) > 1 & abs(beta)>1)
fp[q] = sum((imp.p1+imp.p2+imp.p3) == 3 & abs(beta)<1)
fp2[q] = sum((imp.p1+imp.p2+imp.p3) >1 & abs(beta)<1)
fn[q] = sum((imp.p1+imp.p2+imp.p3) < 3 & abs(beta)>1)

}
table(correct)/nMC*100
table(fp)/nMC*100
table(c2)/nMC*100
table(fp2)/nMC*100


### UNCERTAINTY ANALYSIS: CASES
mu = 3
x = seq(0,15,length=1000)
sens = c(1,0.4)
sd = c(0.4,1)
ymax = mu*max(sens)+2*max(sd)

layout(matrix(c(1,2,3,
                1,4,5,
                6,7,7),3,3,byrow=TRUE),
       widths=c(0.12,0.44,0.44),height=c(0.44,0.44,0.12))
#layout.show(7)

#par(mfrow=c(2,2))
par(mar=c(0,0,0,0))
plot.new() #(0,0,type='n',bty='n',ann=FALSE)
text(0.75,0.25,"LOW",cex=2,srt=90)
text(0.75,0.78,"HIGH",cex=2,srt=90)
text(0.25,0.5,"SENSITIVITY",cex=3,srt=90)
par(bty='l')
par(mar=c(2,2,1,1))
for(i in 1:2){
  for(j in 1:2){
    y=sens[i]*x
    SDy = sd[j]*sens[i]
    plot(x,y,type='l',lwd=4,xlim=c(0,mu+2*max(sd)),ylim=c(0,ymax),cex.axis=1.2)           ## sensitivity
    lines(x,dnorm(x,mu,sd[j]),col=2,lwd=3)   ## input distribution
    lines(dnorm(y,sens[i]*mu,SDy),y,col=3,lwd=3)  ## output distribution
    ##intervals
    lines(c(mu,mu,0),c(0,mu*sens[i],mu*sens[i]),lty=2,col="grey")
    lines(c(mu+1.96*sd[j],mu+1.96*sd[j],0),c(0,mu*sens[i]+1.96*SDy,mu*sens[i]+1.96*SDy),lty=2,col="grey")
    lines(c(mu-1.96*sd[j],mu-1.96*sd[j],0),c(0,mu*sens[i]-1.96*SDy,mu*sens[i]-1.96*SDy),lty=2,col="grey")
    text(1,0.9*ymax,expression(SD[Y]),cex=2)
    text(1.36,0.9*ymax,paste("=",SDy),cex=2,pos=4)
  }
}
plot.new()
par(mar=c(0,0,0,0))
plot.new()
text(0.25,0.75,"LOW",cex=2)
text(0.8,0.78,"HIGH",cex=2)
text(0.5,0.25,"PARAMETER UNCERTAINTY",cex=3)

###############  TAYLOR SERIES: Michaelis-Menten ########################
require(mvtnorm)

V = 100
V.sd = 10
V.scale = 125
k.scale = 250
k = 30
k.sd = 10
Vk.cov = 15
x = 15
SIGMA = matrix(c(V.sd^2,Vk.cov,Vk.cov,k.sd^2),2,2)
V.seq = seq(V-3*V.sd,V+3*V.sd,length=1000)
k.seq = seq(k-3*k.sd,k+3*k.sd,length=1000)
x.seq = 0:60

par(mfrow=c(3,1),mar=c(4,5,1,1))
## Vmax
Y.V = V.seq*x/(k+x)
Y.V.mu = V*x/(k+x)
Y.V.sd = x/(k+x)*V.sd
Y.V.seq = seq(min(Y.V),max(Y.V),length=1000)
plot(V.seq,Y.V,xlab="V",ylab="Y",cex=1.2,cex.axis=1.2,cex.lab=1.6,
     type='l',lwd=3,bty='l',ylim=range(Y.V)) ## sensitivity
lines(V.seq,V.scale*dnorm(V.seq,V,V.sd)+min(Y.V),col=2,lwd=3)               ## input distribution
lines(V.scale*dnorm(Y.V.seq,Y.V.mu,Y.V.sd)+min(V.seq),Y.V.seq,col=3,lwd=3)  ## output distribution
lines(c(V,V,0),c(0,Y.V.mu,Y.V.mu),lty=2,col="grey")                         ## mean

## k
Y.k = V*x/(k.seq+x)
Y.k.mu = V*x/(k+x) + V*x/(k+x)^3*k.sd^2
Y.k.sd = V*x/(k+x)^2*k.sd
Y.rng = c(min(min(Y.k),Y.k.mu-2.5*Y.k.sd),max(max(Y.k),Y.k.mu+2.5*Y.k.sd))
Y.k.seq = seq(Y.rng[1],Y.rng[2],length=1000)
plot(k.seq,Y.k,xlab="k",ylab="Y",cex=1.2,cex.axis=1.2,ylim=Y.rng,cex.lab=1.6
     ,type='l',lwd=3,bty='l')
lines(k.seq,k.scale*dnorm(k.seq,k,k.sd)+Y.rng[1],col=2,lwd=3)               ## input distribution
lines(k.scale*dnorm(Y.k.seq,Y.k.mu,Y.k.sd)+min(k.seq),Y.k.seq,col=3,lwd=3)  ## output distribution
lines(c(k,k,0),c(0,Y.V.mu,Y.V.mu),lty=2,col="grey",lwd=1)                         ## naive mean
lines(c(k/2,0),c(Y.k.mu,Y.k.mu),lty=2,col=3,lwd=2)                         ## corrected mean


## Confidence Interval
mu.Y  = V*x.seq/(k+x.seq) + V*x.seq/(k+x.seq)^3*SIGMA[2,2]-x.seq/(k+x.seq)^2*SIGMA[1,2]
var.Y = (x.seq/(k+x.seq))^2*SIGMA[1,1]+(V*x.seq)^2/(k+x.seq)^4*SIGMA[2,2]-V*x.seq^2/(k+x.seq)^3*SIGMA[1,2]
CI = data.frame(
  fit = mu.Y,
  upr = mu.Y+1.96*sqrt(var.Y),
  lwr = mu.Y-1.96*sqrt(var.Y)
  )
plot(x.seq,V*x.seq/(k+x.seq),xlab="X",ylab="Y",cex=1.2,cex.axis=1.2,cex.lab=1.6
     ,type='l',lwd=3,bty='l',ylim=range(CI))
lines(x.seq,CI$fit,col=3,lwd=3,lty=2)
lines(x.seq,CI$lwr,col=3,lwd=3,lty=2)
lines(x.seq,CI$upr,col=3,lwd=3,lty=2)

## Partial Variance
plot(x.seq,100*(x.seq/(k+x.seq))^2*SIGMA[1,1]/var.Y
     ,xlab="X",ylab="% Var",cex=1.2,cex.axis=1.2,cex.lab=1.6
     ,type='l',lwd=3,bty='l',ylim=c(0,100))
lines(x.seq,100*(V*x.seq)^2/(k+x.seq)^4*SIGMA[2,2]/var.Y,col=2,lwd=3,lty=2)
lines(x.seq,100*abs(-V*x.seq^2/(k+x.seq)^3*SIGMA[1,2])/var.Y,col=3,lwd=3,lty=3)
legend(0,80,legend=c("V","k","Cov"),col=1:3,lty=1:3,cex=1.6,lwd=2)

## X
par(mfrow=c(1,1))
k = 2
V = 40
x.a = 2
x.b = 0.4
x.mu = x.a/x.b
x.sd = sqrt(x.a/x.b^2)
x.samp=rgamma(100000,x.a,x.b)
y.samp=V*x.samp/(k+x.samp)
y.dist = density(y.samp,n=1024)
Y.x.mu = V*x.mu/(k+x.mu) -2*V/(k+x)^2*(1-x/(k+x))*k.sd^2
m = sqrt(V/(k+x.mu)*(1-x.mu/(k+x.mu)))
Y.x.sd = m*k.sd

x.seq=seq(-10,20,length=1000)
x.scale=150
y.rng = c(-10,45)
Y.x.seq = seq(y.rng[1],y.rng[2]*1.5,length=1000)

Y = V*x.seq/(k+x.seq)
Y[x.seq <= -k] = NA
plot(x.seq,Y,xlab="X",ylab="Y",cex=1.2,cex.axis=1.2,cex.lab=1.6
     ,type='l',lwd=4,bty='l',ylim=y.rng)
lines(x.seq,x.scale*dgamma(x.seq,x.a,x.b)+y.rng[1],col=2,lwd=3)                             ## input distribution - exact
lines(x.seq,x.scale*dnorm(x.seq,x.mu,x.sd)+y.rng[1],col=2,lwd=3,lty=2)                      ## input distribution - approx
lines(x.scale*y.dist$y+min(x.seq),y.dist$x,col=3,lwd=3)                                     ## output distribution - mc
lines(x.scale*dnorm(Y.x.seq,Y.x.mu,Y.x.sd)+min(x.seq),Y.x.seq,col=3,lwd=3,lty=2)  ## output distribution
#lines(c(x.mu,x.mu,y.rng[1]),c(y.rng[1],Y.x.mu,Y.x.mu),lty=2,col="grey",lwd=1)                         ## naive mean
#lines(c(min(x.seq),0),rep(mean(y.samp),2),col=3,lwd=2)                         ## corrected mean
#abline(h=0,lty=3)
#abline(v=0,lty=3)
abline(V*x.mu/(k+x.mu)-m*x.mu,m,lty=2)

#CI.y = quantile(y.samp,c(0.025,0.975))
#lines(c(-10,-5),rep(CI.y[1],2),col=3)
#lines(c(-10,-5),rep(CI.y[2],2),col=3)
#lines(c(-10,-5),rep(qnorm(0.025,Y.x.mu,Y.x.sd),2),col=3,lty=2)
#lines(c(-10,-5),rep(qnorm(0.975,Y.x.mu,Y.x.sd),2),col=3,lty=2)
1-pnorm(V,Y.x.mu,Y.x.sd)

### Regression MC uncertainty propagation

xlim = c(0,10) ## range of data
Btrue = c(10,-1.5) ## true regression parameters
Strue = 2 ## true RMSE

# Pseudo-data
N = 10
X = runif(N,xlim[1],xlim[2])
Y = rnorm(N,Btrue[1]+Btrue[2]*X,Strue)
plot(X,Y)
fit = lm(Y~X); abline(fit,col=2)
B = coef(fit)
S = vcov(fit)

## MC uncertainty
library(mvtnorm)
Nmc = 10000
Nbin = 20
xseq = seq(xlim[1],xlim[2],length=Nbin)
MC = matrix(0.0,Nmc,Nbin)
Bmc = matrix(0,Nmc,2)
for(i in 1:Nmc){
  # draw regression parameters
  Bi = as.vector(rmvnorm(1,B,S))
  Bmc[i,] = Bi
  MC[i,] = Bi[1]+Bi[2]*xseq  
}

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
CI = apply(MC,2,quantile,c(0.025,0.5,0.975))

layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
#layout.show(3)
plot(xseq,CI[2,],type='n',ylim=range(CI),xlab="X",ylab="Y",cex.lab=1.5)
ciEnvelope(xseq,CI[1,],CI[3,],col="lightblue")
#lines(xseq,CI[2,],col=4,lty=2)
abline(Bmc[6,])
points(xseq,MC[6,],pch=18)
abline(v=xseq[14],lty=3,lwd=3)

plot(Bmc,cex=0.5,pch=18,xlab="Intercept",ylab="Slope",cex.lab=1.5)

hist(MC[,14],main=" ",xlab="Y",probability=TRUE,cex.lab=1.5)


### Power Analysis
source("Power.R")

## example data
nmax <-50
data <- list(Y=c(4,6),site=1:2,trt=c(0,0),ghs=c(0,0),obs.prec=c(1,1),n=c(6,4))		
data2 <- cbind(data$Y,data$site,data$trt,data$ghs,data$obs.prec,data$n)
colnames(data2) <- c("Y","site","trt","ghs","obs.prec","n")
data2 <- as.data.frame(data2)


## Case 1: Mean & Sin, nonparameteric bootstrap
SEprime2 <- power.MA(nmax,c(5,1),1,data=data2,method='boot',nrep=2000)

## Case 2: Transform to model output space: bootstrap
gfun <- splinefun(1:10,pnorm(1:10,7,3))
SEprime4 <- power.MA(nmax,c(5,1),1,data=data2,method='ana',nrep=2000,gfun = gfun)
SEprime5 <- power.MA(nmax,c(5,1),1,data=data2,method='boot',nrep=2000,gfun = gfun)

## Case 3: Marginal optimization
SE4opt <- list(theta1=SEprime4$Y,theta2=SEprime5$Y)
SEopt <- sample.optim(SE4opt)

## Case 4: Marginal cost optimization
cost <- c(1,15)
SEopt2 <- sample.optim(SE4opt,cost)


par(mfrow=c(2,2))
plot(SEprime2,xlab="Sample Size",ylab="Parameter Standard Error",lwd=3,type='l',cex=1.2,cex.lab=1.4)

plot(SE4opt$theta1,ylim=range(sapply(SE4opt,range)),type='l',col=1,lwd=3,xlab="Sample Size",ylab="Model Standard Error",cex.lab=1.4)
lines(SE4opt$theta2,lwd=3,col=2,lty=2)
legend("topright",legend=c("var1","var2"),col=c(1,2), lty=1:2,lwd=3)

sel = 2:50
plot(SEopt$SE[sel],col=SEopt$varID[sel],pch=17+SEopt$varID[sel],xlab="Sample Size",ylab="Model Standard Error",type='p',cex.lab=1.4)
legend("topright",legend=c("var1","var2"),col=c(1,2),pch=18:19)

sel = 2:65
plot(cumsum(SEopt2$cost[sel]),SEopt2$SE[sel],col=SEopt2$varID[sel],
     pch=17+SEopt2$varID[sel],xlim=c(0,150),xlab="Cumulative Cost ($)",
     ylab="Model Error",type='p',cex.lab=1.4)
legend("topright",legend=c("var1 $1","var2 $15"),col=c(1,2),pch=18:19)


