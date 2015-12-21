## Conceptual model of how DA can work in Paleo

## 1 - models project one paleo time step ahead (e.g. ~ 50 yrs) using an ensemble that samples over sources of uncertainty (state, drivers, params)

par(mfrow=c(1,3))

par(lwd=2,cex=1.3)

Fseq = seq(0,1,length=40)
Fprior = dbeta(Fseq,2,3)

Cseq = seq(-100,100,by=5)
Cprior = dnorm(Cseq,-20,40)


#plot(Cseq,Cprior,type = 'l')
#plot(Fseq,Fprior,type='l')

## cflux given fconif and climate

Fsamp = rbeta(100,2,3)
Csamp = rnorm(100,200*Fsamp - 100,20)

plot(Fsamp,Csamp,xlab="Fraction Conifer",ylab="Carbon flux (g/m2/yr)")
##plot(density(Fsamp,from=0,to=1))
##plot(density(Csamp))

## 2 - some fraction of conifers is infered with uncertainty from the pollen data
Fseq = seq(0,1,length=100)
Flike = dbeta(Fseq,10,6)
Fprior = density(Fsamp,from=0,to=1,n=100)

Fpost = Flike * Fprior$y
Fpost = Fpost/sum(Fpost)/diff(Fseq)[1]

#plot(Fseq,Fprior$y,type='l',xlab="Fraction Conifer",ylab="Probability",lwd=4)

plot(Fseq,Flike,type='l',xlab="Fraction Conifer",ylab="Probability",lwd=4,ylim=range(c(Flike,Fprior$y,Fpost)))
lines(Fseq,Fprior$y,lty=2,lwd=4)
lines(Fseq,Fpost,col=2,lwd=4)
legend("topleft",legend=c("prior","data","posterior"),lty=c(2,1,1),col=c(1,1,2),lwd=4)

## 3 - preferentially weight samples compatible with data, update carbon PDF
Cprior = density(Csamp)
Cseq   = Cprior$x
Cprior = Cprior$y/sum(Cprior$y)/diff(Cseq)[1]
Fsamp2 = sample(Fseq,10000,prob=Fpost,replace=TRUE)
Cpost = density(rnorm(10000,200*Fsamp2 - 100,20),from=min(Cseq),to=max(Cseq))
Cpost = Cpost$y/sum(Cpost$y)/diff(Cseq)[1]

#plot(Cseq,Cprior,type='l',xlab="Carbon flux (g/m2/yr)",ylab="Probability",lwd=4)

plot(Cseq,Cprior,type='l',ylim=range(c(Cprior,Cpost)),xlab="Carbon flux (g/m2/yr)",ylab="Probability",lwd=4,lty=2)
lines(Cseq,Cpost,col=2,lwd=4)
legend("topleft",legend=c("prior","posterior"),lty=c(2,1),col=c(1,2),lwd=4)
