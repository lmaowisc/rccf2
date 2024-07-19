

#Compute KM estimates
KM=function(d,y){
l=length(d)
s=rep(NA,l)
p=1
for (i in 1:l){
Lambda=d[i]/y[i]
s[i]=p*(1-Lambda)
p=s[i]
	}
return(s)
}

#compute S(t-)
sminus=function(s){
return(c(1,s[-length(s)]))
}


#compute the mean frequency function miu(t) and Lambda(t)
meanfreq=function(s,d,e,y){

l=length(s)
du=s*e/y
dL=d/y
u=rep(NA,l)
L=rep(NA,l)
for (k in 1:l){
u[k]=sum(du[1:k])
L[k]=sum(dL[1:k])
	}

return(list(du=du,u=u,dL=dL,L=L))
}


#compute the Psi function

Psifun=function(s,y,dM,dMD,du){

n=nrow(dM)
l=ncol(dM)
dPsi=matrix(NA,n,l)
Psi=matrix(NA,n,l)
Psi1=matrix(NA,n,l)
for (k in 1:l){
	for (i in 1:n){
		dPsi[i,k]=n*s[k]/y[k]*dM[i,k]-n*sum(dMD[i,1:k]/y[1:k])*du[k]
		Psi[i,k]=sum(dPsi[i,1:k])
			}
	}



Xi=apply(Psi^2,2,mean)

return(list(dPsi=dPsi,Psi=Psi,Xi=Xi))
}


# compute the test statistic

test.stat=function(y1,y2,du1,du2,dL1,dL2,n1,n2,tau){


KLR=y1*y2*(n1+n2)/((y1+y2)*(n1*n2))
if(any(y1*y2==0)){
tau=min(tau,which(y1*y2==0)[1]-1)
}

QLR=sum(KLR[1:tau]*(du1-du2)[1:tau])
QD=sum(KLR[1:tau]*(dL1-dL2)[1:tau])

return(list(KLR=KLR,QLR=QLR,QD=QD))
}


#compute the variance function

sigmafun=function(KLR,dPsi1,dPsi2,dMD1,dMD2,y1,y2,tau){
if(any(y1*y2==0)){
tau=min(tau,which(y1*y2==0)[1]-1)
}
n1=nrow(dPsi1)
n2=nrow(dPsi2)
n=n1+n2
V1=matrix(NA,2,n1)
V2=matrix(NA,2,n2)

for (i in 1:n1){
	V1[1,i]=sum((KLR*dPsi1[i,])[1:tau])
	V1[2,i]=sum((n1*KLR/y1*dMD1[i,])[1:tau])
	}
for (i in 1:n2){
	V2[1,i]=sum((KLR*dPsi2[i,])[1:tau])
	V2[2,i]=sum((n2*KLR/y2*dMD2[i,])[1:tau])
	}
Sigma=n2/n*(1/n1)*V1%*%t(V1)+n1/n*(1/n2)*V2%*%t(V2)
return(Sigma)
}




