#main function

recurr=function(id,time,status,trt){

#sample size
n1=length(unique(id[trt==1]))
n2=length(unique(id[trt==2]))
n=n1+n2

tmax=1e9

#process data
r=process.data(id=id,time=time,status=status,trt=trt)

#extract group-level data
t=r$data1[,"t"]
if (any(t>tmax)){
tau=which(t>tmax)[1]-1
}else{tau=length(t)}
d1=r$data1[,"d1"]
d2=r$data1[,"d2"]
e1=r$data1[,"e1"]
e2=r$data1[,"e2"]
y1=r$data1[,"y1"]
y2=r$data1[,"y2"]

#extract subj-level data
dM1=r$dM1
dM2=r$dM2
dMD1=r$dMD1
dMD2=r$dMD2
Y1=r$Y1
Y2=r$Y2

#compute KM estimates
s1p=KM(d=d1,y=y1)
s2p=KM(d=d2,y=y2)
#compute S(t-)
s1=sminus(s=s1p)
s2=sminus(s=s2p)

#compute the mean frequency function miu(t) and death hazard Lambda^D(t)
m1=meanfreq(s=s1,d=d1,e=e1,y=y1)
m2=meanfreq(s=s2,d=d2,e=e2,y=y2)
u1=m1$u
u2=m2$u
dL1=m1$dL
dL2=m2$dL
du1=m1$du
du2=m2$du

#compute the Psi function
a=Psifun(s=s1,y=y1,dM=dM1,dMD=dMD1,du=du1)
Psi1=a$Psi
dPsi1=a$dPsi
Xi1=a$Xi
b=Psifun(s=s2,y=y2,dM=dM2,dMD=dMD2,du=du2)
Psi2=b$Psi
dPsi2=b$dPsi
Xi2=b$Xi

#compute the test statistic
q=test.stat(y1=y1,y2=y2,du1=du1,du2=du2,dL1=dL1,dL2=dL2,n1=n1,n2=n2,tau=tau)
KLR=q$KLR
QLR=q$QLR
QD=q$QD

stat=sqrt(n1*n2/n)*c(QLR,QD)
#compute the variance function
S=sigmafun(KLR=KLR,dPsi1=dPsi1,dPsi2=dPsi2,dMD1=dMD1,dMD2=dMD2,y1=y1,y2=y2,tau=tau)

return(list(stat=stat,S=S,u1=u1,u2=u2,t=t))

}





