

#trt=1,2
#status=0 (censoring), 1 (death), 2 (rec event)
process.data=function(id,time,status,trt){



#event/death level data
data=as.matrix(cbind(id,time,status,trt))


data.e=data[status ==1|status ==2,]

var=paste(data.e[,3],data.e[,4],sep="")
q=length(var)
for (i in 1:q){
if (var[i]=="11"){var[i]="d1"}
else{if (var[i]=="21"){var[i]="e1"}
	else{if (var[i]=="12"){var[i]="d2"}
		else{if (var[i]=="22"){var[i]="e2"}
				}
			}
		}
	}
data1=as.matrix(table(data.e[,2],var))
t=as.numeric(rownames(data1))
data1=cbind(t,data1)
rownames(data1)=NULL

l=length(t)
#to get y

Ct1=time[status!=2&trt==1]
Ct2=time[status!=2&trt==2]

y1=rep(NULL,l)
y2=rep(NULL,l)
for (k in 1:l){
	
	y1[k]=sum(Ct1>=t[k])
	y2[k]=sum(Ct2>=t[k])
	}

data1=cbind(data1,y1,y2)

#subject level data

n1=length(unique(id[trt==1]))
n2=length(unique(id[trt==2]))

subj=subj.lev(data=data,data1=data1,group=1)
Y1=subj$Y
dM1=subj$dM
dMD1=subj$dMD
subj=subj.lev(data=data,data1=data1,group=2)
Y2=subj$Y
dM2=subj$dM
dMD2=subj$dMD

result=list(data1=data1,Y1=Y1,Y2=Y2,dM1=dM1,dM2=dM2,dMD1=dMD1,dMD2=dMD2)
return(result)
}





subj.lev=function(data,data1,group){
data=data[trt==group,]
l=dim(data1)[1]
id=data[,"id"]
time=data[,"time"]
status=data[,"status"]
t=data1[,"t"]
if (group==1){
y=data1[,"y1"]
d=data1[,"d1"]
e=data1[,"e1"]
}else{
y=data1[,"y2"]
d=data1[,"d2"]
e=data1[,"e2"]
	}


uid=unique(id)
n=length(uid)#number of subject
dM=matrix(NA,n,l)
dMD=matrix(NA,n,l)
Y=matrix(NA,n,l)

time.data=cbind(id,time)
termin=time.data[status!=2,]
event=time.data[status==2,]
death=time.data[status==1,]

for (i in 1:n){

	for (k in 1:l){

		Y[i,k]=as.numeric(termin[,2][termin[,1]==uid[i]]>t[k]-(1e-12))
		if (sum(event[,1]==uid[i])==0){
			delta=0} else{
		delta=any(abs(event[,2][event[,1]==uid[i]]-t[k])<1e-12)}
		dM[i,k]=delta-Y[i,k]/y[k]*e[k]

		if (sum(death[,1]==uid[i])==0){
			gamma=0} else{
		gamma=(abs(death[,2][death[,1]==uid[i]]-t[k])<1e-12)}
		dMD[i,k]=gamma-Y[i,k]/y[k]*d[k]
		}

}

result=list(Y=Y,dM=dM,dMD=dMD)

return(result)
}








