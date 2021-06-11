#R code for mobile window
system('sudo apt-get install libssl-dev libxml2-dev libcurl4-openssl-dev -y
R')
# Load data of ACE from Colorado State University
package=c('XML')
for(i in 1:length(package)){if(!require(package[i],character.only=T)) {install.packages(package[i],dependencies=TRUE); require(package[i], character.only=T)}}
web ='http://tropical.atmos.colostate.edu/Realtime/index.php?arch&loc=northatlantic'
ace=readHTMLTable(web,as.data.frame =TRUE, colClasses='numeric')$info
str(ace)
#Variables, choose ordering or ranking, edit red fields
s=ts(ace$'Accumulated Cyclone Energy',start=min(ace$Year)) # time series 
names(s)=time(s) #labels
#lab = 'Order'; x=sort(s)	# Ordering
lab = 'Rank'; x=sort(s,decreasing=TRUE) #ranking
lx = log(x)
ex = exp(x)
N = length(x)
n = 1:N; ln=log(n) # sorting cardinality
m = max(n)+min(n)-n; lm=log(m) # convolute cardinality
title = paste('Annual ACE from ', min(time(s)),' to ', max(time(s)))
variable = 'ACE'
slab=expression('ACE ('*10^4*kt^2*')')
credits=c('Elio Roca Flores elioroca@gmail.com 2020 PCT/CCA/UNAM. Data: Colorado State University')
credits=c('')
# Plot time series
par(mar=c(5.1,5.1,4.1,1.1))
plot(s,type='b',xlab='Year',ylab=slab, main=title,las=1)
mtext(credits,adj=1,col='darkgray',cex=.7)
#Hurst
#El exponente de Hurst H permite determinar la correlación o memoria de los datos (Simonsen, 1998; Jones, 1996). De 0.5 a 1 indica persistencia o correlación positiva, de 0 a 0.5 alternancia o correlación negativa y 0.5 implica aleatoriedad o no correlación (Olvera, 2009). 
package=c('pracma');for(i in 1:length(package)){if(!require(package[i],character.only=T)) {install.packages(package[i],dependencies=TRUE); require(package[i], character.only=T)}}
hurstexp(s)
#Simple R/S Hurst estimation:         0.5855893 
#Corrected R over S Hurst exponent:   0.6920908 
#Empirical Hurst exponent:            0.2930231 
#Corrected empirical Hurst exponent:  0.2715825 
#Theoretical Hurst exponent:          0.5597524 
#Plot shorting
layout(matrix(1:2,1,2));plot(x, las=1);plot(lx,ln,las=1)
title(title, line = -3, outer = TRUE)
mtext(credits,adj=1,col='darkgray',cex=.7)
#Mobile plot
#plot
running.plot=function(c,tu='years',title='',subtitle='Mobile window',legend=slab)
{
require(fields)
image.plot(as.matrix(c),xlab=paste('Time (',tu,')'),ylab=paste('Window size (',tu,')'), main=c(title,'\n',subtitle),xaxt= 'n', yaxt= 'n',legend.lab=legend, legend.line=3)
axis(1,labels=round(seq(min(time(s)),max(time(s)),length.out =11),0),at=seq(0,1,length.out =11), las=2, cex.axis=.7)
#axis(2,labels=round(seq(3,length(s),length.out =6),0),at=seq(2/length(s),1,(1-2/length(s))/5), las=2,cex.axis=1)
mtext(credits,adj=1,col='gray',cex=1)
}
#Mobile mean
require(gtools)
c=running(s,width=1)
if(length(s)%%2==1){H=length(s)/2} else{H=length(s)/2-1} #remove 1 for even length(s)
for(n in 1:H)
{
c=data.frame(c,c(rep(NA,n),running(s,width=n*2+1),rep(NA,n)))
}
#m=as.matrix(c) # data in matrix form
running.plot(c,tu='years',title=title,subtitle='Running mean',legend=slab)
setEPS()
postscript("Fig4a.eps")
dev.off()
#Mobile standard deviation
c=running(s,fun=sd,width=1)
if(length(s)%%2==1){H=length(s)/2} else{H=length(s)/2-1} #remove 1 for even length(s)
for(n in 1:H)
{
c=data.frame(c,c(rep(NA,n),running(s,fun=sd,width=n*2+1),rep(NA,n)))
}
#sd=as.matrix(c) # data in matrix form
#plot
running.plot(c,tu='years',title=title,subtitle='Standard deviation',legend=slab)
setEPS()
postscript("Fig4b.eps")
running.plot(c,tu='years',title='',subtitle='',legend=slab)
dev.off()
#Mobile linear fit parameters, choose fit in red, resolves the last uncommented. Nonlinear fit do no works for mobile window.
par(mar=c(5.1,5.1,3.1,1.1),oma=c(0,1,0,2))
layout(matrix(1:4,2))
#power a
require(gtools)
beta.a=function(X){
N=length(X)
x=sort(X,decreasing=TRUE) #ranking
n=1:N
fit=lm( log(x) ~ log(n) + log(N+1-n) )		#PMP x=cn^a(N+1-n)^b 
#ln2=(log(1:length(x)))^2; fit=lm( log(x) ~ ln2 + log(N+1-n) )	#NMP #x=cexp(a*log(n)^2)(N+1-n)^b
return(fit$coefficients[2])
}
c=running(s,fun=beta.a, width=1)
if(length(s)%%2==1){H=length(s)/2} else{H=length(s)/2-1} #remove 1 for even length(s)
for(n in 1:H)
{
c=data.frame(c,c(rep(NA,n),running(x,fun=beta.a,width=n*2+1),rep(NA,n))) # odd length as 2021 for AN
#c=cbind(c,c(rep(NA,n),running(x,fun=beta.a,width=n*2+1),rep(NA,n))) # even length
}
#plot
running.plot(c,tu='years',title=title,subtitle='Parameter a',legend='')
setEPS()
postscript("Fig4c.eps")
running.plot(c,tu='years',title='',subtitle='',legend='')
dev.off()

#power b
beta.b=function(X){
N=length(X)
x=sort(X,decreasing=TRUE) #ranking
n=1:N
fit=lm( log(x) ~ log(n) + log(N+1-n) )		#PMP x=cn^a(N+1-n)^b 
#ln2=(log(1:length(x)))^2; fit=lm( log(x) ~ ln2 + log(N+1-n) )	#NMP #x=cexp(a*log(n)^2)(N+1-n)^b
return(fit$coefficients[3])
}
require(gtools)
c=running(s,fun=beta.b, width=1)
for(n in 1:floor(length(x)/2) )
{
c=cbind(c,c(rep(NA,n),running(x,fun=beta.b,width=n*2+1),rep(NA,n)))
}
#plot
running.plot(c,tu='years',title=title,subtitle='Parameter b',legend='')
setEPS()
postscript("Fig4d.eps")
running.plot(c,tu='years',title='',subtitle='',legend='')
dev.off()
#Mobile correlation
if(!require('data.table')){install.packages('data.table',dependencies=TRUE);library('data.table')}
#mei data
url='https://psl.noaa.gov/enso/mei/data/meiv2.data'
data<-fread(url, skip=1)
#annual mean
mei=data.frame(data, mean=round(rowMeans(data[,2:13]),2))
smei=ts(mei[,13],start=min(mei[,1])); names(smei)=time(smei)
#subset data
s=ts(s[(length(s)-length(smei)+1):length(s)],start=min(time(smei)) )
c=running(s,smei,fun=cor,width=1)
if(length(s)%%2==1){H=length(smei)/2} else{H=length(smei)/2-1} #remove 1 for even length(s)
for(n in 1:H)
{
c=data.frame(c,c(rep(NA,n),running(s,smei,fun=cor,width=n*2+1),rep(NA,n)))
}
#plot
running.plot(c,tu='year',title=title,subtitle='Pearson’s correlation',legend='Correlation')

Fig. Correlación móvil entre ACE anual para el Atlántico Norte y MEI 1979-2020.
