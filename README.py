# ranking
Python code to rank annual ACE

Appendix. Python code
#install packages sudo apt-get install python3-scipy python3-matplotlib
import os
os.system('sudo apt-get install python-scipy python-matplotlib python-pandas') 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
#data 
import pandas as pd
web ='http://tropical.atmos.colostate.edu/Realtime/index.php?arch&loc=northatlantic'
#web ='http://tropical.atmos.colostate.edu/Realtime/index.php?arch&loc=northeastpacific'
table = pd.read_html(web, header=0)
print(table[0])
data=table[0] #select first table
X=data["Accumulated Cyclone Energy"]
year=data["Year"]
plt.scatter(year,X,facecolors='none', edgecolors='k')
plt.plot(year,X,  '-k', linewidth=1)
plt.xlabel('Year')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.show(block=False)
# Rank
x=sorted(X, reverse=True)
N=len(x)
n=np.linspace(1, N, N)
#plot
plt.figure()
#linear plot
plt.subplot(1,2,1)
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')
plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.title('Linear plot')
#logarithmic plot
plt.subplot(1,2,2)
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')
plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.title('Log-log plot')
plt.xscale('log');plt.yscale('log') #comment for linear plot
plt.show(block=False)
#Cocho fit 
def PMP(n,a,b,c): 
	return c*n**a*(N+1-n)**b

def NMP(n, a, b, c):
 return c * np.exp( a*(np.log(n))**2  ) * (N+1-n)**b

#remove first 26 years 1851-1877
#x=sorted(X[26:N], reverse=True);N=len(x);n=np.linspace(1, N, N)

#fit
popt, pcov = curve_fit(PMP, n, x) 	#popt[0] = a , popt[1] = c, nonlinear least squares
np.corrcoef(x,PMP(n, *popt) ) 
chisquare(x, f_exp=PMP(n, *popt))
plt.figure()
#linear plot
plt.subplot(1,2,1)
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')
popt, pcov = curve_fit(NMP, n, x)
plt.plot(n, NMP(n, *popt), 'r.', label=r'$c*exp( a*(np.log(n))^2)*(N+1-n)^b$')
popt, pcov = curve_fit(PMP, n, x)
plt.plot(n, PMP(n, *popt), 'b-', label=r'$c*n^a*(N+1-n)^b$')
plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.legend()

#logarithmic plot
plt.subplot(1,2,2)
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')
popt, pcov = curve_fit(NMP, n, x)
plt.plot(n, NMP(n, *popt), 'r.', label=r'$c*exp( a*(np.log(n))^2)*(N+1-n)^b$')
popt, pcov = curve_fit(PMP, n, x)
plt.plot(n, PMP(n, *popt), 'b-', label=r'$c*n^a*(N+1-n)^b$')
plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.legend()
plt.xscale('log');plt.yscale('log')
plt.show(block=False)

#Fit both tails. Arbitrary sing parameters
def PMP(n, a, b, c):
 return c * n**a  * (N+1-n)**b

def LSL(n, a, b, c):
 return c  + a*np.log(n)  + b*np.log(N+1-n)

def ESL(n, a, b, c):
 return c * np.exp( -a*n )  + b*np.log(N+1-n)

def NMP(n, a, b, c):
 return c * np.exp( a*(np.log(n))**2  ) * (N+1-n)**b

def NML(n, a, b, c):			#2 parameters
 return c * np.exp( a*(np.log(n))**2  ) * np.log(N+1-n) 

def ZMP(n, a, b, c):
 return c * np.log((N+1)/n)**a  * (N+1-n)**b 

def ZSL(n, a, b, c):
 return c * np.log((N+1)/n)**a  + b*np.log(N+1-n)

function = {'PMP':PMP,'LSL':LSL,'ESL':ESL,'NMP':NMP,'NML':NML,'ZMP':ZMP,'ZSL':ZSL}
name = ['PMP','LSL','ESL','NMP','NML','ZMP','ZSL']
#Kolmogorov-Smirnov
from scipy.stats import ks_2samp
for i in range(len(name)):
 f=function[name[i]]
 popt, pcov= curve_fit(f, n, x)
 name[i],round(ks_2samp(x, f(n, *popt))[0],3),round(ks_2samp(x, f(n, *popt))[1],3)

def fit(f,n,x,line,name):
 popt, pcov = curve_fit(f, n, x)
 plt.plot(n, f(n, *popt), line, label=name)
 #print(popt)
 for m in range(len(popt)):
  print(round(popt[m],2))
 return np.corrcoef(x,f(n, *popt) )		#Pearson correlation
# return chisquare(x, f_exp=f(n, *popt))
# return np.square(np.subtract(x,f(n, *popt))).mean() #MSE

#plot
line=['k:','k-','r:','b:','b-','y:','y-','b--','k-.']
name=['PMP','LSL','ESL','NML','NML','ZMP','ZSL']
#name=[r'$c * n^a  * (N+1-n)^b$',r'$c  + a*np.log(n)  + b*np.log(N+1-n)$',r'$c * np.exp( -a*n )  + b*np.log(N+1-n)$',r'$c * np.exp( a*(np.log(n))^2  )  * (N+1-n)^b$',r'$c * np.exp( a*(np.log(n))^2  )  *np.log(N+1-n)$',r'$ c * np.log((N+1)/n)^a  * (N+1-n)^b$',r'$c * np.log((N+1)/n)^a  + b*np.log(N+1-n)$']

plt.figure()
#plt.subplot(1,2,1)
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')
				#corr 1851-2018	corr 1851-2019	coef a	b	c
fit(PMP,n,x,line[0],name[0])		#.9907		.9916		-0.2	0.6	13.44
fit(LSL,n,x,line[1],name[1])		#.9946		.9952		-48.03	9.96	246.61
fit(ESL,n,x,line[2],name[2])		#.9965		.9961		0.03	12.58	178.09
fit(NMP,n,x,line[3],name[3])		#.9983		.9984		-0.05	0.35	43.0
fit(NML,n,x,line[4],name[4])		#.9978		.9980		-0.06	1.0	51.13
fit(ZMP,n,x,line[5],name[5]) 		#.9954		.9960		0.75	-0.29	365.86
fit(ZSL,n,x,line[6],name[6])		#.9955		.9961		0.87	7.73	59.85
plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.legend()
plt.title('Linear plot')
plt.show(block=False)

#plt.subplot(1,2,2)
plt.figure()
plt.scatter(n,x,facecolors='none', edgecolors='k', label='Experimental')

fit(PMP,n,x,line[0],name[0])		#cor .9907
fit(LSL,n,x,line[1],name[1])		#.9946
fit(ESL,n,x,line[2],name[2])		#.9965
fit(NMP,n,x,line[3],name[3])		#.9983
fit(NML,n,x,line[4],name[4])		#.9978 2nd place with only two parameters
fit(ZMP,n,x,line[5],name[5]) 		#.9954
fit(ZSL,n,x,line[6],name[6])		#.9955

plt.xlabel('Rank')
plt.ylabel(r'$ACE (10^4 kt^2)$')
plt.legend()
plt.title('Log-log plot')
plt.xscale('log');plt.yscale('log')
plt.show(block=False)

#Mobile full window
#in console install: pip install more-itertools
from more_itertools import substrings
#X=range(6);N=6  #test mode
w=list(substrings(X))
#function
u=list(range(len(w)))
for i in range(len(w)):
 u[i]=np.mean(w[i])

#fill gaps with np.nan
for i in range(N-1):
 j=0
 while j<(i+1):
  u.insert((2+i)*N-(1+i),np.nan)
  j+=1

#down side matrix
m=[]
while u != []:
  m.append(u[:N])
  u = u[N:]

# upside matrix
#m=[]
#while u != []:
#  m.insert(0,u[:N])
#  u = u[N:]

#plot
import matplotlib.pyplot as plt 
plt.xlabel('Year')
plt.xticks(range(0,N,20), range(year.min(),year.max(),20)) #add xlabel
plt.ylabel('Window(years)')
plt.title('Mobile window\n Mean')
plt.imshow(m); 
ax = plt.gca(); ax.set_ylim(ax.get_ylim()[::-1]) #invert y axis
plt.colorbar()
plt.show(block=False)


#Mobile odd window. mean
import os
os.system('sudo apt-get install python-more-itertools') 
from more_itertools import substrings
w=list(substrings(X))
u=list(range(len(w)))
for i in range(len(w)):
 u[i]=np.mean(w[i])

#list of odd elements centered

s=range(N/2-1)
for i in s:
 j=1
 while j<(N-2*i):
  u.pop((1+i)*N)
  j+=1
 
 j=1
 while j<i+2:
  u.insert((1+i)*N,np.nan)
  j+=1
 
 j=1
 while j<i+2:
  u.insert((2+i)*N-1-i,np.nan)
  j+=1

u.pop((s[len(s)-1]+2)*N)
u.pop((s[len(s)-1]+2)*N)
u.pop((s[len(s)-1]+2)*N)

# upside matrix
m=[]
while u != []:
  m.insert(0,u[:N])
  u = u[N:]

#plot
import matplotlib.pyplot as plt 
plt.figure(figsize=(12,5))
plt.xlabel('Year')
plt.xticks(np.linspace(0,N,num=10), np.linspace(year.min(),year.max(),num=10).round(0).astype(int) ) #add xlabel
plt.ylabel('Window(years)')
plt.yticks(np.linspace(0,N/2,num=10), np.linspace(N-2,1,num=10).round(0).astype(int) ) #add ylabel
plt.title('Mobile window\n Mean')
plt.imshow(m); 
#ax = plt.gca(); ax.set_ylim(ax.get_ylim()[::-1]) #invert y axis
plt.colorbar()
plt.show(block=False)

#Mobile odd window. Fit coefficient
import os
os.system('sudo apt-get install python-more-itertools') 
#os.system('pip install more_itertools')
from more_itertools import substrings

#coef a
def fa(f,X):
 N=len(X)
 x=sorted(X, reverse=True)
 n=np.linspace(1, N, N) 
 popt, pcov = curve_fit(f, n, x)
 return popt[0]	#0 for a, 1 for b,2 for c and last for k

w=list(substrings(x))
u=list(range(len(w)))
for i in range(len(w)):
 try:
  u[i]=fa(LSL,w[i])
 except:
  u[i]=np.nan

#list of odd elements centered

s=range(N/2-1)
for i in s:
 j=1
 while j<(N-2*i):
  u.pop((1+i)*N)
  j+=1
 
 j=1
 while j<i+2:
  u.insert((1+i)*N,np.nan)
  j+=1
 
 j=1
 while j<i+2:
  u.insert((2+i)*N-1-i,np.nan)
  j+=1

u.pop((s[len(s)-1]+2)*N)
u.pop((s[len(s)-1]+2)*N)
u.pop((s[len(s)-1]+2)*N)

# upside matrix
m=[]
while u != []:
  m.insert(0,u[:N])
  u = u[N:]

#plot
import matplotlib.pyplot as plt 
plt.figure(figsize=(12,4))
plt.xlabel('Year')
plt.xticks(np.linspace(0,N-1,num=10).round(0), np.linspace(year.min(),year.max(),num=10).round(0).astype(int) ) #add xlabel
plt.ylabel('Window(years)')
plt.yticks(np.linspace(0,N/2,num=10), np.linspace(N-2,1,num=10).round(0).astype(int) ) #add ylabel
#plt.title(titulo)
plt.imshow(m); 		#full matrix
#plt.imshow(m[0:(len(m)-10)]); 	#without small subsets
#ax = plt.gca(); ax.set_ylim(ax.get_ylim()[::-1]) #invert y axis
plt.colorbar()
plt.show(block=False)
#Fig 6 b vs a
from more_itertools import substrings
 
f=PMP
#f=LSL
titulo='Mobile window\n {} coef a'.format(f.__name__)
#coef a
def fa(f,X):
 N=len(X)
 x=sorted(X, reverse=True)
 n=np.linspace(1, N, N) 
 popt, pcov = curve_fit(f, n, x)
 return popt[0]		#0 for a, 1 for b,2 for c and last for k

w=list(substrings(X))
ua=list(range(len(w)))
for i in range(len(w)):
 try:
  ua[i]=fa(f,w[i])
 except:
  ua[i]=np.nan

#coef b
def fb(f,X):
 N=len(X)
 x=sorted(X, reverse=True)
 n=np.linspace(1, N, N) 
 popt, pcov = curve_fit(f, n, x)
 return popt[1]		#0 for a, 1 for b,2 for c and last for k

w=list(substrings(X))
ub=list(range(len(w)))
for i in range(len(w)):
 try:
  ub[i]=fb(f,w[i])
 except:
  ub[i]=np.nan

#plot a vs b
plt.figure(figsize=(7,5))
plt.rcParams['font.size'] = '15'
plt.plot(ua,ub,'k.')
#plt.plot(ua[1000:],ub[1000:],'k.')
#plt.plot(ua[-1000:],ub[-1000:],'k.')
plt.scatter(ua[-1],ub[-1],c='red', marker='o')
plt.xlabel('a')
plt.ylabel('b')
plt.show(block=False)
plt.savefig('f6.eps', bbox_inches='tight')
plt.close()
#Correlation
from scipy.stats import pearsonr
#nas = np.logical_or(np.isnan(ua), np.isnan(ub)) 
#pearsonr(ua[~nas], ub[~nas])
 
pearsonr(ua[1000:],ub[1000:])[0]
pearsonr(ua[-1000:],ub[-1000:])[0]
 
corr=list(range(len(w)))
for i in range(len(w)):
 try:
  corr[i]=pearsonr(ua[i:],ub[i:])[0]
 except:
  corr[i]=np.nan
 
#enter
plt.plot(corr,'k.')
#plt.plot(corr[1000:-3],'k.')
plt.xlabel('Window index')
plt.ylabel('Correlation')
plt.show(block=False)
ab
cd
e
#Fig. PMP function parameters a and b for North Atlantic annual ACE 1851 to 2020 time series. a) All windows for parameter a non positive c*(len(n)+1-n)**b * n**a. b) Parameter a positive: c*(len(n)+1-n)**b * n**-a. c) Windows bigger than 6 years ((170+164)*(170-164)/2=1002) ua[1000:] with Pearson coef of -0.7. d) Biggest 1000 windows or window bigger than 125=N-45 years ((1+45)*(45-1)/2=1012) ua[-1000:] with Pearson coef of -0.96. e) Correlation between parameter windows.
 

#Fig. LSL function parameters a and b for North Atlantic annual ACE 1851 to 2020 time series.
#appendix. Build step a step odd window
u=list(substrings(range(6)))
j=1
while j<N:
 u.pop(N)
 j+=1

j=1
while j<2:
 u.insert(N,np.nan)
 j+=1

j=1
while j<2:
 u.insert(2*N-1,np.nan)
 j+=1

j=1
while j<(N-2):
 u.pop(2*N)
 j+=1

j=1
while j<3:
 u.insert(2*N,np.nan)
 j+=1

j=1
while j<3:
 u.insert(3*N-2,np.nan)
 j+=1

u.pop(3*N)
#Appendix. Table comparison of coefficients
			a		b		c
PMP (Cocho)
Linear fit R		-0.2368       0.5091		23.71244         
Non linea rfit R	-0.2020  	0.6045 		13.4409
Non linear fit Python	-0.20201392   0.60447737  13.44087708
NMP (Lognormal)
Linear fit R		-0.03861      0.46561        21.39228
Non linear fit R	-0.04952  	0.34792 	43.00159  
Non linear fit Python	-0.04951566   0.34792209  43.00162643

#Differences in R between linear (black line) and nonlinear fit (red dots) for PMP Cocho DGBD function for the ranking of Atlantic ACE 1851-2019 (circles).
#Appendix. Redefining and varying coefficients.
from numpy import *
from matplotlib.pyplot import *

#function
def PMP(n,a,b,c): 
	return c*(len(n)+1-n)**b/n**a

def LSL(n, a, b, c):
 return c  + b*np.log(len(n)+1-n) - a*np.log(n) 

def ESL(n, a, b, c):
 return c * np.exp( -b*n )  + a*np.log(len(n)+1-n)

def PME(n, a, b, c):
 return c * np.exp( -b*n ) * a*np.log(len(n)+1-n)

def PMS(n, a, b, c):
 return c * np.exp( -n**b ) * a*np.log(len(n)+1-n)

def NMP(n, a, b, c):
 return c * np.exp( -1*b*(np.log(n))**2  ) * (len(n)+1-n)**a

def NML(n, a, b, c):
 return c * np.exp( -1*b*(np.log(n))**2  ) * np.log(len(n)+1-n)  #dos parámetros

def ZMP(n, a, b, c):
 return c * np.log((len(n)+1)/n)**a / (len(n)+1-n)**b 

def ZSL(n, a, b, c):
 return c * np.log((len(n)+1)/n)**a + b*np.log(len(n)+1-n)

f=NMP
#plot
figure()
a=1; b=1;c=10;n=arange(0,100,1)	
subplot(1,2,1)
title('Variando a y b=1')
for a in arange(0,10,0.5):
 x=f(n,a,b,c)
 plot(n,x)

a=1; b=1;c=10;n=arange(0,100,1)	
subplot(1,2,2)
title('Variando b y a=1')
for b in arange(0,10,0.5):
 x=f(n,a,b,c)
 plot(n,x)

show()
#a)b)
#Fig. Función Cocho variando los parámetros a (izquierda) y b (derecha). Para las expresiones a) c*n**a*(len(n)+1-n)**b y b) c*(len(n)+1-n)**b/n**a.

#Fig. Función NMP= c * np.exp( 1*a*(np.log(n))**2  ) * (len(n)+1-n)**b variando los parámetros a (izquierda) y b (derecha).

#Fig. Función NMP= c * np.exp( -1*b*(np.log(n))**2  ) * (len(n)+1-n)**a, variando los parámetros a (izquierda) y b (derecha).
#NMP
#LSL
#ESL
#PME#PMS
#NML
#ZMP
#ZSL
#Functional families
# f(x)= {log L ln(n), Gusein Zade Z log((N+1)/n)a, Power P na, geometric G an, exp E exp(bn), lognormal N exp(-bln(n)²), Weibull W exp(-anb), Gamma or factorial F an!} {n,N+1-n}
