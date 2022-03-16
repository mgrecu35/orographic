import numpy as np
import matplotlib.pyplot as plt

nbin=80
Dm_b=np.exp(np.linspace(np.log(2),np.log(200),nbin+1))
Dm=0.5*(Dm_b[1:]+Dm_b[:-1])
rhow=1.0
m=np.pi*(Dm*1e-4)**3/6.0*rhow  # in grams
Dm0=10
m0=np.pi*(Dm0*1e-4)**3/6.0*rhow # in grams
L=1
Ng=3*L*(m/m0)**2*np.exp(-m/m0)

plt.semilogx(Dm,Ng)

Nc=Ng*np.log(Dm_b[1:]/Dm_b[:-1])
V=1e-3
Nd=Nc/m*V
b=1500 # in cm^3/gram
def kern_g(i,j,m,b,V):
    k=b*(m[i]+m[j])/V*1e-6 # to transform in s-1
    return k

K=np.zeros((nbin,nbin),float)
for i in range(nbin):
    for j in range(nbin):
        K[i,j]=kern_g(i,j,m,b,V)


from numba import jit
ijL=[]
for i in range(nbin):
    for j in range(i):
        ijL.append([i,j])
    ijL.append([i,i])
@jit(nopython=True)
def calcprob(S,K,Nd,nbin):
    ic=0
    s=0
    for i in range(nbin):
        for j in range(i):
            s=s+Nd[i]*Nd[j]*K[i,j]
            S[ic]=s
            ic+=1
        if Nd[i]>1:
            S[ic]=s+Nd[i]*(Nd[i]-1)/2*K[i,j]
        else:
            S[ic]=s
        ic+=1
    #print(ic)

S=np.zeros((int(nbin*(nbin+1)/2)),float)

Nd0=Nd.copy()

it=0
t=0

from bisectm import *
NdL=[]
kL=[]
tL=[]
nL=len(ijL)
calcprob(S,K,Nd,nbin)
St=S[-1]
while it<500*3600:
    if it%10==0:
        calcprob(S,K,Nd,nbin)
        St=S[-1]
        S=S/S[-1]
    r1=np.random.random()
    r2=np.random.random()
    dt=-np.log(r1)/St
    
    ni1,ni2=bisect2(S,nL,r2)
    i1,i2=ijL[ni2]
    if Nd[i1]>1 and Nd[i2]>1:
        dm_new=(Dm[i1]**3+Dm[i2]**3)**(1/3.0)
        if dm_new<Dm_b[nbin]:
            n1,n2=bisect2(Dm_b,nbin+1,dm_new)
        else:
            n1=nbin-1
        
        Nd[i1]-=1
        Nd[i2]-=1
        Nd[n1]+=(dm_new)**3/Dm[n1]**3

    t+=dt
    if int(t)%60==0:
        k=int(t)%60
        if k not in kL:
            NdL.append(Nd)
            kL.append(k)
    it+=1

