V=100
nbins=80
from netCDF4 import Dataset
import numpy as np
def kern(i,j,d,vfall,ecoal,V):
    k=np.pi*((d[i]+d[j])/2)**2*np.abs(vfall[i]-vfall[j])/V
    return k

dD=0.05
D=np.arange(nbins)*dD+dD/2.0
Dm_b=np.exp(np.linspace(np.log(0.1),np.log(16),nbins+1))
Dm=0.5*(Dm_b[1:]+Dm_b[:-1])

def readScatProfR(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    bscat=fh['bscat'][:]*4*np.pi
    Deq=10*(mass*1e3*6/np.pi)**(0.333) # in mm
    ext=fh['ext'][:]
    vfall=fh['fall_speed'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    #print(fh)
    #stop
    return temp,mass,bscat,Deq,ext,scat,g,vfall

freq=13.8
fnameRain='../scatter-1.1/test/scattering_input/liquid-water_%06.2f-GHz_scat.nc'%freq

    
temp_r,mass_r,bscat_r,Deq_r,\
    ext_r,scat_r,g_r,vfall_r=readScatProfR(fnameRain)

vfall=np.interp(Dm,Deq_r[:],vfall_r[:])

K=np.zeros((nbins,nbins),float)
ecoal=1
Dm_m=Dm*1e-3
for i in range(nbins):
    for j in range(nbins):
        K[i,j]=kern(i,j,Dm_m,vfall,ecoal,V)


def dm_lwc(nw,rwc,rho):
    dm=(4**4*rwc/(np.pi*1e6*Nw))**0.25*1e3 # in mm
    return dm

mu=2.0
Nw=0.08e8
lwc=1  #g/m^3
rhow=1 #g/cm^3

dm=dm_lwc(Nw,lwc,rhow) # in mm
lambd=(4+mu)/dm
from scipy.special import gamma as gam
def fmu(mu):
    return 6/4**4*(4+mu)**(mu+4)/gam(mu+4)

f_mu=fmu(mu)
dD=Dm_b[1:]-Dm_b[:-1]
Nd=Nw*f_mu*np.exp(-lambd*Dm)*(Dm/dm)**mu*dD*1e-3 #(m^3)
W=(Nd*(0.1*Dm)**3*np.pi/6*rhow).sum()
Nd=Nd*V
a=np.nonzero(K>0)

from bisectm import *


it=0
iint=0

Nd0=Nd.copy()

from numba import jit

#@jit(nopython=True)
def calcprob(S,K,Nd,nbins):
    ic=0
    s=0
    for i in range(nbins):
        for j in range(i):
            s=s+Nd[i]*Nd[j]*K[i,j]
            S[ic]=s
            ic+=1
        if Nd[i]>1:
            S[ic]=s+Nd[i]*(Nd[i]-1)/2*K[i,i]
        else:
            S[ic]=s
        ic+=1
    #print(ic)
    

S=np.zeros((int(nbins*(nbins+1)/2)),float)

Nd0=Nd.copy()

it=0
t=0
ijL=[]
for i in range(nbins):
    for j in range(i):
        ijL.append([i,j])
    ijL.append([i,i])
    
from bisectm import *
NdL=[]
kL=[]
tL=[]
nL=len(ijL)
calcprob(S,K,Nd,nbins)
St=S[-1]
while it<180*3600 and t<1800:
    if it%10==0:
        calcprob(S,K,Nd,nbins)
        St=S[-1]
        S=S/S[-1]
    r1=np.random.random()
    r2=np.random.random()
    dt=-np.log(r1)/St
    if dt>10:
        stop
    ni1,ni2=bisect2(S,nL,r2)
    i1,i2=ijL[ni2]
    if Nd[i1]>1 and Nd[i2]>1:
        dm_new=(Dm[i1]**3+Dm[i2]**3)**(1/3.0)
        if dm_new<Dm_b[nbins]:
            n1,n2=bisect2(Dm_b,nbins+1,dm_new)
        else:
            n1=nbins-1
        
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

