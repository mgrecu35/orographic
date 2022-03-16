V=100
from netCDF4 import Dataset
import numpy as np
def kern(i,j,d,vfall,ecoal,V):
    k=np.pi*((d[i]+d[j])/2)**2*np.abs(vfall[i]-vfall[j])/V
    return k

dD=0.05
D=np.arange(160)*dD+dD/2.0
Dm_b=np.exp(np.linspace(np.log(D[0]),np.log(D[-1]),41))
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

K=np.zeros((40,40),float)
ecoal=1
Dm_m=Dm*1e-3
for i in range(40):
    for j in range(40):
        K[i,j]=kern(i,j,Dm_m,vfall,ecoal,V)


def dm_lwc(nw,rwc,rho):
    dm=(4**4*rwc/(np.pi*1e6*Nw))**0.25*1e3 # in mm
    return dm

mu=0.0
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
ijL=[]
kL=[]
for i1,i2 in zip(a[0],a[1]):
    ijL.append([i1,i2])
    kL.append(K[i1,i2])

from bisectm import *

kL=np.array(kL)
#ind=np.argsort(np.array(kL))
#kL=kL[ind]
ijL=np.array(ijL)
#kL=kL.cumsum()/kL.sum()

k0=K.sum()
t=0
it=0
iint=0
S=np.zeros((40,40),float)

Nd0=Nd.copy()

from numba import jit
@jit(nopython=True)
def calcprob(S,K):
    for i in range(40):
        for j in range(40):
            S[i,j]=Nd[i]*Nd[j]*K[i,j]
            
while it<300*3600:
    calcprob(S,K)
    Srow=S.sum(axis=0)
    Srow=Srow.cumsum()
    Srow_sum=Srow[-1]
    Srow/=Srow_sum

    r1=np.random.random()
    r2=np.random.random()
    dt=-np.log(r1)/Srow_sum
    ni1,ni2=bisect2(Srow,40,r2)
    dr=r2-Srow[ni1]
    v=S[ni2,:].cumsum()/Srow_sum
    nj1,nj2=bisect2(v,40,dr)
    i1=ni2
    i2=nj2
    if Nd[i1]>1 and Nd[i2]>1:
        Nd[i1]-=1
        Nd[i2]-=1
        dm_new=(Dm[i1]**3+Dm[i2]**3)**(1/3.0)
        if dm_new<Dm_b[40]:
            n1,n2=bisect2(Dm_b,41,dm_new)
        else:
            n1=40
            #print((dm_new)**3/Dm[n1]**3)
        #stop
        Nd[n1]+=(dm_new)**3/Dm[n1]**3
        #print(n1,n2)
    #    iint+=1
            #print(ij0)
    t+=dt
    it+=1
