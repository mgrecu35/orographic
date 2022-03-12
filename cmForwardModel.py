import combAlg as cmb
cmb.mainfortpy()
cmb.initp2()
from netCDF4 import Dataset


from scipy.special import gamma as gam
import numpy as np


def nw_lambd(swc,nc,mu):
    rhow=1e6
    lambd=(nc*rhow*np.pi*gam(4+mu)/gam(1+mu)/6.0/swc)**(0.333)  # m-1
    n0=nc*lambd/gam(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    return n0,lambd

from numba import jit

@jit(nopython=False)
def get_Zn(w,nw,lambd,W,Z,att,dm,dm_out,Deq,bscat,ext,vfall,mu,wl):
    dD=0.05
    rhow=1 #gcm-3
    Dint=np.arange(160)*dD+dD/2.0
    bscatInt=np.interp(Dint,Deq,bscat)
    extInt=np.exp(np.interp(Dint,Deq,np.log(ext)))  #m^2
    vfallInt=np.interp(Dint,Deq,vfall)
    fact=1e3/np.pi**5/0.93*wl**4
   
    nP=W.shape[0]
    f_mu=6/4**4*(4+mu)**(mu+4)/gam(mu+4)
    
    for j in range(nP):
        vdop=0
        nc0=0
        Vol=0
        zray=0.0
        for i in range(160):
            d=dD*i+dD/2
            Nd=f_mu*np.exp(-lambd[j]*d)*(d/dm[j])**mu*dD #(mm)
            W[j]=W[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow #(g/m3)
            dm_out[j]=dm_out[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow*(d) #(g/m3)
            Z[j]=Z[j]+nw[j]*Nd*bscatInt[i]
            vdop=vdop+nw[j]*Nd*bscatInt[i]*vfallInt[i]
            att[j]=att[j]+nw[j]*Nd*extInt[i]*4.343 #(/km)1
            nc0=nc0+nw[j]*Nd
            Vol=Vol+nw[j]*Nd*(1e-3*d)**3*np.pi/6
        Z[j]=np.log10(Z[j]*fact)*10
        dm_out[j]=dm_out[j]/(W[j]+1e-9)

        
fnameIce='/home/grecu/scatter-1.1/ice-self-similar-aggregates_13-GHz_scat.nc'
fnameRain='/home/grecu/scatter-1.1/liquid-water_13-GHz_scat.nc'

fnameIce35='/home/grecu/scatter-1.1/ice-self-similar-aggregates_35-GHz_scat.nc'
fnameRain35='/home/grecu/scatter-1.1/liquid-water_35-GHz_scat.nc'

def readScatProf(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    fraction=fh['fraction'][:]
    bscat=fh['bscat'][:]*4*np.pi
    Deq=10*(mass*1e3*6/np.pi)**(0.333) # in mm
    ext=fh['ext'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    vfall=fh['fall_speed'][:]
    return temp,mass,fraction,bscat,Deq,ext,scat,g,vfall

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
    return temp,mass,bscat,Deq,ext,scat,g,vfall,fh

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,fh=readScatProfR(fnameRain)
wl=fh['wavelength'][:]*1000
fact=1e6/np.pi**5/0.93*wl**4
import matplotlib.pyplot as plt
#plt.semilogy(bscat_r[9,:]*fact)
#plt.semilogy(Deq_r**6)
#stop
tempKa,massKa,fractionKa,bscatKa,DeqKa,extKa,scatKa,gKa,vfallKa=readScatProf(fnameIce35)
tempKa_r,massKa_r,bscatKa_r,DeqKa_r,extKa_r,scatKa_r,gKa_r,vfallKa_r,fh=readScatProfR(fnameRain35)

freq=13.8
freqKa=35.5

wl=300/freq
wlKa=300/freqKa
mu=2
@jit(nopython=True)
def gett_atten(z_att_m,z_m,att_tot,z):
    a=np.nonzero(z_m[0,:,:]>0)
    nz=z_att_m.shape[0]
    for i, j in zip(a[0],a[1]):
        pia_tot=0
        for k in range(nz-1,-1,-1):
            if z_m[k,i,j]>-10:
                pia_tot+=att_tot[k,i,j]*(z[k+1,i,j]-z[k,i,j])*4.343
                z_att_m[k,i,j]-=pia_tot
                pia_tot+=att_tot[k,i,j]*(z[k+1,i,j]-z[k,i,j])*4.343
            else:
                z_att_m[k,i,j]-=pia_tot

                
def calcZ(rwc,swc,gwc,ncr,ncs,ncg,z,Deq,ext,bscat,scat,g,vfall,\
          Deq_r,ext_r,bscat_r,scat_r,g_r,vfall_r,wl):
    print(wl)
    att_total=rwc.copy()*0
    z_total=rwc.copy()*0.0

    a=np.nonzero(rwc>0.01)
    ncr[ncr<0.001]=0.001
    #ncr[a]=10
    nw_r,lambd_r=nw_lambd(rwc[a],ncr[a],mu)
    nwr=rwc.copy()*0.0
    w_r=rwc[a].copy()*0.0
    z_r=rwc[a].copy()*0.0
    att_r=rwc[a].copy()*0.0
    dm_r=rwc[a].copy()*0.0
    print(len(a[0]))
    get_Z(rwc[a],nw_r,lambd_r,w_r,z_r,att_r,dm_r,Deq_r,bscat_r[9,:],ext_r[9,:],vfall_r,mu,wl)
    nwr[a]=np.log10(nw_r)
    
    z_total[a]+=10.**(0.1*z_r)
    att_total[a]+=att_r

    a=np.nonzero(swc>0.01)
    ncs[ncs<0.001]=0.001
    nw_s,lambd_s=nw_lambd(swc[a],ncs[a],mu)
    nws=rwc.copy()*0.0
    
    w_s=swc[a].copy()*0.0
    z_s=swc[a].copy()*0.0
    att_s=swc[a].copy()*0.0
    dm_s=swc[a].copy()*0.0
    get_Z(swc[a],nw_s,lambd_s,w_s,z_s,att_s,dm_s,Deq[12,:],bscat[-1,12,:],ext[-1,12,:],\
          vfall[12,:],mu,wl)
    nws[a]=np.log10(nw_s)
    z_total[a]+=10.**(0.1*z_s)
    att_total[a]+=att_s
    
    a=np.nonzero(gwc>0.01)
    nw_g,lambd_g=nw_lambd(gwc[a],ncg[a],mu)
    nwg=rwc.copy()*0.0
    
    w_g=gwc[a].copy()*0.0
    z_g=gwc[a].copy()*0.0
    att_g=gwc[a].copy()*0.0
    dm_g=gwc[a].copy()*0.0
    get_Z(gwc[a],nw_g,lambd_g,w_g,z_g,att_g,dm_g,Deq[14,:],bscat[-1,14,:],ext[-1,14,:],\
          vfall[14,:],mu,wl)
    nwg[a]=np.log10(nw_g)
        
    z_total[a]+=10.**(0.1*z_g)
    att_total[a]+=att_g
    
    z_total=10*np.log10(z_total+1e-9)
    z_m=np.ma.array(z_total,mask=z_total<-10)
    z_att_m=z_m.copy()
    gett_atten(z_att_m,z_m,att_total,z)
    return z_m,z_att_m, att_total, nwr,nws,nwg, w_r, nw_r, z_r
    
import matplotlib.pyplot as plt
#plt.hist(np.log10(nw_s/0.08))

#z_m,z_att_m=calcZ(rwc,swc,gwc,ncr,ncs,ncg,z,Deq,ext,scat,g,vfall,\
#                  Deq_r,ext_r,scat_r,g_r,vfall_r,wl)

mu=2.
rwc=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1,1.5,2.,2.5,3.])
nw_dm=np.loadtxt("NwDm.txt")

def calcZkuR(rwc,Deq_r,bscat_r,ext_r,vfall_r,mu,wl,nw_dm):
    Nw=rwc.copy()*0+0.08e8
    dm=(4**4*rwc/(np.pi*1e6*Nw))**0.25*1e3 # in mm
    for i in range(dm.shape[0]):
        inw=int((dm[i]-0.5)/0.04)
        if inw>59:
            inw=59
        dnw=0.5*(nw_dm[inw,1]-6.9)+0.0125*np.random.randn()
        Nw[i]=10**(dnw+6.9)
    nwr=rwc.copy()*0.0+Nw
    w_r=rwc.copy()*0.0
    z_r=rwc.copy()*0.0
    att_r=rwc.copy()*0.0
    dm_r=rwc.copy()*0.0
    w_s=rwc.copy()*0.0
    z_s=rwc.copy()*0.0
    att_s=rwc.copy()*0.0
    dm_s=rwc.copy()*0.0
    lambd=(4+mu)/dm
    get_Zn(rwc,nwr,lambd,w_r,z_r,att_r,dm,dm_r,Deq_r,bscat_r[9,:],ext_r[9,:],vfall_r,mu,wl)
    return nwr,z_r,att_r

def calcZkuS(rwc,T,Deq,bscat,ext,vfall,mu,wl):
    Nw=rwc.copy()*0+0.08e8*8
    a=np.nonzero(T>-10)
    b=np.nonzero(T[a]<0)
    Nw[a][b]=0.08e8*(8-7*(T[a][b]+10)/10.)
    a=np.nonzero(T>0)
    Nw[a]=0.08e8
    dm=(4**4*rwc/(np.pi*1e6*Nw))**0.25*1e3 # in mm
    nwr=rwc.copy()*0.0+Nw
    w_r=rwc.copy()*0.0
    z_r=rwc.copy()*0.0
    att_r=rwc.copy()*0.0
    dm_r=rwc.copy()*0.0
    w_s=rwc.copy()*0.0
    z_s=rwc.copy()*0.0
    att_s=rwc.copy()*0.0
    dm_s=rwc.copy()*0.0
    lambd=(4+mu)/dm
    get_Zn(rwc,nwr,lambd,w_s,z_s,att_s,dm,dm_s,Deq[14,:],bscat[3,14,:],ext[3,14,:],vfall[3,:],mu,wl)
    return nwr,z_s,att_s

fname="../../Data/wrfout_d03_2018-06-25_02:00:00"
fname="cm1out_000008.nc"
import numpy as np
def read_wrf(fname,it):
    f=Dataset(fname)
    qv=f['qv'][it,:,:,:]    # water vapor
    qr=f['qr'][it,:,:,:]     # rain mixing ratio
    qs=f['qs'][it,:,:,:]     # snow mixing ratio
    qc=f['qc'][it,:,:,:]    # cloud mixing ratio
    qg=f['qg'][it,:,:,:]  # graupel mixing ratio
    ncr=f['ncr'][it,:,:,:]     # rain mixing ratio
    ncs=f['ncs'][it,:,:,:]     # snow mixing ratio
    ncg=f['ncg'][it,:,:,:]   # graupel mixing ratio
    #z=f['z_coords'][:]/1000.             # height (km)
    th=f['th'][it,:,:,:]+300    # potential temperature (K)
    prs=f['prs'][it,:,:,:]
    # pressure (Pa)
    T=th*(prs/100000)**0.286  # Temperature
    t2c=T-273.15
    z=np.arange(81)*.25
    #stop
  
    R=287.058  #J*kg-1*K-1
    rho=prs/(R*T)
    return qr,qs,qg,ncr,ncs,ncg,rho,z,T,f
it=-1
qr,qs,qg,ncr,ncs,ncg,rho,z,T,fh=read_wrf(fname,it)


swc=rho*(qs)*1.05e3
gwc=rho*qg*1.05e3
rwc=rho*qr*1.05e3
ncr=ncr*rho*1
ncg=ncg*rho*1
ncs=(ncs)*rho*1
zt=swc.copy()*0+1e-9
a=np.nonzero(rwc>0.001)
nwr_a,z_r_a,att_r_a=calcZkuR(rwc[a],Deq_r,bscat_r,ext_r,vfall_r,mu,wl,nw_dm)
zt[a]+=10**(0.1*z_r_a)
a=np.nonzero(swc>0.001)
nws_a,z_s_a,att_s_a=calcZkuS(swc[a],T[a],Deq,bscat,ext,vfall,mu,wl)
a=np.nonzero(gwc>0.001)
nwg_a,z_g_a,att_g_a=calcZkuS(gwc[a],T[a],Deq,bscat,ext,vfall,mu,wl)

a=np.nonzero(swc>0.001)
ncs[ncs<0.001]=0.001
nw_s,lambd_s=nw_lambd(swc[a],ncs[a],mu)
nws=swc.copy()*0.0
w_s=swc[a].copy()*0.0
z_s=swc[a].copy()*0.0
att_s=swc[a].copy()*0.0
dm_s=swc[a].copy()*0.0
@jit(nopython=False)
def get_Z(w,nw,lambd,W,Z,att,dm,Deq,bscat,ext,vfall,mu,wl):
    dD=0.05
    rhow=1 #gcm-3
    Dint=np.arange(160)*dD+dD/2.0
    bscatInt=np.interp(Dint,Deq,bscat)
    extInt=np.exp(np.interp(Dint,Deq,np.log(ext)))  #m^2
    vfallInt=np.interp(Dint,Deq,vfall)
    fact=1e6/np.pi**5/0.93*wl**4
    print(W.shape)
    nP=W.shape[0]
    #print(nP,mu,fact)
    print(fact)
    #print(bscatInt)
    for j in range(nP):
        vdop=0
        nc0=0
        Vol=0
        zray=0.0
        for i in range(160):
            d=dD*i+dD/2
            Nd=np.exp(-lambd[j]*d*0.1)*(0.1*lambd[j]*d)**mu*dD #(mm)
            W[j]=W[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow #(g/m3)
            dm[j]=dm[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow*(0.1*d) #(g/m3)
            Z[j]=Z[j]+nw[j]*Nd*bscatInt[i]
            vdop=vdop+nw[j]*Nd*bscatInt[i]*vfallInt[i]
            att[j]=att[j]+nw[j]*Nd*extInt[i]*1e3 #(/km)1
            nc0=nc0+nw[j]*Nd
            Vol=Vol+nw[j]*Nd*(1e-3*d)**3*np.pi/6
        Z[j]=np.log10(Z[j]*fact)*10
        dm[j]=dm[j]/(W[j]+1e-9)
        nw_d=4.0**4/np.pi/1e1*W[j]/dm[j]**4
        nw[j]=nw_d
get_Z(swc[a],nw_s,lambd_s,w_s,z_s,att_s,dm_s,Deq[24,:],bscat[-1,24,:],ext[-1,24,:],\
      vfall[24,:],mu,wl)

a=np.nonzero(swc>0.001)
zt[a]+=10**(0.1*z_s_a)
a=np.nonzero(gwc>0.001)
zt[a]+=10**(0.1*z_g_a)
zt=np.log10(zt)*10.

ztm=np.ma.array(zt,mask=zt<0)
stop
for i in range(150,190,5):
    plt.figure()
    plt.pcolormesh(range(zt.shape[-1]),z[:-1,0,0],ztm[:,i,:],vmin=0,cmap='jet',\
                   vmax=50)
    plt.ylim(1,15)
    plt.colorbar()

cfad=np.zeros((55,60),float)
hgrid=1+np.arange(60)*0.25

@jit(nopython=False)
def makecfad(z_m,z,cfad,hgrid):
    a=np.nonzero(z_m[0,:,:]>30)
    for i, j in zip(a[0],a[1]):
        zm1=np.interp(hgrid,0.5*(z[:-1,i,j]+z[1:,i,j]),z_m[:,i,j])
        for k in range(60):
            i0=int(zm1[k])
            if i0>=0 and i0<50:
                cfad[i0,k]+=1

makecfad(zt,z,cfad,hgrid)
