from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import combAlg as sdsu
#from bisectm import *


sdsu.mainfortpy()
sdsu.initp2()


import glob
import numpy as np

fs=sorted(glob.glob("wrfout*.aveg.nc"))
nt=0
#zKuL=[]
R=287
h=(125/2.+np.arange(152)*125)/1e3
h1=(0.+np.arange(153)*125)/1e3
Rd=287
cfadKa=np.zeros((90,50),float)
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_filter


piamL=[]
zmL=[]
sfcRainL=[]
import xarray as xr
for f in fs:
    fh=Dataset(f)
    qr=fh["QRAIN"][:]
    nx,ny,nz=qr.shape
    zka3d=np.zeros((nx,ny,152),float)
    zku3d=np.zeros((nx,ny,152),float)
    zkat3d=np.zeros((nx,ny,152),float)
    zkut3d=np.zeros((nx,ny,152),float)
    rwc3d=np.zeros((nx,ny,152),float)
    swc3d=np.zeros((nx,ny,152),float)
    w3d=np.zeros((nx,ny,152),float)
    t3d=np.zeros((nx,ny,152),float)
    prate2d=np.zeros((nx,ny),float)
    piaka_2d=np.zeros((nx,ny),float)
    piaku_2d=np.zeros((nx,ny),float)
    nx,ny,nz=qr.shape
    dn2d=np.random.randn(prate2d.shape[0],prate2d.shape[1])
    from scipy.ndimage import gaussian_filter
    dn2d=gaussian_filter(dn2d,sigma=3)*4
    qv=fh["QVAPOR"][:]
    w=fh["W"][:]
    qh=fh["QICE"][:]
    qs=fh["QSNOW"][:]*0.85
    qg=fh["QGRAUP"][:]*0.85
    ph=fh["PH"][:]+fh["PHB"][:]
    zh=ph/9.81e3
    press=fh["P"][:]+fh["PB"][:]
    dbz=fh["REFL_10CM"][:]
    TH=fh["T"][:]+300
    T=TH*(press/1e5)**(0.287)
    rho=press/T/Rd
    rwc=qr*rho*1e3
    swc=(qs)*rho*1e3
    gwc=(qg)*rho*1e3
    hwc=(qh)*rho*1e3
    a=np.nonzero(dbz[:,:,0]>-45)
    rwcL=[]
    swcL=[]
    gwcL=[]
    hwcL=[]
    zKaL=[]
    for i1,i2 in zip(a[0],a[1]):
        rwc1=np.interp(h,zh[i1,i2,:],rwc[i1,i2,:])
        swc1=np.interp(h,zh[i1,i2,:],swc[i1,i2,:])
        gwc1=np.interp(h,zh[i1,i2,:],gwc[i1,i2,:])
        hwc1=np.interp(h,zh[i1,i2,:],hwc[i1,i2,:])
        w1=np.interp(h,zh[i1,i2,:],w[i1,i2,:])
        temp=np.interp(h,zh[i1,i2,:],T[i1,i2,:])
        temp1=np.interp(h1,zh[i1,i2,:],T[i1,i2,:])
        press1=np.interp(h,zh[i1,i2,:],press[i1,i2,:])
        wv1=np.interp(h,zh[i1,i2,:],qv[i1,i2,:])
        sl=2.5
        dr=0.125
        rwcL.append(rwc1)
        swcL.append(swc1)
        gwcL.append(gwc1)
        hwcL.append(hwc1)
        dn1=np.zeros((152),float)-0.5+dn2d[i1,i2]
        zka_m ,zka_t, attka, piaka, \
            kext_ka,salb_ka,asym_ka,kext_ka_,salb_ka_,asym_ka_,pRateKa\
            =sdsu.reflectivity_ka(rwc1,swc1+hwc1+gwc1,wv1,dn1,temp,press1,dr,sl)
        zku_m ,zku_t, attku, piaku, \
            kext_ku,salb_ku,asym_ku,kext_ku_,salb_ku_,asym_ku_,pRateKu\
            =sdsu.reflectivity_ku2(rwc1,swc1+hwc1+gwc1,wv1,dn1,temp,press1,dr,sl)
        zKaL.append(zka_m)
        rwc3d[i1,i2,:]=rwc1
        swc3d[i1,i2,:]=swc1+hwc1+gwc1
        prate2d[i1,i2]=pRateKa[0]
        piaka_2d[i1,i2]=piaka
        piaku_2d[i1,i2]=piaku
        w3d[i1,i2,:]=w1
        t3d[i1,i2,:]=temp
        dr=0.125
        noms=0
        alt=400.
        freqKa=35.5
        nonorm=0
        theta=0.35/1.0
       
        if zka_m.max()>25:
            zms = sdsu.multiscatterf(kext_ka[::-1],salb_ka[::-1],asym_ka[::-1],\
                                     zka_t[::-1],dr,noms,alt,\
                                     theta,freqKa,nonorm)
            
            zka_m=zms[::-1]
        zka3d[i1,i2,:]=zka_m
        zku3d[i1,i2,:]=zku_m
        zkat3d[i1,i2,:]=zka_t
        zkut3d[i1,i2,:]=zku_t
        for k in range(-90):
            idbz=int(zka_m[k])
            if idbz>=0 and idbz<50 and k<90 and zka_m[k]>4:
                cfadKa[k,idbz]+=1
    pwc1=np.array(swcL).mean(axis=0)+np.array(gwcL).mean(axis=0)\
        +np.array(hwcL).mean(axis=0)
    pwc2=np.array(rwcL).mean(axis=0)
    plt.subplot(121)
    plt.plot(np.array(rwcL).mean(axis=0),h)
    plt.plot(np.array(swcL).mean(axis=0),h)
    plt.plot(np.array(gwcL).mean(axis=0),h)
    plt.plot(np.array(hwcL).mean(axis=0),h)
    plt.plot(pwc1,h)
    plt.legend(['r','s','g','h'])
    temp=np.interp(h,zh[i1,i2,:],T[i1,i2,:])
    temp1=np.interp(h1,zh[i1,i2,:],T[i1,i2,:])
    press1=np.interp(h,zh[i1,i2,:],press[i1,i2,:])
    wv1=np.interp(h,zh[i1,i2,:],qv[i1,i2,:])
    #sl=1.5
    dr=0.125
    dn1=np.zeros((152),float)
    #zka_m ,zka_t, attka, piaka, \
    #    kext_ka,salb_ka,asym_ka,kext_ka_,salb_ka_,asym_ka_,pRateKa\
    #    =sdsu.reflectivity_ku2(pwc2,pwc1,wv1,dn1,temp,press1,dr,sl)
    #plt.subplot(122)
    #plt.plot(zka_m,h)
    #plt.plot(np.array(zKaL).mean(axis=0),h)
    #plt.plot(zka_t,h)
    #plt.ylim(0,15)

    a=np.nonzero(dbz[:,:,0]>15)
    di=[i for i in range(-1,2) for j in range(-1,2)]
    dj=[j for i in range(-1,2) for j in range(-1,2)]

    w1=np.exp(-1.0*np.array(di)**2)
    w2=np.exp(-1.0*np.array(di)**2)
    for i1,i2 in zip(a[0],a[1]):
        if i1>0 and i2>0 and i1<nx-1 and i2<ny-1:
            zm=np.zeros((152),float)
            piam=0
            sfcrain_m=0
            wt=0
            for di1, di2, w11,w12 in zip(di,dj,w1,w2):
                zm+=10**(0.1*zka3d[i1+di1,i2+di2,:])*w11*w12
                piam+=10**(-0.1*piaka_2d[i1+di1,i2+di2])*w11*w12
                sfcrain_m+=prate2d[i1+di1,i2+di2]*w11*w12
                wt+=w11*w12
            zm=np.log10(zm/wt)*10
            piam=-10*np.log10(piam/wt)
            sfcrain_m/=wt
            if sfcrain_m>20 or sfcrain_m<1.:
                continue
            zmL.append(zm[0:94:2][::-1])
            piamL.append(piam)
            sfcRainL.append(sfcrain_m)
            for k in range(90):
                idbz=int(zm[k])
                if idbz>=0 and idbz<50 and k<90 and zm[k]>4:
                    cfadKa[k,idbz]+=1
    plt.figure()
    plt.pcolormesh(cfadKa,norm=matplotlib.colors.LogNorm(),cmap='jet')
    plt.xlim(0,40)
    
    import pickle
    pickle.dump({"zmKa":np.array(zmL),"piaKa":np.array(piamL),\
                 "sfcRain":np.array(sfcRainL)},open("11june2014_05.pklz","wb"))

    zka3d_xr=xr.DataArray(zka3d)
    zku3d_xr=xr.DataArray(zku3d)
    zkat3d_xr=xr.DataArray(zkat3d)
    zkut3d_xr=xr.DataArray(zkut3d)
    piaka_2d_xr=xr.DataArray(piaka_2d)
    piaku_2d_xr=xr.DataArray(piaku_2d)
    prate2d_xr=xr.DataArray(prate2d)
    w3d_xr=xr.DataArray(w3d)
    rwc3d_xr=xr.DataArray(rwc3d)
    swc3d_xr=xr.DataArray(swc3d)
    temp3d_xr=xr.DataArray(t3d)
    dn2d_xr=xr.DataArray(dn2d)
    ds=xr.Dataset({"zka":zka3d_xr,"zku":zku3d_xr,\
                   "zkat":zkat3d_xr,"zkut":zkut3d_xr,\
                   "piaKa":piaka_2d_xr,"prate2d":prate2d_xr,\
                   "w3d":w3d_xr,"rwc3d":rwc3d_xr,"swc3d":swc3d_xr,"temp3d":temp3d_xr,\
                   "piaKu":piaku_2d_xr,"dn2d":dn2d_xr})
    fout=f.replace("aveg","aveg_grid")
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(fout, encoding=encoding)

    stop
    for i1 in range(nx):
        for i2 in range(ny):
            if dbz[i1,i2,:].max()<10:
                continue
            rwc1=np.interp(h,zh[i1,i2,:],rwc[i1,i2,:])
            swc1=np.interp(h,zh[i1,i2,:],swc[i1,i2,:])
            gwc1=np.interp(h,zh[i1,i2,:],gwc[i1,i2,:])
            hwc1=np.interp(h,zh[i1,i2,:],hwc[i1,i2,:])
            temp=np.interp(h,zh[i1,i2,:],T[i1,i2,:])
            temp1=np.interp(h1,zh[i1,i2,:],T[i1,i2,:])
            press1=np.interp(h,zh[i1,i2,:],press[i1,i2,:])
            wv1=np.interp(h,zh[i1,i2,:],qv[i1,i2,:])
            sl=1
            dr=0.125
            dn1=np.zeros((152),float)
            zku_m ,zku_t, attku, piaku, \
                kext,salb,asym,kext_,salb_,asym_,pRate\
                =sdsu.reflectivity_ku(rwc1,swc1,gwc1,hwc1,wv1,dn1,temp,press1,dr,sl)
            zka_m ,zka_t, attka, piaka, \
                kext_ka,salb_ka,asym_ka,kext_ka_,salb_ka_,asym_ka_,pRateKa\
                =sdsu.reflectivity_ka(rwc1,swc1+gwc1+hwc1,wv1,dn1,temp,press1,dr,sl)
            a=np.nonzero(temp<273.15)
            hz=h[a[0][0]]
            iz0=int(h[a[0][0]]/dr)
            di0=30-iz0
            di0=0
            for k in range(90):
                if k+di0>=0:
                    idbz=int(zka_t[k])
                    if idbz>=0 and idbz<50 and k+di0<90 and zku_m[k]>4:
                        cfadKa[k+di0,idbz]+=1
