fname="../orographic2/wrfout_d04_2014-06-11_18:36:00.aveg_grid.nc"

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import combAlg as sdsu
#from bisectm import *
import numpy as np

#sdsu.mainfortpy()
#sdsu.initp2()

from sklearn import neighbors
n_neighbors=50
from sklearn.neighbors import KNeighborsRegressor
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
knnC = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
from sklearn.model_selection import train_test_split

z1L=[]
pRateL=[]
z140L=[]
pRate40L=[]
xtfL_1=[]
xtfL_2=[]
xtfL_40_1=[]
xtfL_40_2=[]
    
def readData(fname,z1L,pRateL,z140L,pRate40L,\
             xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2):
    fh=Dataset(fname)
    zka=fh["zka"][:]
    zkat=fh["zkat"][:]
    temp=fh["temp3d"][:]
    a=np.nonzero(zka[:,:,25]>0)
    b=np.nonzero(a[0]==105)
    piaKa=fh["piaKa"][:]
    prate2d=fh["prate2d"][:]
    zRL=[]
    z40RL=[]
    attL=[]
    attZL=[]
    pia2dL=[]
    att40L=[]
    attZ40L=[]
    pia2d40L=[]
    
    zs=np.array([8.21868427, 8.27254699, 8.30685336, 8.30634989, 8.26425975,
                 8.23223687, 8.13703432, 8.08034568, 7.97254592, 7.86304235,
                 7.73122532, 7.59266912, 7.46987865, 7.30355459, 7.18684333,
                 7.00568015, 6.82390663, 6.63781607, 6.42879146, 6.25028939,
                 6.08009431, 6.13512346, 6.17414461, 6.226422  , 6.26672013,
                 6.29542539, 6.32231197, 6.33430186, 6.3356904 , 6.36359843,
                 6.21911478, 6.00297172, 5.82091536, 5.,5.,5.,5.,5.])
    zm=np.array([17.88353431, 18.05536806, 18.21657803, 18.39863995, 18.61310613,
                 18.84283269, 19.12744192, 19.39869986, 19.72046197, 20.05707027,
                 20.41630877, 20.78765812, 21.15194279, 21.53537528, 21.8916305 ,
                 22.27120987, 22.64107185, 23.00349155, 23.36079428, 23.69367785,
                 24.0106508 , 24.23204546, 24.43642963, 24.60436233, 24.7196606 ,
                 24.75354306, 24.64879862, 24.35431599, 23.81914873, 23.07190774,
                 22.48808443, 22.33656934, 22.31730397, 22.,22.,22.,22.,22.])
    
    ic40=0
    for i1,i2 in zip(a[0],a[1]):
        if piaKa[i1,i2]<40:
            attZL.append((10**(0.1*zka[i1,i2,34:]*0.72)).sum())
            attL.append(zkat[i1,i2,33]-zka[i1,i2,33])
            pia2dL.append(piaKa[i1,i2])
            zRL.append(zka[i1,i2,:38])
            z1=zka[i1,i2,:38].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(((z1-zm)/zs)[:32])
            z1l.append((pia2dL[-1]+np.random.randn()*2-6)/8.0)
            xtfL_1.append(z1l[:32][::-1])
            xtfL_2.append(z1l[-1])
            z1L.append(z1l)
            pRateL.append(prate2d[i1,i2])
        if piaKa[i1,i2]>40:
            ic40+=1
            attZ40L.append((10**(0.1*zka[i1,i2,34:]*0.72)).sum())
            att40L.append(zkat[i1,i2,33]-zka[i1,i2,33])
            pia2d40L.append(piaKa[i1,i2])
            z1=zka[i1,i2,23:75].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(z1[:])
            #z1l.append()
            z140L.append(z1l)
            xtfL_40_1.append(z1l[:32][::-1])
            xtfL_40_2.append(z1l[-1])
            pRate40L.append(prate2d[i1,i2])
    return z1L,pRateL,z140L,pRate40L,xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2

import glob
fnames=glob.glob("../orographic2/*.aveg_grid.nc")
def readPrecipData(fname,knn,knnC):
    fh=Dataset(fname)
    prate2d=fh["prate2d"][:]
    zka=fh["zka"][:]
    piaKa=fh["piaKa"][:]
    nx,ny,nz=zka.shape
    zs=np.array([8.21868427, 8.27254699, 8.30685336, 8.30634989, 8.26425975,
                 8.23223687, 8.13703432, 8.08034568, 7.97254592, 7.86304235,
                 7.73122532, 7.59266912, 7.46987865, 7.30355459, 7.18684333,
                 7.00568015, 6.82390663, 6.63781607, 6.42879146, 6.25028939,
                 6.08009431, 6.13512346, 6.17414461, 6.226422  , 6.26672013,
                 6.29542539, 6.32231197, 6.33430186, 6.3356904 , 6.36359843,
                 6.21911478, 6.00297172, 5.82091536, 5.,5.,5.,5.,5.])
    zm=np.array([17.88353431, 18.05536806, 18.21657803, 18.39863995, 18.61310613,
                 18.84283269, 19.12744192, 19.39869986, 19.72046197, 20.05707027,
                 20.41630877, 20.78765812, 21.15194279, 21.53537528, 21.8916305 ,
                 22.27120987, 22.64107185, 23.00349155, 23.36079428, 23.69367785,
                 24.0106508 , 24.23204546, 24.43642963, 24.60436233, 24.7196606 ,
                 24.75354306, 24.64879862, 24.35431599, 23.81914873, 23.07190774,
                 22.48808443, 22.33656934, 22.31730397, 22.,22.,22.,22.,22.])
    pRate2D_Ret=prate2d.copy()*0
    z1L=[]
    z1CL=[]
    ijL=[]
    ijCL=[]
    p1=[]
    p2=[]
    zKaL=[]
    for i1 in range(nx):
        for i2 in range(ny):
            if zka[i1,i2,:32].max()>0:
                z1=zka[i1,i2,:38].copy()
                z1[z1<0]=0
                z1+=np.random.randn()
                z1l=list(((z1-zm)/zs)[:32])
                z1l.append((piaKa[i1,i2]+np.random.randn()*2-6)/8.0)
                if piaKa[i1,i2]<40:
                    zKaL.append(zka[i1,i2,:])
                    ijL.append([i1,i2])
                    z1L.append(z1l)
                    p1.append(prate2d[i1,i2])
                    #pRet=knn.predict([z1l])
                    #pRate2D_Ret[i1,i2]=pRet[0]
                else:
                    #z1l[-1]=5
                    z1l=zka[i1,i2,23:75].copy()
                    p2.append(prate2d[i1,i2])
                    #pRetC=knn.predict([z1l])
                    #pRate2D_Ret[i1,i2]=pRetC[0]
                    ijCL.append([i1,i2])
                    z1CL.append(z1l)
    pRate=knn.predict(np.array(z1L))
    print(np.mean(p1),pRate.mean())
    pRateC=knnC.predict(np.array(z1CL))
    print(np.mean(p2),pRateC.mean())
    for i,p1 in enumerate(pRate):
        pRate2D_Ret[ijL[i][0],ijL[i][1]]=p1
    for i,p1 in enumerate(pRateC):
        pRate2D_Ret[ijCL[i][0],ijCL[i][1]]=p1
    return prate2d, pRate2D_Ret,zKaL

import pickle
plt.rcParams.update({'font.size': 11})
knnDict=pickle.load(open("knnKa.pklz","rb"))

knn=knnDict["knn"]
knnC=knnDict["knnC"]
fnames=sorted(fnames)

iplot=0
if iplot==1:
    fh=Dataset(fnames[0])
    zka=fh["zka"][:]
    zku=fh["zku"][:]
    piaKa=fh["piaKa"][:]
    zkam=np.ma.array(zka,mask=zka<0)
    zkum=np.ma.array(zku,mask=zku<0)
    plt.figure(figsize=(8,8))
    plt.subplot(211)
    plt.pcolormesh(np.arange(zka.shape[1]),np.arange(zka.shape[2])*0.125,zkum[100,:,:].T,vmin=0,vmax=45,\
                   cmap='jet')
    plt.title('Ku-band')
    plt.ylim(0,15)
    plt.xlim(0,300)
    plt.ylabel('Height (km)')
    cb=plt.colorbar()
    cb.ax.set_title('dBZ')
    plt.subplot(212)
    plt.pcolormesh(np.arange(zka.shape[1]),np.arange(zka.shape[2])*0.125,zkam[100,:,:].T,vmin=0,vmax=40,\
                   cmap='jet')
    plt.title('Ka-band')
    plt.ylim(0,15)
    plt.xlim(0,300)
    plt.ylabel('Height (km)')
    plt.ylabel('Distance (km)')
    cb=plt.colorbar()
    cb.ax.set_title('dBZ')
    plt.savefig('radarObs.11June2014.png')

pRate2D,pRate2D_ret,zKaL=readPrecipData(fnames[0],knn,knnC)

stop

from numba import jit
nx,ny=pRate2D.shape
pRate2D_av=np.zeros((nx-25,ny-25),float)
pRate2D_ret_av=np.zeros((nx-25,ny-25),float)
pIntv=np.array([-1,0.01,1,2,4,6,8,10,20,30,50,100,150])
hist2D=np.zeros((5,13),float)
@jit(nopython=False)
def distrib(pRate2D,pRate2D_av,hist2D,pIntv):
    nx,ny=pRate2D.shape
    for i in range(nx-25):
        for j in range(ny-25):
            pRate2D_av[i,j]=pRate2D[i:i+25,j:j+25].mean()
            


distrib(pRate2D,pRate2D_av,hist2D,pIntv)
distrib(pRate2D_ret,pRate2D_ret_av,hist2D,pIntv)
a=np.nonzero(pRate2D_av>0.01)
pIntv=np.array([.01,1,2,4,6,8,10])

p1=pRate2D_ret_av[a]
p_ref=pRate2D_av[a]

for intv0,intv1 in zip(pIntv[:-1],pIntv[1:]):
    a=np.nonzero((p_ref-intv0)*(p_ref-intv1)<0)
    rms=(((p1[a[0]]-p_ref[a[0]])**2).mean())**0.5/p_ref[a[0]].mean()
    print("%6.2f %6.2f %6.2f %6.2f"%(intv0,intv1,\
                                     (-1+p1[a[0]].mean()/p_ref[a[0]].mean())*100,rms*100))
    
#stop

for fname in fnames:
    z1L,pRateL,z140L,pRate40L,\
        xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2=readData(fname,z1L,pRateL,z140L,pRate40L,\
                                                   xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2)
    
x_train,x_test,\
    y_train, y_test = train_test_split(np.array(z1L)[:,:], np.array(pRateL),\
                                       test_size=0.5, random_state=42)

knn.fit(x_train,y_train)
yp=knn.predict(x_test)

pIntv=np.array([.01,1,2,4,6,8,10,20,30,50,100,150])
for intv0,intv1 in zip(pIntv[:-1],pIntv[1:]):
    a=np.nonzero((y_test-intv0)*(y_test-intv1)<0)
    rms=(((yp[a[0]]-y_test[a[0]])**2).mean())**0.5/y_test[a[0]].mean()
    print("%6.2f %6.2f %6.2f %6.2f"%(intv0,intv1,\
                                     (-1+yp[a[0]].mean()/y_test[a[0]].mean())*100,rms*100))

#stop
x_train40,x_test40,\
    y_train40, y_test40 = train_test_split(np.array(z140L), np.array(pRate40L),\
                                       test_size=0.5, random_state=42)

knnC.fit(x_train40,y_train40)
yp=knnC.predict(x_test40)


import pickle

pickle.dump({"knn":knn,"knnC":knnC},open("knnKa.pklz","wb"))
