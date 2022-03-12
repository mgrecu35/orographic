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
n_neighbors=30
from sklearn.neighbors import KNeighborsRegressor
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
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
    zku=fh["zku"][:]
    zkut=fh["zkut"][:]
    temp=fh["temp3d"][:]
    a=np.nonzero(zku[:,:,25]>0)
    b=np.nonzero(a[0]==105)
    piaKu=fh["piaKu"][:]
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
        if zkut[i1,i2,33]-zku[i1,i2,33]>0 and zkut[i1,i2,33]>10 and piaKu[i1,i2]<40:
            attZL.append((10**(0.1*zku[i1,i2,34:]*0.72)).sum())
            attL.append(zkut[i1,i2,33]-zku[i1,i2,33])
            pia2dL.append(piaKu[i1,i2])
            zRL.append(zku[i1,i2,:38])
            z1=zku[i1,i2,:38].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(((z1-zm)/zs)[:32])
            z1l.append((pia2dL[-1]+np.random.randn()*3-6)/3.0)
            xtfL_1.append(z1l[:32][::-1])
            xtfL_2.append(z1l[-1])
            z1L.append(z1l)
            pRateL.append(prate2d[i1,i2])
        else:
            ic40+=1
            attZ40L.append((10**(0.1*zku[i1,i2,34:]*0.72)).sum())
            att40L.append(zkut[i1,i2,33]-zku[i1,i2,33])
            pia2d40L.append(piaKu[i1,i2])
            z1=zku[i1,i2,:38].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(((z1-zm)/zs)[:32])
            z1l.append(5)
            z140L.append(z1l)
            xtfL_40_1.append(z1l[:32][::-1])
            xtfL_40_2.append(z1l[-1])
            pRate40L.append(prate2d[i1,i2])
    return z1L,pRateL,z140L,pRate40L,xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2

import glob
fnames=glob.glob("../orographic2/*.aveg_grid.nc")
for fname in fnames:
    z1L,pRateL,z140L,pRate40L,\
        xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2=readData(fname,z1L,pRateL,z140L,pRate40L,\
                                                   xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2)
    
x_train,x_test,\
    y_train, y_test = train_test_split(np.array(z1L)[:,:], np.array(pRateL),\
                                       test_size=0.5, random_state=42)

knn.fit(x_train,y_train)
yp=knn.predict(x_test)

for intv in [[1,2],[2,4],[4,6],[6,8],[8,10],[10,20],[20,30]]:
    a=np.nonzero((y_test-intv[0])*(y_test-intv[1])<0)
    rms=(((yp[a[0]]-y_test[a[0]])**2).mean())**0.5/y_test[a[0]].mean()
    print("%6.2f %6.2f %6.2f %6.2f"%(intv[0],intv[1],\
                                     (-1+yp[a[0]].mean()/y_test[a[0]].mean())*100,rms*100))


x_train40,x_test40,\
    y_train40, y_test40 = train_test_split(np.array(z140L), np.array(pRate40L),\
                                       test_size=0.5, random_state=42)

knn.fit(x_train40,y_train40)
yp40=knn.predict(x_test40)
