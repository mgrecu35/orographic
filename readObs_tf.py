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
from sklearn.model_selection import train_test_split

def readData(fname):
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
    z1L=[]
    pRateL=[]
    
    att40L=[]
    attZ40L=[]
    pia2d40L=[]
    z140L=[]
    pRate40L=[]
    xtfL_1=[]
    xtfL_2=[]
    xtfL_40_1=[]
    xtfL_40_2=[]
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
        if zkat[i1,i2,33]-zka[i1,i2,33]>0 and zkat[i1,i2,33]>10 and piaKa[i1,i2]<40:
            attZL.append((10**(0.1*zka[i1,i2,34:]*0.72)).sum())
            attL.append(zkat[i1,i2,33]-zka[i1,i2,33])
            pia2dL.append(piaKa[i1,i2])
            zRL.append(zka[i1,i2,:38])
            z1=zka[i1,i2,:38].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(((z1-zm)/zs)[:32])
            z1l.append((pia2dL[-1]+np.random.randn()*3-6)/8.0)
            xtfL_1.append(z1l[:32][::-1])
            xtfL_2.append(z1l[-1])
            z1L.append(z1l)
            pRateL.append(prate2d[i1,i2])
        else:
            ic40+=1
            attZ40L.append((10**(0.1*zka[i1,i2,34:]*0.72)).sum())
            att40L.append(zkat[i1,i2,33]-zka[i1,i2,33])
            pia2d40L.append(piaKa[i1,i2])
            z1=zka[i1,i2,:38].copy()
            z1[z1<0]=0
            z1+=np.random.randn()
            z1l=list(((z1-zm)/zs)[:32])
            z1l.append(5)
            z140L.append(z1l)
            xtfL_40_1.append(z1l[:32][::-1])
            xtfL_40_2.append(z1l[-1])
            pRate40L.append(prate2d[i1,i2])
    return z1L,pRateL,z140L,pRate40L,zRL,xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2

z1L,pRateL,z140L,pRate40L,zRL,xtfL_1,xtfL_2,xtfL_40_1,xtfL_40_2=readData(fname)

x_train1,x_test1,\
    ind_train, ind_test = train_test_split(np.array(xtfL_1)[:,:], np.arange(len(pRateL)),\
                                       test_size=0.5, random_state=42)

y_test=np.array(pRateL)[ind_test]
y_train=np.array(pRateL)[ind_train]

x_test2=np.array(xtfL_2)[ind_test]
x_train2=np.array(xtfL_2)[ind_train]

from combModel import lstm_model, tf

model=lstm_model(1,1)
x_train1=x_train1[:,:,np.newaxis]
x_test1=x_test1[:,:,np.newaxis]

itrain=1
if itrain==1:
    model.compile(
        optimizer=tf.keras.optimizers.Adam(),  \
        loss='mse',\
        metrics=[tf.keras.metrics.MeanSquaredError()])
    history = model.fit([x_train1,x_train2], y_train[:], \
                        batch_size=32,epochs=40,
                        validation_data=([x_test1,x_test2], y_test[:]))
else:
    model=tf.keras.models.load_model("radarProfilingKa.h5")

model.save("radarProfilingKa.h5")

yp=model.predict([x_test1,x_test2])[:,0]
for intv in [[1,2],[2,4],[4,6],[6,8],[8,10],[10,20]]:
    a=np.nonzero((y_test-intv[0])*(y_test-intv[1])<0)
    rms=(((yp[a[0]]-y_test[a[0]])**2).mean())**0.5/y_test[a[0]].mean()
    print("%6.2f %6.2f %6.2f %6.2f"%(intv[0],intv[1],\
                                     (-1+yp[a[0]].mean()/y_test[a[0]].mean())*100,rms*100))
