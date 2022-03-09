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
zKuL=[]
R=287
h=(250/2.+np.arange(76)*250)/1e3
h1=(0.+np.arange(77)*250)/1e3

for f in fs:
    fh=Dataset(f)
    qr=fh["QRAIN"][:]
    qv=fh["QVAPOR"][:]
    qh=fh["QHAIL"][:]
    qs=fh["QSNOW"][:]
    qg=fh["QGRAUPLE"][:]
    ph=fh["PH"][:]+fh["PHB"][:]
    press=fh["P"][:]+fh["PB"][:]
    T=fh["T"]
