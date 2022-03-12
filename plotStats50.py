import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11})

import matplotlib.pyplot as plt
l1=open("scores50.txt","r").readlines()
stats=[]
for l in l1[1:7]:
    stats.append([float(v) for v in l.split()])
for l in l1[8:]:
    stats.append([float(v) for v in l.split()])

stats=np.array(stats)
plt.figure()
plt.bar(np.arange(1,6)-0.225,stats[1:6,2],0.45,label='Bias')
plt.bar(np.arange(1,6)+0.225,stats[1:6,3],0.45,label='NRMS')

#plt.xticks(range(7),('1.0-2.0','2.0-4.0','4.0-6.0','6.0-8.0','8.0-10.0','10.0-20.0','20.0-30.0'))
plt.xticks(range(1,6),('1.0-2.0','2.0-4.0','4.0-6.0','6.0-8.0','8.0-10.0'))
plt.ylabel('%')
plt.legend()
plt.gca().yaxis.grid()
plt.xlabel('mm/h')
plt.title('Ku-band inclined\n50 km resolution')
plt.savefig('kuBand_inclined_50km.png')

plt.figure()

plt.bar(np.arange(1,6)-0.225,stats[7:14,2],0.45,label='Bias')
plt.bar(np.arange(1,6)+0.225,stats[7:14,3],0.45,label='NRMS')

#plt.xticks(range(7),('1.0-2.0','2.0-4.0','4.0-6.0','6.0-8.0','8.0-10.0','10.0-20.0','20.0-30.0'))
plt.xticks(range(1,6),('1.0-2.0','2.0-4.0','4.0-6.0','6.0-8.0','8.0-10.0'))
plt.ylabel('%')
plt.legend()
plt.gca().yaxis.grid()
plt.xlabel('mm/h')
plt.title('Ka-band polar\n50 km resolution')
plt.savefig('kaBand_polar_50km.png')


l1=open("scores.txt","r").readlines()
stats=[]
for l in l1[1:7]:
    stats.append([float(v) for v in l.split()])


stats=np.array(stats)
plt.figure()
plt.bar(np.arange(0,6)-0.225,stats[0:6,2],0.45,label='Bias')
plt.bar(np.arange(0,6)+0.225,stats[0:6,3],0.45,label='NRMS')
plt.xticks(range(6),('1.0-2.0','2.0-4.0','4.0-6.0','6.0-8.0','8.0-10.0','10-20'))
plt.ylabel('%')
plt.legend()
plt.gca().yaxis.grid()
plt.xlabel('mm/h')
plt.title('Ka-band polar\nProfile resolution')
plt.savefig('kaBand_inclined_profile.png')
