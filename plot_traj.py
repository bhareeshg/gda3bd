import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

datsec=pd.read_csv(sys.argv[1],skiprows=2);
esec=datsec.e;
isec=datsec.i*180/np.pi;
osec=datsec.omega*180/np.pi;
Osec=datsec.Omega*180/np.pi;

cmap = plt.cm.get_cmap('RdYlBu');

fig, axs = plt.subplots(2,2);

axs[0][0].plot(datsec.time,esec,markersize=0.01,linestyle='--',alpha=0.6);
axs[0][1].plot(datsec.time,isec,markersize=0.01,linestyle='--',alpha=0.6);
axs[1][0].plot(datsec.time,osec,markersize=0.01,linestyle='--',alpha=0.6);
axs[1][1].plot(datsec.time,Osec,markersize=0.01,linestyle='--',alpha=0.6,label='Secular');
axs[1][1].legend();

axs[1][0].set_xlabel(r't(years)')
axs[1][1].set_xlabel(r't(years)')

axs[0][0].set_ylabel(r'$e$')
axs[0][1].set_ylabel(r'$i$')
axs[1][0].set_ylabel(r'$\omega$')
axs[1][1].set_ylabel(r'$\Omega$')


fig.subplots_adjust(wspace=0.39)

sp=[axs[0][0],axs[0][1]]
plt.setp([a.get_xticklabels() for a in sp], visible=False)


filename="trajectory.pdf";
plt.savefig(filename,bbox_inches='tight')

