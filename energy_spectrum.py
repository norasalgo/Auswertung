import matplotlib.pyplot as plt
import numpy as np
from math import erfc as erfc
import scipy.optimize

# ---Load Data ---------------

path_data = './data_notext'


y_pmt1_nogate = np.loadtxt(path_data + '/PMT1_energy_nogate.spe')[40:]
y_pmt2_nogate = np.loadtxt(path_data + '/PMT2_energy_nogate.spe')[40:]
y_pmt1_withgate = np.loadtxt(path_data + '/PMT1_energy1275_withgate.spe')
y_pmt2_withgate = np.loadtxt(path_data + '/PMT2_energy511_withgate.spe')

fig, axes = plt.subplots(nrows=2, ncols=2)
#PMT1 No gate
axes[0 ,0 ].scatter(np.arange(len(y_pmt1_nogate)),y_pmt1_nogate, s=0.6)
axes[0 ,0 ].set_xlabel('Channel number')
axes[0 ,0 ].set_ylabel('Event counts')

#PMT2 no gat
axes[0 , 1].scatter(np.arange(len(y_pmt2_nogate)),y_pmt2_nogate,s=0.8)
axes[0 , 1].set_xlabel('Channel number')
axes[0 , 1].set_ylabel('Event counts')

#PMT1 withgate, 1254keV
axes[1,0].scatter(np.arange(len(y_pmt1_withgate)),y_pmt1_withgate, s=0.6)
axes[1,0].set_xlabel('Channel number')
axes[1,0].set_ylabel('Event counts')

#PMT2 withgate, 511keV
axes[1,1].scatter(np.arange(len(y_pmt2_withgate)),y_pmt2_withgate, s= 0.8)
axes[1,1].set_xlabel('Channel number')
axes[1,1].set_ylabel('Event counts')


fig.tight_layout()
plt.show()
fig.savefig('energy_spectrum.png', dpi = 600)

