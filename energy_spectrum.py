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

fig = plt.figure()

ax_1 = plt.subplot(221)
ax_1.scatter(np.arange(len(y_pmt1_nogate)),y_pmt1_nogate, s=0.8)
ax_2 = plt.subplot(222)
ax_2.scatter(np.arange(len(y_pmt2_nogate)),y_pmt2_nogate,s=0.8)
ax_3 = plt.subplot(223)
ax_3.scatter(np.arange(len(y_pmt1_withgate)),y_pmt1_withgate, s=0.8)
ax_4 = plt.subplot(224)
ax_4.scatter(np.arange(len(y_pmt2_withgate)),y_pmt2_withgate, s= 0.8)
# plt.yscale('log')
plt.show()