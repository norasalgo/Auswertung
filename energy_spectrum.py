import matplotlib.pyplot as plt
import numpy as np
from math import erfc as erfc
import scipy.optimize

# ---Load Data ---------------

path_data = './data_notext'


y_data_1 = np.loadtxt(path_data + '/PMT1_energy_nogate.spe')[40:]
y_data_2 = np.loadtxt(path_data + '/PMT2_energy_nogate.spe')[40:]
x_data = np.arange(len(y_data_1))


plt.plot(x_data,y_data_2)
plt.yscale('log')
plt.show()