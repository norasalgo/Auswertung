import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc as erfc

import scipy.optimize


# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/timecalibration.spe')
y_data = np.array(y_data)
x_data = np.arange(len(y_data))


#peaks--------------------

p_1 = [900,1000]


#plot----------

fig = plt.figure()
plt.plot(x_data,y_data/10e4)
plt.xlabel('Channel number')
plt.ylabel(r'Event counts $[10^4]$')

plt.plot(x_data[p_1[0]:p_1[1]],y_data[p_1[0]:p_1[1]],'r')
plt.show()
#fig.savefig('time_calibration_data.png', dpi = 600)
