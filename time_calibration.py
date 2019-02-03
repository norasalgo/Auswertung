import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc as erfc

import scipy.optimize


# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/timecalibration.spe')
y_data = np.array(y_data)
y_data = y_data/10e4
x_data = np.arange(len(y_data))

#peaks--------------------
cuts = [700,1000,1300,1700,2000,2300,2650,3000,3300,3650,4000]

weighted_means = []
for i in range(len(cuts)-1):
    wm = sum(y_data[cuts[i]:cuts[i+1]]*x_data[cuts[i]:cuts[i+1]])/ sum(y_data[cuts[i]:cuts[i+1]])
    weighted_means.append(wm)

var_wm = []
for i in range(len(cuts)-1 ):
    var_i = 0
    for j in range(cuts[i],cuts[i+1]):
        var_i += (x_data[j] - weighted_means[i])**2*y_data[j]
    var_wm.append(var_i/sum(y_data[cuts[i]:cuts[i+1]]))

print(weighted_means)
print(var_wm)
#Plot every peak to make sure we have the right cuts
plt.plot(x_data,y_data)
plt.xlabel('Channel number')
plt.ylabel(r'Event counts $[10^4]$')
for i in range(len(cuts)-1):
    plt.plot(x_data[cuts[i]:cuts[i+1]], y_data[cuts[i]:cuts[i+1]])
    plt.axvline(x=weighted_means[i], color='b', linestyle='--', linewidth = 0.5)
plt.show()


#calculate alpha factor (bins/ns)

peak_distances = np.array([weighted_means[i+1]-weighted_means[i] for i in range(len(weighted_means)-1)])
mean_distance = sum(peak_distances)/len(peak_distances)
var_dist = sum([(p-mean_distance)**2 for p in peak_distances]) /len(peak_distances)
alpha = mean_distance/8
err_alpha = np.sqrt(var_dist)/8
print('alpha',alpha)
print('alpha_err', err_alpha)

#plot----------
#
# fig = plt.figure()
# plt.plot(x_data,y_data)
# plt.xlabel('Channel number')
# plt.ylabel(r'Event counts $[10^4]$')
#
# #plt.plot(x_data[p_1[0]:p_1[1]],y_data[p_1[0]:p_1[1]],'r')
#
#
#
# plt.show()
# #fig.savefig('time_calibration_data.png', dpi = 600)
