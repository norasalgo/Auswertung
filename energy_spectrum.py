import matplotlib.pyplot as plt
import numpy as np
from math import erfc as erfc
import scipy.optimize
# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/Al_time_20min.spe')
y_data = np.array(y_data)
x_data = np.arange(len(y_data))[1000:]
y_data = y_data[1000:]
sigma_i = np.sqrt(y_data)

# ------fitshit-----


def gaussian(x, mu, sig):
    return 1 / (np.sqrt(2*np.pi) * sig) * np.exp(-1 / 2 * ((x - mu) / sig) ** 2)


def ex_gaussian(x, mu, sig, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sig ** 2 - 2 * x)) * erfc(
        (mu + lam * sig ** 2 - x) / (np.sqrt(2) * sig))

#
# def f_al (x, a, b, mu, sig):
#     return a * b*gaussian(x, mu, sig)

def f_al (x, *params):
    a, b, mu, sig = params
    return a * b*gaussian(x, mu, sig)



fit_params, cov = scipy.optimize.curve_fit(f=f_al,xdata=x_data,ydata=y_data, p0=(1,1,1900,10))
# a_fit, b_fit, mu_fit, sig_fit = fit_params
y_fit = f_al(x_data,*fit_params)
print(fit_params)

#plotshit------------
fig = plt.figure()
ax = plt.subplot(111)
ax.scatter(x_data,y_data, marker = '.')
ax.scatter(x_data, y_fit,marker = ".")
fig.savefig('al_time_fit.png', dpi = 300)
plt.show()
