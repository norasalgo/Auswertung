import matplotlib.pyplot as plt
import numpy as np
from math import erfc as erfc
import scipy.optimize
# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/Al_time_20min.spe')
y_data = np.array(y_data)
x_data = np.arange(len(y_data))
sigma_i = np.sqrt(y_data)

# ------fitshit-----


def gaussian(x, mu, sig):
    return 1 / (np.sqrt(2*np.pi) * sig) * np.exp(-1 / 2 * ((x - mu) / sig) ** 2)


def ex_gaussian(x, mu, sig, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sig ** 2 - 2 * x)) * erfc(
        (mu + lam * sig ** 2 - x) / (np.sqrt(2) * sig))


def f_al (x, a, b, mu, sig):
    return a * b*gaussian(x, mu, sig)

def chi2():
    pass

y_fit, cov = scipy.optimize.curve_fit(f=f_al,xdata=x_data,ydata=y_data)


#plotshit------------
