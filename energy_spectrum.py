import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.special.erfc as erfc

# ---Load Data ---------------

path_data = './data_notext'

data = np.loadtxt(path_data + '/Al_time_20min.spe')


# ------Fitshit-----


def gaussian(x,mu, sigma):
    return 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-1/2*((x-mu)/sigma)**2)


def ex_gaussian(x, mu, sigma, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sigma ** 2 - 2 * x)) * erfc(
        (mu + lam * sigma ** 2 - x) / (np.sqrt(2) * sigma))


x = range(len(data))
plt.scatter(x, data)
plt.show()
print(data[20:40])
