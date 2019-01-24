import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc as erfc

import scipy.optimize


# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/POM_time_10hours.spe')
y_data = np.array(y_data)
x_data = np.arange(len(y_data))
sigma_i = []
for y_i in y_data:
    if y_i == 0:
        sigma_i.append(1)
    else:
        sigma_i.append(np.sqrt(y_i))



# ------fitshit-----


def gaussian(x, mu, sig):
    return 1 / (np.sqrt(2 * np.pi) * sig) * np.exp(-1 / 2 * ((x - mu) / sig) ** 2)


def ex_gaussian(x, mu, sig, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sig ** 2 - 2 * x)) * erfc((mu + lam * (sig ** 2) - x) / (np.sqrt(2) * sig))

x= np.linspace(0,100,1000)
y = ex_gaussian(x,50,5, 3)
plt.plot(x,y)
def f_plas_t_1(x, *params):
    a, b, c, mu, sig, lam = params
    return a * b * gaussian(x, mu, sig) + c * ex_gaussian(x,mu,sig,lam)


def f_plas_t_2(x, *params):
    a, b, c, d, mu, sig, lam1, lam2 = params
    return a * b * gaussian(x, mu, sig) + c * ex_gaussian(x,mu,sig,lam1)  + d*ex_gaussian(x,mu,sig, lam2)

# start_param_1 = (1,100000,100000,1750,10,3)
#
# fit_params, cov = scipy.optimize.curve_fit(f=f_plas_t_1, xdata=x_data[1000:], ydata=y_data[1000:], p0=start_param_1)
# a_fit, b_fit,c_fit, mu_fit, sig_fit, lam_fit = fit_params
# x_fit = np.linspace(x_data[1000],x_data[-1],1000)
# y_fit = f_plas_t_1(x_fit, *fit_params)
# print(fit_params)
#
# # # plotshit------------
# fig = plt.figure()
# ax = plt.subplot(111)
# ax.scatter(x_data[1000:], y_data[1000:], marker='.')
# ax.plot(x_fit, y_fit, 'r', linewidth=1)
# # fig.savefig('plastic_time_fit', dpi=300)
# plt.show()
