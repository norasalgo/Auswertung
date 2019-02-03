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

alpha = 41.15952880418127  # time conversion


# ------fitshit-----


def gaussian(x, mu, sig):
    return 1 / (np.sqrt(2 * np.pi) * sig) * np.exp(-1 / 2 * ((x - mu) / sig) ** 2)


def ex_gaussian(x, mu, sig, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sig ** 2 - 2 * x)) * erfc(
        (mu + lam * (sig ** 2) - x) / (np.sqrt(2) * sig))


# x= np.linspace(0,100,1000)
# y = ex_gaussian(x,50,5, 0.1)
# plt.plot(x,y)
# plt.show()
def f_plas_t_1(x, *params):
    a, b, c, mu, sig, lam = params
    return a * b * gaussian(x, mu, sig) + c * ex_gaussian(x, mu, sig, lam)


def f_plas_t_2(x, *params):
    a, b, c, d, mu, sig, lam1, lam2 = params
    return a * b * gaussian(x, mu, sig) + c * ex_gaussian(x, mu, sig, lam1) + d * ex_gaussian(x, mu, sig, lam2)


# def f_plas_t_3(x, *params):
#     a, b, c, d ,e ,mu, sig, lam1, lam2 = params
#     lam3 = (lam1 + lam2)
#     return a * b * gaussian(x, mu, sig) + c * ex_gaussian(x, mu, sig, lam1) + d * ex_gaussian(x, mu, sig, lam2) + e * ex_gaussian(x, mu, sig, lam3)


start_param_1 = (1, 100000, 100000, 1750, 10, 0.1)
start_param_2 = (40, 100000, 1000000, 1000000, 1750, 10, 0.1, 0.5)
# start_param_3 = (40, 100000, 1000000, 1000000, -10000, 1750, 10, 0.1, 0.5)

bounds2 = ([-np.inf,-np.inf,-np.inf, 0 ,-np.inf,-np.inf,-np.inf,-np.inf],
        [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
fit_params_1, cov_1 = scipy.optimize.curve_fit(f=f_plas_t_1, xdata=x_data[1000:], sigma=sigma_i[1000:],
                                               ydata=y_data[1000:], p0=start_param_1)
a_fit, b_fit, c_fit, mu_fit, sig_fit, lam_fit = fit_params_1
err_1 = np.sqrt(np.diag(cov_1))
print('params1=', fit_params_1)
print('err1=', err_1)

fit_params_2, cov_2 = scipy.optimize.curve_fit(f=f_plas_t_2, xdata=x_data[1000:], ydata=y_data[1000:], p0=start_param_2,
                                               sigma=sigma_i[1000:], bounds=bounds2)
a_fit_2, b_fit_2, c_fit_2, d_fit_2, mu_fit_2, sig_fit_2, lam_1_fit_2, lam_2_fit_2 = fit_params_2
err_2 = np.sqrt(np.diag(cov_2))
print('params_2=', fit_params_2)
print('err2=', err_2)
tau_1 = 1 / (lam_1_fit_2 * alpha)
tau_2 = 1 / (lam_2_fit_2 * alpha)
print('ta1&2 fit2=', tau_1, tau_2)

# fit_params_3, cov_3 = scipy.optimize.curve_fit(f=f_plas_t_3, xdata= x_data[1000:], ydata=y_data[1000:], p0= start_param_3, sigma=sigma_i[1000:])
# a_fit_3, b_fit_3, c_fit_3, d_fit_3, e_fit_3, mu_fit_3, sig_fit_3, lam_1_fit_3, lam_2_fit_3 = fit_params_3
# print(1/(lam_1_fit_3*alpha),1/(lam_2_fit_3*alpha))

# # plotshit---------------------------------------------------------------------
x_start = 1600
x_stop = 2000
x_fit_1 = np.linspace(x_start, x_stop, 1000)
y_fit_1 = f_plas_t_1(x_fit_1, *fit_params_1)

x_fit_2 = np.linspace(x_start, x_stop, 1000)
y_fit_2 = f_plas_t_2(x_fit_2, *fit_params_2)

# x_fit_3 = np.linspace(x_start, x_stop, 1000)
# y_fit_3 = f_plas_t_3(x_fit_3, *fit_params_3)


fig = plt.figure()
# ax_1 = plt.subplot(311)
# ax_1.set_xlabel('chanel number')
# ax_1.set_ylabel('counter')
# ax_1.scatter(x_data[x_start:x_stop], y_data[x_start:x_stop], marker='.', s=0.8, label='data')
# ax_1.plot(x_fit_1, y_fit_1, 'r', linewidth=0.5, label='fit')
# ax_1.plot(x_fit_1, b_fit * gaussian(x_fit_1, mu_fit, sig_fit), 'g--', linewidth=0.5, label='gaussian')
# ax_1.plot(x_fit_1, len(x_fit_1) * [a_fit], 'b--', linewidth=0.5, label='Background')
# ax_1.plot(x_fit_1, c_fit * ex_gaussian(x_fit_1, mu_fit, sig_fit, lam_fit), 'm--', linewidth=0.5, label='Convolution')
# ax_1.legend(loc=1)

ax_2 = plt.subplot(312)
ax_2.set_xlabel('chanel number')
ax_2.set_ylabel('counter')
ax_2.scatter(x_data[x_start:x_stop], y_data[x_start:x_stop], marker='.', s=0.8, label='Data')
ax_2.plot(x_fit_2, y_fit_2, 'r', linewidth=0.5, label='Fit (Total)')
ax_2.plot(x_fit_2, b_fit_2 * gaussian(x_fit_2, mu_fit_2, sig_fit_2), 'g--', linewidth=0.5, label='Fit: Gaussian')
ax_2.plot(x_fit_2, len(x_fit_2) * [a_fit_2], 'b--', linewidth=0.5, label='Fit: Background')
ax_2.plot(x_fit_2, c_fit_2 * ex_gaussian(x_fit_2, mu_fit_2, sig_fit_2, lam_1_fit_2), 'm--', linewidth=0.5,
          label='Fit: Convolution1')
ax_2.plot(x_fit_2, d_fit_2 * ex_gaussian(x_fit_2, mu_fit_2, sig_fit_2, lam_2_fit_2), 'm--', linewidth=0.5,
          label='Fit: Convolution2')
ax_2.legend(loc=1)

# #Fit mit drei convolutions
# ax_3 = plt.subplot(313)
# ax_3.set_xlabel('chanel number')
# ax_3.set_ylabel('counter')
# ax_3.scatter(x_data[x_start:x_stop], y_data[x_start:x_stop], marker='.', s=0.8, label='Data')
# ax_3.plot(x_fit_3, y_fit_3, 'r', linewidth=0.5, label='Fit (Total)')
# ax_3.plot(x_fit_3, b_fit_3 * gaussian(x_fit_2, mu_fit_2, sig_fit_2), 'g--', linewidth=0.5, label='Fit: Gaussian')
# ax_3.plot(x_fit_3, len(x_fit_3) * [a_fit_3], 'b--', linewidth=0.5, label='Fit: Background')
# ax_3.plot(x_fit_3, c_fit_3 * ex_gaussian(x_fit_3, mu_fit_3, sig_fit_3, lam_1_fit_3), 'm--', linewidth=0.5,
#           label='Fit: Convolution1')
# ax_3.plot(x_fit_3, d_fit_3 * ex_gaussian(x_fit_3, mu_fit_3, sig_fit_3, lam_2_fit_3), 'm--', linewidth=0.5,
#           label='Fit: Convolution2')
# ax_3.plot(x_fit_3, e_fit_3 * ex_gaussian(x_fit_3, mu_fit_3, sig_fit_3, (lam_1_fit_3 + lam_2_fit_3)), 'm--', linewidth=0.5,
#           label='Fit: Convolution3')
# ax_3.legend(loc=1)

plt.show()
# fig.savefig('plastic_time_fit_both.png', dpi=600)
