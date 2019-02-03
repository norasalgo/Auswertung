import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc as erfc
import scipy.optimize

# ---Load Data ---------------

path_data = './data_notext'

y_data = np.loadtxt(path_data + '/Al_time_20min.spe')
y_data = np.array(y_data)
x_data = np.arange(len(y_data))[1000:]
y_data = y_data[1000:]
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
        (mu + lam * sig ** 2 - x) / (np.sqrt(2) * sig))


def f_al_1(x, *params):
    """
fit function with one gauss
    """
    a, b, mu, sig = params
    return a + b * gaussian(x, mu, sig)

def f_al_2(x, *params):
    """
fit function for one gauss and one convoluted gauss + exp decay
    :param x:
    :param params: 6 tuple : a,b,c,mu,sig,lam
    """
    a, b, c, mu, sig, lam = params
    return a + b * gaussian(x,mu,sig) + c*ex_gaussian(x,mu,sig,lam)

#fit gauss
start_params_1 = (1, 1, 1900, 10)
fit_params_1, cov_1 = scipy.optimize.curve_fit(f=f_al_1, xdata=x_data, ydata=y_data, p0=start_params_1)
a_fit_1, b_fit_1, mu_fit_1, sig_fit_1 = fit_params_1
err_1 = np.sqrt(np.diag(cov_1))
print('fit params1',fit_params_1)
print('err1', err_1)

x_fit_1 = np.linspace(start=x_data[0], stop=x_data[-1], num=10 * len(x_data))
y_fit_1 = f_al_1(x_fit_1, *fit_params_1)

#fit gauss + convol
start_params_2 = (1, -70000, 10000, 1700, 10, 0.1)
fit_params_2, cov_2 = scipy.optimize.curve_fit(f=f_al_2, xdata=x_data, ydata=y_data, p0=start_params_2)
a_fit_2, b_fit_2, c_fit_2, mu_fit_2, sig_fit_2, lam_fit_2 = fit_params_2
err_2 = np.sqrt(np.diag(cov_2))
tau_fit = 1/(lam_fit_2*alpha)

print('fit params2',fit_params_2)
print('err2', err_2)
print('tau fit', tau_fit)

x_fit_2 = np.linspace(start=x_data[0], stop=x_data[-1], num=10 * len(x_data))
y_fit_2 = f_al_2(x_fit_2, *fit_params_2)




# plotshit------------
fig = plt.figure()
ax_1 = plt.subplot(211)
ax_1.set_xlabel('Channel number')
ax_1.set_ylabel('Counter')
ax_1.scatter(x_data[625:900], y_data[625:900], marker='.', s=0.7)
ax_1.plot(x_fit_1[6250:9000], y_fit_1[6250:9000], 'r', linewidth=1)

ax_2 = plt.subplot(212)
ax_2.set_xlabel('Channel number')
ax_2.set_ylabel('Counter')
ax_2.scatter(x_data[625:900], y_data[625:900], marker='.', s= 0.7)
ax_2.plot(x_fit_2[6250:9000], y_fit_2[6250:9000], 'r', linewidth=1)
# fig.savefig('al_time_fit.png', dpi=300)
# plt.show()