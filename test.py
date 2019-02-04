import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc as erfc
from scipy.integrate import quad


def gaussian(x, mu, sig):
    return 1 / (np.sqrt(2 * np.pi) * sig) * np.exp(-1 / 2 * ((x - mu) / sig) ** 2)

def ex_gaussian(x, mu, sig, lam):
    return lam / 2 * np.exp(lam / 2 * (2 * mu + lam * sig ** 2 - 2 * x)) * erfc(
        (mu + lam * (sig ** 2) - x) / (np.sqrt(2) * sig))

params = [5.29225661e+00, 7.06560109e+05, 3.49827036e+05, 3.35635147e+05,
 1.71655407e+03, 1.07879049e+01, 4.78193371e-02, 1.12908754e-02]

b = params[1]
c = params[2]
mu = params[4]
sig = params[5]
lam1 = params[6]
lam2 = params[7]


int_gauss = b*quad(gaussian, 0, 4000, args=(mu,sig))[0]
int_ex_gauss = c*quad(ex_gaussian,0,4000, args=(mu, sig, lam1))[0]
print(int_ex_gauss)
x = np.linspace(0,4000,1000)
y = gaussian(x, mu,sig)
plt.plot(x,y )
plt.show()