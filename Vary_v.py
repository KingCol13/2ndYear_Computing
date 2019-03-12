#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:44:29 2019

@author: kingsley
"""

import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as spo

#file = "nb300ncol2000T_vary.npy"
file = "nb300ncol2000cont_vary10-110.npy"

data = sp.load(file)

P = data[0]
V = data[1]
T = data[2]
KE = data[3]
v = data[4]

y = P
x = V

plt.grid()

#Optimize curve fit:

x_fit = sp.linspace(sp.amin(x), sp.amax(x), 1000)

#function to fit parameters for
def fit_func(x,a):
    return a/x

#firt guesses
p0 = [1e23]

fit, cov = spo.curve_fit(fit_func, x, y, p0)

print(fit)
print(sp.sqrt(cov))

plt.plot(x_fit, fit_func(x_fit, *fit), label=r"$\frac{a}{x}$ Fit", color="blue")


#plotting

ax = plt.gca()

#label graph
xlabel, ylabel, title = r"Volume ($m^{2}$)", r"Pressure (Pa m)", "Reciprocal Dependence of P on V"

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)

#set limits
#plt.xlim(0,1)
#plt.ylim(0,1)

plt.plot(x, y, 'x',  color='green', label='Data Points')

plt.legend()
plt.show()