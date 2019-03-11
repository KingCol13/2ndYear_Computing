#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:44:29 2019

@author: kingsley
"""

import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as spo

#file = "nb300ncol2000cont_vary10-110.npy"
file = "nb300ncol2000no_bvary10-310.npy"
#file = "nb300ncol2000v_vary1-30.npy"

data = sp.load(file)

xerr=0
yerr=0

P = data[0]
V = data[1]
T = data[2]
KE = data[3]
v = data[4]
x = data[5]

y = (P*V)/(x*T)

plt.grid()

#Optimize curve fit:

x_fit = sp.linspace(sp.amin(x), sp.amax(x), 1000)

#function to fit parameters for
def fit_func(x,a, c):
    return a/x + c

#firt guesses
p0 = [8, 80]

fit, cov = spo.curve_fit(fit_func, x, y, p0)

print(fit)

plt.plot(x_fit, fit_func(x_fit, *fit), label=r"$\frac{A}{x} + c$ Fit", color="blue")


#plotting

ax = plt.gca()

#label graph
xlabel, ylabel, title = "$N", r"$\frac{PV}{NT}$ $(JK^{-1})$", "Convergence of Simulation to Ideal Gas Law"

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)

#set limits
#plt.xlim(0,1)
#plt.ylim(0,1)

plt.plot(x, y, 'x',  color='green', label='Data Points')

plt.legend()
plt.show()