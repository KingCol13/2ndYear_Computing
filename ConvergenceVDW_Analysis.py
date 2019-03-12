#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:44:29 2019

@author: kingsley
"""

import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as spo

file = "nb300ncol2000VDWn_vary.npy"

k_b = 1.38064852e-23

data = sp.load(file)

xerr=0
yerr=0

P = data[0]
V = data[1]
T = data[2]
KE = data[3]
v = data[4]
x = data[5]

b = sp.pi*0.001**2

y = (x*k_b*T)/(V-x*b) - P

x = x/(V[0])

plt.grid()

#Optimize curve fit:

x_fit = sp.linspace(sp.amin(x), sp.amax(x), 1000)

#function to fit parameters for
def fit_func(x,a):
    return a*x**2

#firt guesses
p0 = [-8]

fit, cov = spo.curve_fit(fit_func, x, y, p0)

print(fit)
print("a = " + str(fit[0]*V[0]*V[0]))

plt.plot(x_fit, fit_func(x_fit, *fit), label=r"$ax^2$ Fit", color="blue")


#plotting

ax = plt.gca()

#label graph
xlabel, ylabel, title = r"$\frac{N}{V}$ $(m^{-2})$", r"$\frac{Nk_bT}{V-Nb} - P (Pa m)$", "Finding First Constant in VDW State Equation"

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)

#set limits
#plt.xlim(0,1)
#plt.ylim(0,1)

plt.plot(x, y, 'x',  color='green', label='Data Points')

plt.legend()
plt.show()