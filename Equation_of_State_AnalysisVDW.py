#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:44:29 2019

@author: kingsley
"""

import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as spo

file = "br0.1nb300ncol2000T_vary.npy"

data = sp.load(file)

P = data[0]
V = data[1]
T = data[2]
KE = data[3]
v = data[4]

b = sp.pi*0.1**2

y = P
x = (300*T)/V
#x = (300*T)/(V-300*b)

plt.grid()

#Optimize curve fit:

x_fit = sp.linspace(sp.amin(x), sp.amax(x), 1000)

#function to fit parameters for
def fit_func(x,a, c):
    return a*x - c

#firt guesses
p0 = [8, 80]

fit, cov = spo.curve_fit(fit_func, x, y, p0)

print(fit)
print(sp.sqrt(cov))

print(fit[1]*(V[0]/300)**2)

plt.plot(x_fit, fit_func(x_fit, *fit), label=r"$ax-c$ Fit", color="blue")


#plotting

ax = plt.gca()

#label graph
xlabel, ylabel, title = r"$\frac{NT}{V-Nb} (Km^{-2}$)", r"$P$ (Pa m)", "Van der Waals Equation of State"

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title)

#set limits
#plt.xlim(0,1)
#plt.ylim(0,1)

plt.plot(x, y, 'x',  color='green', label='Data Points')

plt.legend()
plt.show()