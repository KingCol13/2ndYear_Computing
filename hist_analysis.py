#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:09:52 2019

@author: kingsley
"""

import scipy as sp
import matplotlib.pyplot as plt
import ballsmodule3 as bm

vel_data = sp.load("rb0.1vel_hist.npy")

T = 6.321275270823683e+21

def hist_predict(x):
    x = sp.floor(x*39.41327649611157)/39.41327649611157
    a = 1/(k_b * T)
    b = 1/(2*k_b*T)
    dx = sp.amax(vel_data)/50
    return (a/(2*b))*(sp.exp(-b*x**2)-sp.exp(-b*(x+dx)**2))

plt.hist(vel_data, bins=50, label="Histogram of Simulation Data")

x = sp.linspace(0,sp.amax(vel_data), 1000)

plt.plot(x, 1000*300*hist_predict(x), label="Prediction Based on Temperature of System")

plt.xlabel(r"Velocity ($ms^{-1}$)")
plt.ylabel("Frequency (no units)")
plt.title("Histogram of Velocities of Balls")
plt.legend()