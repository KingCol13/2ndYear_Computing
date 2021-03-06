# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:11:56 2019

@author: ktc17
"""

import scipy as sp
import ballsmodule3 as bm
import time
import matplotlib.pyplot as plt

k_b = 1.38064852e-23

no_balls = 300
no_collisions = 2000
animate = False
histos = False

start = time.time()
print("Started simulation at: " + str(time.ctime(start)))


P = []
V = []
T = []
KE = []
v = []

no_b = sp.linspace(10, 310, 51)
b_r = sp.linspace(0.0001, 1.0001, 51)
c_r = sp.linspace(10, 110, 51)
t_vary = sp.sqrt(sp.linspace(1, 100, 100))

sims = 31

j=0

# used ball_radius = 0.001 for IGL data when not varying it

for i in c_r:
    sim = bm.Simulation(no_balls, ball_radius=0.1, mag_vel=i, cont_rad=10, spacing=0.7)
    j+=1
    print("Running simulation " + str(j) + " of " + str(len(t_vary)) + ".")
    sim.run(no_collisions , animate=animate, histos=histos)
    
    P.append(sim.pressure())
    V.append(sim.volume())
    T.append(sim.temperature())
    KE.append(sim.ke())
    v.append(sim.av_velocity())



P = sp.array(P)
V = sp.array(V)
T = sp.array(T)
KE = sp.array(KE)
v = sp.array(v)

end = time.time()
print("Time elapsed: " + str(end-start))

data = [P, V, T, KE, v, no_b]

sp.save("br0.1nb300ncol2000V_vary.np", data)