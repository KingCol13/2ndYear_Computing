"""

Where I actually run the code.

"""

import scipy as sp
import ballsmodule3 as bm
import time

#sim = ob.Simulation()
#sim.run(50, animate=True)

k_b = 1.38064852e-23

no_balls = 10
ball_radius = 0.1
mag_vel = 1
cont_rad = 10
spacing = 0.5

no_collisions = 2000
animate = True
histos = False

start = time.time()
print("Started simulation at: " + str(time.ctime(start)))

sim = bm.Simulation(no_balls, mag_vel=mag_vel, ball_radius=ball_radius, cont_rad=cont_rad, spacing=spacing)
print("Initial KE = " + str(sim.ke()))
print("Mean initial velocity = " + str(sim.av_velocity()))

sim.run(no_collisions , animate=animate, histos=histos)
print("Final KE = " + str(sim.ke()))

if histos:
    sim.dist_histo()
    sim.vel_histo()
    sim.rel_dist_histo()

P = sim.pressure()
V = sim.volume()
T = sim.temperature()

print("Pressure = " + str(P))
print("Temperature = " + str(T))
print("Volume = " + str(V))

end = time.time()

print("Time elapsed: " + str(end-start))

print("PV = " + str(P*V))
print("NK_bT = " + str(no_balls * k_b * T))


mfpl_predicted = V/(sp.sqrt(2)*4*no_balls*ball_radius)

print("mfpl_predicted = " + str(mfpl_predicted))

sim.mfpl_histo(50)