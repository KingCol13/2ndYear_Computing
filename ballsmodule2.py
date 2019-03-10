# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:05:57 2019

A module containing my classes and functions for balls

@author: ktc17
"""

import scipy as sp
import pylab as pl
import matplotlib.pyplot as plt
import random as rnd

class Ball:
    """
    Ball object containing position r, velocity v, mass m and radius R.
    Conains private variables to be modified by methods.
    Holds a matplotlib circle patch.
    
    """
    def __init__(self, r, v, m, R):
        """
        Initialise the ball object, using 2D arrays for position and velocity
        use floats for mass and ball radius
        Initialises a patch and member variables
        """
        if len(r) != 2 or len(v)!= 2:
            raise Exception("Two dimensions only for r and v please!")
        self.__m = m
        self.__R = R
        self.__r = sp.array(r, dtype="float64")
        self.__v = sp.array(v, dtype="float64")
        
        if self.__R > 0:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=True, color='r')
        else:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=False, color = "blue")
            
        self.__time = 0
    
    def setpatch(self, color="red", fill=1):
        """
        Set the ball patch to the ball's current location
        """
        #If it's a container, change it's patch so we can see inside.
        if self.__R > 0:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=True, color='r')
        else:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=False, color = "blue")
    def setr(self, r):
        """
        Set the ball's position to a 2D array r
        """
        if len(r) != 2:
            raise Exception("Two dimensions for r please!")
        self.__r = sp.array(r)
    def setv(self, v):
        """
        Set the ball's velocity (e.g. after collision) to new velocity v
        """
        if len(v) != 2:
            raise Exception("Two dimensions for v please!")
        self.__v = sp.array(v)
    def move(self, dt):
        """
        Move the ball through time by dt to a new position depending on it's velocity.
        Also increments the ball's time variable
        """
        self.__r+=self.__v*dt
        self.__time+=dt
        self.setpatch()
    
    def pos(self):
        """
        Return the ball's position as a 2D array.
        """
        return self.__r
    def vel(self):
        """
        Return the ball's velocity as a 2D array.
        """
        return self.__v
    def rad(self):
        """
        Return the ball's radius as a float.
        """
        return self.__R
    def mass(self):
        """
        Return the ball's mass as a float.
        """
        return self.__m
    def get_patch(self):
        """
        Return the ball's pylab patch as a 2D array.
        """
        return self.__patch
    
    def time_to_collision(self, other):

        """
        Calculate the time to collision between ball self and ball other
        this code is called 2N times (where N is number of balls) per collision
        so must be highly optimised.
        """
        
        #TODO: Optimise this a lot
        
        t_threshold = 1e-6
        v12 = self.__v - other.__v
        r12 = self.__r - other.__r
        R_collide = other.__R + self.__R
        a = sp.dot(v12, v12)
        if a==0:
            return sp.inf
        b = 2*sp.dot(r12, v12)
        c = sp.dot(r12,r12) - R_collide*R_collide
        discrim = b*b - 4*a*c
        if discrim < 0:
            return sp.inf
        else:
            t_col = [(-b-sp.sqrt(discrim))/(2*a), (-b+sp.sqrt(discrim))/(2*a)]
            if min(t_col)>t_threshold:
                return min(t_col)
            elif max(t_col)>t_threshold:
                return max(t_col)
            else:
                return sp.inf
        
    
    def collide(self, other, t_col):
        
        """
        Move two balls to touch each other, given the time of collision, t_col
        then change their respective velocities thereby simulating a collision
        """
        
        
        self.move(t_col)
        other.move(t_col)

        #save the differences in vectors to speed things up
        r12 = other.__r - self.__r
        u12 = other.__v - self.__v
        m1 = self.__m
        m2 = other.__m
        
        
        #Vector 2D elastic collision equations shamelessly ripped from
        #wikipedia elastic collision article
        v1 = self.__v + (((2*m2)*sp.dot(u12, r12))/((m1+m2)*(sp.dot(r12, r12))))*r12
        v2 = other.__v - (((2*m1)*sp.dot(u12, r12))/((m1+m2)*(sp.dot(r12, r12))))*r12
        
        #TODO: Neaten this hack to stop container moving
        if isinstance(other, Container):
            v2 = [0,0]
        if isinstance(self, Container):
            v1 = [0,0]
        self.setv(v1)
        other.setv(v2)



class Container(Ball):
    """
    Container inherits from ball and is special as it has a negative radius
    for optimisation purposes. May also be forced to not move. Given a high mass.
    Also drawn differently.
    """
    def __init__(self, m, R):
        """
        Initialise a container by calling ball's init and set up patch for container.
        """
        Ball.__init__(self, [0,0], [0,0], m ,R)
        self._patch =  pl.Circle([0,0], radius = -R, fill=False)

class Simulation():
    """
    Simulation class is a composition of container and all the balls.
    Contains methods for initialising, running and evaluating the simulation.
    """
    def __init__(self, no_balls):
        """
        Init creates the container, balls and places them.
        It also initialises the triangular matrix of collision times.
        Checks that the average momentum of the balls is 0.
        """
        cont_rad = -10
        self._container = Container(5000, cont_rad)
        self._ball_list = []
        self._ball_list.append(self._container)
        
        mag_vel = 1
        ball_radius = 0.1
        ball_mass = 1
        bpr = 1 #balls per ring
        ring_no = 0
        j=0
        r_r = 0
        #TODO: Automate r_gap and arc_gap to spread balls out
        r_gap, arc_gap = 0.5, 0.5 #arbitrary gap between the balls in radius and theta space
        for i in range(0,no_balls):

            vel = mag_vel * [rnd.random()-0.5, rnd.random()-0.5]
            if r_r >= -cont_rad-ball_radius:
                raise Exception("Ball position has exceeded maximum, tweak settings.")
            r_theta = (2*sp.pi/bpr)*j
            ball_pos = r_r*sp.array([sp.cos(r_theta), sp.sin(r_theta)])
            self._ball_list.append(Ball(ball_pos, vel, ball_mass, ball_radius))
            j+=1
            if j%bpr == 0:
                j=0
                ring_no+=1
                r_r = 2*ring_no*ball_radius + r_gap*ring_no
                bpr = sp.floor(2*sp.pi/(2*sp.arctan(ball_radius/r_r) + arc_gap/r_r))
        
        #find the total momentum of all the balls
        p_sum = 0
        for i in self._ball_list:
            p_sum += i.vel()*i.mass()
        p_average = p_sum/no_balls
        print("p_av = " + str(p_average))
        
        #remove any resultant momentum, should be small anyway (random walk)
        for i in self._ball_list:
            i.setv(i.vel()-p_average/i.mass())
        
        #add one to account for container (treated as a ball)
        self._no_balls = no_balls + 1
        
        
        #Create a time to collision array t_cols, initialised with infinities
        self._t_cols = sp.ones([self._no_balls,self._no_balls], dtype="float64")*sp.inf
        
        #Fill the t_cols array with times to collision (make it upper triangular since it would be symmetric)
        for i in range(0, self._no_balls):
            for j in range(i, self._no_balls):
                if(i!=j):
                    t_col = (self._ball_list[i].time_to_collision(self._ball_list[j]))
                    self._t_cols[i,j] = t_col
    def next_collision(self):
#        """
#        Find the next collision from the time to collision matrix and move all the balls
#        to the time of the next collision. Also update the time to collision matrix with
#        new times to next collision. Will call time_to_collision approximately 2N times.
#        """
        
        #Find time to nearest collision t_min in matrix
        t_min = sp.amin(self._t_cols)
        #Find position in matrix of t_min
        col_pos = sp.where(self._t_cols == t_min)
        #Convert matrix position to ball indices
        i_min, j_min = int(col_pos[0]), int(col_pos[1])
        
        #Collide the balls
        self._ball_list[i_min].collide(self._ball_list[j_min], t_min)
        
        #Move balls that didn't collide normally
        for i in range(0, self._no_balls):
            if i!=i_min and i!=j_min:
                self._ball_list[i].move(t_min)
        
        #Subtract t_min from t_cols for next time update (as we have moved forwards in time now compared to when t_cols was calculated)
        self._t_cols-=t_min
        
        #Update t_cols array keeping it upper triangular
        self._t_cols[i_min, j_min] = (self._ball_list[i_min].time_to_collision(self._ball_list[j_min]))
        for i in range(0, self._no_balls):
            #Make sure we don't check balls with themselves (but we checked against each other before)
            if i!=i_min and i!=j_min:
                #Conditional magic to keep t_cols triangular
                if i>i_min:
                    self._t_cols[i_min, i] = self._ball_list[i].time_to_collision(self._ball_list[i_min])
                else:
                    self._t_cols[i, i_min] = self._ball_list[i].time_to_collision(self._ball_list[i_min])
                
                if i>j_min:
                    self._t_cols[j_min, i] = self._ball_list[i].time_to_collision(self._ball_list[j_min])
                else:
                    self._t_cols[i, j_min] = self._ball_list[i].time_to_collision(self._ball_list[j_min])
                           
        
    def run(self, num_frames, animate=False, histos=False):
        """
        Run the actual simulation for number of collisions, num_frames.
        Choose to animate with animate=True.
        Iterates next_collision method in simulation.
        """
        #Create arrays for histograms
        if histos:
            self.r_histo = []
            self.r_rel_histo = []
            self.v_histo = []
            histo_start = sp.floor(num_frames/2) #frame when program will create histogram
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            for i in range(0, self._no_balls):
                ax.add_patch(self._ball_list[i].get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.01)#0.001)
            if histos and frame > histo_start:
                for i in self._ball_list:
                    self.r_histo.append(sp.sqrt(sp.dot(i.pos(), i.pos())))
                    self.v_histo.append(sp.sqrt(sp.dot(i.vel(), i.vel())))
                    for j in self._ball_list:
                        if i!=j:
                            r12 = i.pos() - j.pos()
                            self.r_rel_histo.append(sp.sqrt(sp.dot(r12, r12)))
        if animate:
            pl.show()
            
    def ke(self):
        """
        Return the total kinetic energy of the simulation by looping over all the balls
        """
        ke = 0
        for i in range(0, self._no_balls):
            ke+=0.5*self._ball_list[i].mass()*sp.dot(self._ball_list[i].vel(), self._ball_list[i].vel())
        return ke
    def vel_histo(self):
        """
        Plot a velocity histogram to show the speed of all the balls at final frame.
        (should approximate a Maxwell-Boltzmann distribution)
        """
        plt.figure()
        plt.title("Velocity Histogram")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Frequency (No Units)")
        plt.hist(self.v_histo, bins=50)
    def dist_histo(self):
        """
        Plot a distance histogram to show the radial position of all the balls at final frame.
        """
        plt.figure()
        plt.title("Radial Distance Histogram")
        plt.xlabel("Distance (m)")
        plt.ylabel("Frequency (No Units)")
        plt.hist(self.r_histo, bins=50)
    def rel_dist_histo(self):
        """
        Plot a histogram of relative distance between each pair of balls
        """
        plt.figure()
        plt.title("Relative Distance Histogram")
        plt.xlabel("Distance (m)")
        plt.ylabel("Frequency (No Units)")
        plt.hist(self.r_rel_histo, bins=50)