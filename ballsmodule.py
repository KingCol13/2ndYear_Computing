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
        if isinstance(self, Container) > 0:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=False, color = "blue")
        else:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=True, color='r')
    def setr(self, r):
        if len(r) != 2:
            raise Exception("Two dimensions for r please!")
        self.__r = sp.array(r)
    def setv(self, v):
        if len(v) != 2:
            raise Exception("Two dimensions for v please!")
        self.__v = sp.array(v)
    def move(self, dt):
        self.__r+=self.__v*dt
        self.__time+=dt
        self.setpatch()
    
    def pos(self):
        return self.__r
    def vel(self):
        return self.__v
    def rad(self):
        return self.__R
    def mass(self):
        return self.__m
    def get_patch(self):
        return self.__patch
    
    def time_to_collision(self, other):
        
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
    def __init__(self, m, R):
        Ball.__init__(self, [0,0], [0,0], m ,R)
        self._patch =  pl.Circle([0,0], radius = 10, fill=False)

class Simulation():
    def __init__(self, no_balls):
        cont_rad = -10
        self._container = Container(5000, cont_rad)
        self._ball_list = []
        self._ball_list.append(self._container)
        #TODO: Fix initial placement
        mag_vel = 1
        ball_radius = 0.1
        ball_mass = 1
        bpr = 6 #balls per ring
        r_gap, theta_gap = 0.1, 0.2 #arbitrary gap between the balls in radius and theta space
        #TODO: Change placement implementation to increase bpr when r_r increases
        for i in range(0,no_balls):
            vel = mag_vel * [rnd.random()-0.5, rnd.random()-0.5]
            r_r = 2*sp.floor(i/bpr+1)*ball_radius + r_gap*sp.floor(i/6+1)
            if r_r >= -cont_rad-ball_radius:
                raise Exception("Ball position has exceeded maximum, tweak settings.")
            #bpr = sp.floor(2*sp.pi/(2*sp.arctan(ball_radius/r_r)+theta_gap))-1
            r_theta = (2*sp.pi/bpr)*(i%bpr)
            ball_pos = r_r*sp.array([sp.cos(r_theta), sp.sin(r_theta)])
            self._ball_list.append(Ball(ball_pos, vel, ball_mass, ball_radius))
        #add one to account for container (treated as a ball)
        self._no_balls = no_balls + 1
        
        
    def next_collision(self):
        #TODO: Implement O(n) optimisation by only checking last collided balls
        t_min = float("inf")
        t_cols = sp.ones([self._no_balls,self._no_balls], dtype="float64")*sp.inf
        for i in range(0, self._no_balls):
            for j in range(i, self._no_balls):
                if(i!=j):
                    t_col = (self._ball_list[i].time_to_collision(self._ball_list[j]))
                    t_cols[i,j] = t_col
                    if(t_col<t_min):
                        i_min, j_min = i, j
                        t_min = t_col
#        print(t_cols)
#        print("t_min = " + str(t_min))
#        print("min(t_cols = " + str(sp.amin(t_cols)))
        self._ball_list[i_min].collide(self._ball_list[j_min], t_min)
        for i in range(0, self._no_balls):
            if i!=i_min and i!=j_min:
                self._ball_list[i].move(t_min)
    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            for i in range(0, self._no_balls):
                ax.add_patch(self._ball_list[i].get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.002)#0.001)
        if animate:
            pl.show()
            
    def ke(self):
        ke = 0
        for i in range(0, self._no_balls):
            ke+=0.5*self._ball_list[i].mass()*sp.dot(self._ball_list[i].vel(), self._ball_list[i].vel())
        return ke
    def vel_histo(self):
       plt.figure()
       x = []
       for i in self._ball_list:
           x.append(sp.sqrt(sp.dot(i.vel(), i.vel())))
       plt.hist(x)
    def dist_histo(self):
       plt.figure()
       x = []
       for i in self._ball_list:
           x.append(sp.sqrt(sp.dot(i.pos(), i.pos())))
       plt.hist(x)