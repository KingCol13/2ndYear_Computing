# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:05:57 2019

A module containing my classes and functions for balls

@author: ktc17
"""

import scipy as sp
import pylab as pl

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
        
        #TODO: Fix this hack making container radius negative
        if isinstance(self, Container) > 0:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=False, color = "blue")
        else:
            self.__patch = pl.Circle(self.__r, radius = self.__R, fill=True, color='r')
        
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
        v12 = self.vel() - other.vel()
        r12 = self.pos() - other.pos()
        if isinstance(other, Container):
            R_collide =  other.rad() - self.rad()
        else:
            R_collide = other.rad() + self.rad()
        coeff = [sp.dot(v12, v12), 2*sp.dot(r12, v12), sp.dot(r12,r12) - R_collide*R_collide]
        t_col = sp.roots(coeff)
        print(t_col)
        if t_col[0] > 0 and sp.isreal(t_col[0]):
            print("Returning t_col[0]")
            return t_col[0]
        elif t_col[1] > 0 and sp.isreal(t_col[1]):
            print("Returning t_col[1]")
            return t_col[1]
        else:
            return -1
        
    
    def collide(self, other, t_col, dt):
        self.move(t_col)
        other.move(t_col)
        normal = (self.pos() - other.pos())
        normal= normal/sp.sqrt(sp.dot(normal, normal))
        m1 = self.mass()
        m2 = other.mass()
        u1perp = sp.dot(self.vel(), normal)
        u2perp = sp.dot(other.vel(), normal)
        v1perp = ((m1-m2)*u1perp + (2*m2)*u2perp)/(m1+m2)
        v2perp = ((m2-m1)*u2perp + (2*m1)*u1perp)/(m1+m2)
        v1 = self.vel() + (v1perp - u1perp)*normal
        v2 = self.vel() + (v2perp - u2perp)*normal
        #TODO: Neaten this hack to stop container moving
        if isinstance(other, Container):
            v2 = [0,0]
        self.setv(v1)
        other.setv(v2)
        self.move(dt-t_col)
        other.move(dt-t_col)



class Container(Ball):
    def __init__(self, m, R):
        Ball.__init__(self, [0,0], [0,0], m ,R)
        self._patch =  pl.Circle([0,0], radius = 10, fill=False)

class Simulation():
    def __init__(self):
        cont_rad = 10
        self._container = Container(5000, cont_rad)
        self._ball = Ball([0.548547,1.8658765], [1.87537843,0.4359834], 1, 1)
    def next_collision(self):
        t_col = self._ball.time_to_collision(self._container)
        if t_col != -1:
            self._ball.collide(self._container, t_col, t_col)
    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_artist(self._container.get_patch())
            ax.add_patch(self._ball.get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.2)#0.001)
        if animate:
            pl.show()
        