#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 19:48:49 2020

@author: duypham

Inspired by https://gist.github.com/tstreamDOTh/4af1d6b5a641deda16641181aa1e9ee8
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

class Particle:
    def __init__(self,x0,num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS):
        self.position_i=[]          # particle position
        self.velocity_i=[]          # particle velocity
        self.pos_best_i=[]          # best position individual
        self.err_best_i=-1          # best error individual
        self.err_i=-1               # error individual
        self.num_relay = num_relay
        self.CT = CT
        self.Ipick = Ipick
        self.Ifault = Ifault
        self.TOP_desired = TOP_desired
        self.tolerance_dn = tolerance_dn
        self.discrimination_time = discrimination_time
        self.min_TDS = min_TDS
        self.max_TDS = max_TDS
        self.min_TMS = min_TMS
        

        for i in range(0,num_dimensions):
            self.velocity_i.append(np.random.uniform(-1,1))
            self.position_i.append(x0[i])

    # evaluate current fitness
    def evaluate(self,costFunc):
        self.err_i=costFunc(self.position_i,self.num_relay,self.CT, self.Ipick, self.Ifault, self.TOP_desired, self.tolerance_dn, \
                                 self.discrimination_time, self.min_TDS, self.max_TDS, self.min_TMS)

        # check to see if the current position is an individual best
        if self.err_i < self.err_best_i or self.err_best_i==-1:
            self.pos_best_i=self.position_i
            self.err_best_i=self.err_i

    # update new particle velocity
    def update_velocity(self,pos_best_g):
        w=0.95       # constant inertia weight (how much to weigh the previous velocity)
        c1=0.5        # cognative constant
        c2=1.5        # social constant

        for i in range(0,num_dimensions):
            r1=np.random.random()
            r2=np.random.random()

            vel_cognitive=c1*r1*(self.pos_best_i[i]-self.position_i[i])
            vel_social=c2*r2*(pos_best_g[i]-self.position_i[i])
            self.velocity_i[i]=w*self.velocity_i[i]+vel_cognitive+vel_social

    # update the particle position based off new velocity updates
    def update_position(self,bounds):
        for i in range(0,num_dimensions):
            self.position_i[i]=self.position_i[i]+self.velocity_i[i]

            # adjust maximum position if necessary
            if self.position_i[i]>bounds[i][1]:
                self.position_i[i]=bounds[i][1]

            # adjust minimum position if neseccary
            if self.position_i[i] < bounds[i][0]:
                self.position_i[i]=bounds[i][0]
                
class PSO():
    def __init__(self,costFunc,x0,bounds,num_particles,maxiter,num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS):
        global num_dimensions

        num_dimensions=len(x0)
        self.err_best_g=-1                   # best error for group
        self.pos_best_g=[]                   # best position for group

        # establish the swarm
        swarm=[]
        for i in range(0,int(num_particles)):
            swarm.append(Particle([x0[0][i],x0[1][i],x0[2][i],x0[3][i]],num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS))

        # begin optimization loop
        self.best_chart = []
        i=0
        while i < maxiter:
            #print i,err_best_g
            # cycle through particles in swarm and evaluate fitness
            for j in range(0,int(num_particles)):
                swarm[j].evaluate(costFunc)

                # determine if current particle is the best (globally)
                if swarm[j].err_i < self.err_best_g or self.err_best_g == -1:
                    self.pos_best_g=list(swarm[j].position_i)
                    self.err_best_g=float(swarm[j].err_i)

            # cycle through swarm and update velocities and position
            for j in range(0,int(num_particles)):
                swarm[j].update_velocity(self.pos_best_g)
                swarm[j].update_position(bounds)
            self.best_chart.append(self.err_best_g)
            i+=1
        
    def output(self):
        return self.pos_best_g, self.best_chart, self.err_best_g
