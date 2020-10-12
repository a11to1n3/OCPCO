#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:58:50 2020

@author: duypham

Inspired by and inherited from Matlab version of Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb Programmed by Xin-She Yang at Cambridge University Programming dates: Nov 2008 to June 2009 Last revised: Dec 2009 (simplified version for demo only)
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

# A d-dimensional objective function
def Cuckoo_top_functions(P,num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS):
  step_TOP = 0.001
  n = (discrimination_time/step_TOP)*2
  TDS = P[0]
  C1 = P[1]
  C2 = P[2]
  C3 = P[3]
  step_TDS = 0.01
  n_TDS = np.random.randn(int(n),1)*(max_TDS-min_TDS)+min_TDS

  fittemp = np.zeros((int(n),1))
  for i in range(int(n)):
    M=Ifault/Ipick
    K=C1/((M**C2)-1)+C3
    OperatingTimeTemp = K*n_TDS[i]
    if (num_relay==1):
        TOP_min = TOP_desired
    else:
        TOP_min = TOP_desired + min_TMS
    if (TOP_min <= OperatingTimeTemp):
        fittemp[i] = OperatingTimeTemp
    else:
        fittemp[i] = 1000

  fit = fittemp.min()
  index = fittemp.argmin()
  OperatingTime = fit
  C0 = n_TDS[index][0]
  position = np.array([C0,C1, C2, C3])
  return position, OperatingTime

# Application of simple constraints
def simplebounds(s,low,up):
  # Apply the lower bound
  ns_tmp = s
  for i in range(len(low)):
    if ns_tmp[i] < low[i]:
      ns_tmp[i] = low[i]
  
  # Apply the upper bounds 
  for i in range(len(up)):
    if ns_tmp[i] > up[i]:
      ns_tmp[i] = up[i]

  # Update this new move 
  s = ns_tmp
  return s

# Replace some nests by constructing new solutions/nests
def empty_nests(nest,low,up,pa):
  # A fraction of worse nests are discovered with a probability pa
  N=len(nest)
  # Discovered or not -- a status vector
  K = (np.random.rand(len(nest)) > pa).reshape(-1,1)
  # In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
  # this cuckoo's egg is less likely to be discovered, thus the fitness should 
  # be related to the difference in solutions.  Therefore, it is a good idea 
  # to do a random walk in a biased way with some random step sizes.  
  # New solution by biased/selective random walks
  stepsize = np.random.rand()*(nest[np.random.permutation(N)]-nest[np.random.permutation(N)])

  new_nest = nest + stepsize*K
  for j in range(len(new_nest)):
    s = new_nest[j]
    new_nest[j] = simplebounds(s,low,up)  
  return new_nest

def get_best_nest(num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_TDS, max_TDS, min_TMS, nest,newnest,fitness):
  # Evaluating all new solutions
  for j in range(len(nest)):
    fnew,fitnessNew = Cuckoo_top_functions(newnest[j],num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_TDS, max_TDS, min_TMS)

    if fitnessNew <= fitness[j]:
      fitness[j] = fitnessNew
      nest[j]=newnest[j]

  # Find the current best
  fmin = fitness.min()
  K = fitness.argmin() 
  best = nest[K]
  return fmin, best, nest

def get_cuckoos(nest,best,low,up):
  # Levy flights
  N = nest.shape[0]
  # Levy exponent and coefficient
  # For details, see equation (2.21), Page 16 (chapter 2) of the book
  # X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
  beta = 3/2
  sigma = (math.gamma(1+beta)*np.sin(np.pi*beta/2)/(math.gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(1/beta)
  for j in range(N):
    s = nest[j]
    # This is a simple way of implementing Levy flights
    # For standard random walks, use step=1;
    # Levy flights by Mantegna's algorithm
    u = np.random.randn(len(s))*sigma
    v = np.random.randn(len(s))
    step = u/abs(v)**(1/beta)
  
    # In the next equation, the difference factor (s-best) means that 
    # when the solution is the best solution, it remains unchanged.     
    stepsize = 0.01*step*(s-best)
    # Here the factor 0.01 comes from the fact that L/100 should the typical
    # step size of walks/flights where L is the typical lenghtscale; 
    # otherwise, Levy flights may become too aggresive/efficient, 
    # which makes new solutions (even) jump out side of the design domain 
    # (and thus wasting evaluations).
    # Now the actual random walks or flights
    s = s + stepsize*np.random.randn(len(s))
    # Apply simple bounds/limits
    nest[j] = simplebounds(s,low,up)
  return nest

def cuckoo_search(N, max_iter, num_relay, CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_top, max_top, min_C1, max_C1, min_C2, max_C2, min_C3, max_C3, min_TDS, max_TDS,min_TMS):
  # Discovery rate of alien eggs/solutions
  pa=0.25
  
  # Simple bounds of the search domain
  low = [min_TDS, min_C1, min_C2, min_C3]
  up = [max_TDS, max_C1, max_C2, max_C3]
  dim = 4
  # Random initial solutions
  nest = np.zeros((int(N),dim))
  for i in range(int(N)):
    nest[i] = np.array(low)+(np.array(up)-np.array(low))*np.random.rand(len(low))

  # Get the current best
  fitness=10**10*np.ones((int(N),1))
  fmin,bestnest,nest = get_best_nest(num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_TDS, max_TDS, min_TMS,nest,nest,fitness)
  N_it = 0
  best_chart = []
  ## Starting iterations
  for it in range(max_iter):
    # Generate new solutions (but keep the current best)
    new_nest = get_cuckoos(nest,bestnest,low,up)   
    fnew,best,nest = get_best_nest(num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_TDS, max_TDS, min_TMS,nest,new_nest,fitness)
    # Update the counter
    N_it= N_it + int(N) 
    # Discovery and randomization
    new_nest = empty_nests(nest,low,up,pa)
    
    # Evaluate this set of solutions
    fnew,best,nest = get_best_nest(num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_TDS, max_TDS, min_TMS,nest,new_nest,fitness)
    # Update the counter again
    N_it = N_it + int(N)
    # Find the best objective so far  
    if fnew < fmin:
      fmin = fnew
      bestnest = best
    best_chart.append(bestnest)
  return bestnest, best_chart, fmin

def operatingTimeFuction(Ifn,Ip,A,B,C,I50n,TDS,plotI50n):
  if (plotI50n):
    if (Ifn<=I50n):
      operatingTime = ((A/((Ifn/Ip)**B-1))+C)*TDS
    else:
      operatingTime = 0
  else:
      operatingTime = ((A/((Ifn/Ip)**B-1))+C)*TDS
  return operatingTime
