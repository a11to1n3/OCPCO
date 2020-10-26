#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 19:48:49 2020

@author: duypham
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

def FireFly_ffa_mincon(fhandle,u0, Lb=[], Ub=[], para=[20, 500, 0.25, 0.20, 1]):
  n=para[0]
  MaxGeneration=para[1]
  alpha=para[2]
  betamin=para[3]
  gamma=para[4]

  # Total number of function evaluations
  NumEval=n*MaxGeneration

  # Check if the upper bound & lower bound are the same size
  if len(Lb) !=len(Ub):
      display('Simple bounds/limits are improper!')
      return

  # Calcualte dimension
  d = len(u0)

  # Initial values of an array
  zn= np.ones((n,1))*10**100
  # ------------------------------------------------
  # generating the initial locations of n fireflies
  ns,Lightn = FireFly_init_ffa(n,d,Lb,Ub,u0)

  # Iterations or pseudo time marching
  for k in range(MaxGeneration):    # start iterations

    # This line of reducing alpha is optional
    alpha = FireFly_alpha_new(alpha,MaxGeneration)

    # Evaluate new solutions (for all n fireflies)
    for i in range(n):
      zn[i]=fhandle(ns[i,:])
      Lightn[i]=zn[i]

    # Ranking fireflies by their light intensity/objectives
    Lightn = np.sort(zn)
    Index = np.argsort(zn)
    ns_tmp=ns

    for i in range(n):
      ns[i,:]=ns_tmp[Index[i],:]

    # Find the current best
    nso=ns
    Lighto=Lightn
    nbest=ns[0,:]
    Lightbest=Lightn[0];

    # For output only
    fbest=Lightbest

    # Move all fireflies to the better locations
    ns = FireFly_ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,Lightbest,alpha,betamin,gamma,Lb,Ub)

  return nbest,fbest,NumEval

def FireFly_init_ffa(n,d,Lb,Ub,u0):
  # if there are bounds/limits
  if len(Lb)>0:
    ns = np.zeros((n,len(Lb)))
    for i in range(n):
      ns[i,:]=Lb+(Ub-Lb)*np.random.rand(1,d)
  else:
    ns = np.zeros(n,d)
    # generate solutions around the random guess
    for i in range(n):
      ns[i,:]=u0+np.random.rand(1,d)

  # initial value before function evaluations
  Lightn = np.ones((n,1))*10**100
  return ns,Lightn

def FireFly_alpha_new(alpha,NGen):
  delta = 1-(10**(-4)/0.9)**(1/NGen)
  alpha = (1-delta)*alpha
  return alpha

# Move all fireflies toward brighter ones
def FireFly_ffa_move(n,d,ns,Lightn,nso,Lighto,nbest,Lightbest,alpha,betamin,gamma,Lb,Ub):
  # Scaling of the system
  scale = abs(Ub-Lb)

  # Updating fireflies
  for i in range(n):
  # The attractiveness parameter beta=exp(-gamma*r)
    for j in range(n):
      r = np.sqrt(sum((ns[i,:]-ns[j,:])**2))
      # Update moves
      if Lightn[i]>Lighto[j]: # Brighter and more attractive
        beta0 = 1
        beta = (beta0-betamin)*np.exp(-gamma*r**2)+betamin
        tmpf = alpha*(np.random.rand(1,d)-0.5)*scale
        ns[i,:]=ns[i,:]*(1-beta)+nso[j,:]*beta+tmpf

  # Check if the updated solutions/locations are within limits
  ns = FireFly_findlimits(n,ns,Lb,Ub)
  return ns

# Make sure the fireflies are within the bounds/limits
def FireFly_findlimits(n,ns,Lb,Ub):
  for i in range(n):
    # Apply the lower bound
    ns_tmp=ns[i,:]
    I = ns_tmp < Lb
    ns_tmp[I] = Lb[I]

    # Apply the upper bounds
    J = ns_tmp > Ub
    ns_tmp[J] = Ub[J] 
    # Update this new move
    ns[i,:]=ns_tmp

  return ns
