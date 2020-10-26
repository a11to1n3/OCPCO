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

def Ant_ALO(N,Max_iter,lb,ub,dim,fobj):
  antlion_position = Ant_initialization(N,dim,ub,lb)
  ant_position = Ant_initialization(N,dim,ub,lb)

  Sorted_antlions = np.zeros((N,dim))
  Elite_antlion_position = np.zeros((1,dim))
  Elite_antlion_fitness = np.inf
  Convergence_curve = np.zeros((Max_iter))
  antlions_fitness = np.zeros((1,N))
  ants_fitness = np.zeros((1,N))

  # Calculate the fitness of initial antlions and sort them
  for i in range(len(antlion_position)):
    antlions_fitness[0,i]=fobj(antlion_position[i,:]) 

  sorted_antlion_fitness = np.sort(antlions_fitness)
  sorted_indexes = np.argsort(antlions_fitness)
  
  for newindex in range(N):
    Sorted_antlions[newindex,:] = antlion_position[sorted_indexes[0][newindex],:]

  Elite_antlion_position = Sorted_antlions[0]
  Elite_antlion_fitness = sorted_antlion_fitness[0][0]

  Current_iter=1
  while Current_iter < Max_iter:
    # This for loop simulate random walks
    for i in range(len(ant_position)):
      # Select ant lions based on their fitness (the better anlion the higher chance of catching ant)
      Rolette_index = Ant_RouletteWheelSelection(1./sorted_antlion_fitness)
      if Rolette_index == -1:  
        Rolette_index = 1
    
      # RA is the random walk around the selected antlion by rolette wheel
      RA = Ant_Random_walk_around_antlion(dim,Max_iter,lb,ub, Sorted_antlions[Rolette_index,:],Current_iter)
      
      # RA is the random walk around the elite (best antlion so far)
      RE = Ant_Random_walk_around_antlion(dim,Max_iter,lb,ub, Elite_antlion_position[0],Current_iter)

      ant_position[i,:] = RA.T[Current_iter,:]+RE.T[Current_iter,:]/2 # Equation (2.13) in the paper          
    
    for i in range(len(ant_position)):   
        # Boundar checking (bring back the antlions of ants inside search
        # space if they go beyoud the boundaries
        Flag4ub = ant_position[i,:] > np.array(ub).T
        Flag4lb = ant_position[i,:] < np.array(lb).T
        ant_position[i] = (ant_position[i]*(~(Flag4ub+Flag4lb)))+np.array(ub)*Flag4ub+np.array(lb)*Flag4lb
        
        ants_fitness[0,i] = fobj(ant_position[i,:])        

    
    # Update antlion positions and fitnesses based of the ants (if an ant 
    # becomes fitter than an antlion we assume it was cought by the antlion  
    # and the antlion update goes to its position to build the trap)
    double_population = []
    double_population.extend(Sorted_antlions)
    double_population.extend(ant_position)

    double_fitness = []
    double_fitness.extend(sorted_antlion_fitness[0])
    double_fitness.extend(ants_fitness[0])
        
    double_fitness_sorted = np.sort(double_fitness)
    I = np.argsort(double_fitness)
    
    double_sorted_population = np.array(double_population)[I]
        
    antlions_fitness = double_fitness_sorted[:N]
    Sorted_antlions = double_sorted_population[:N,:]
        
    # Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness[0]<Elite_antlion_fitness: 
        Elite_antlion_position = Sorted_antlions[0,:]
        Elite_antlion_fitness = antlions_fitness[0]

      
    # Keep the elite in the population
    Sorted_antlions[0,:] = Elite_antlion_position
    antlions_fitness[0] = Elite_antlion_fitness
  
    # Update the convergence curve
    Convergence_curve[Current_iter] = Elite_antlion_fitness

    # Display the iteration and best optimum obtained so far
    if Current_iter % 500 ==0:
        display('At iteration ', str(Current_iter), ' the elite fitness is ', str(Elite_antlion_fitness))

    Current_iter=Current_iter+1

  return Elite_antlion_fitness, Elite_antlion_position, Convergence_curve

def Ant_initialization(SearchAgents_no,dim,ub,lb):

  Boundary_no = len(ub) # numnber of boundaries
  # If the boundaries of all variables are equal and user enter a signle
  # number for both ub and lb
  if Boundary_no == 1:
    X = np.random.rand(SearchAgents_no,dim)*(np.array(ub)-np.array(lb))+np.array(lb)

  # If each variable has a different lb and ub
  if Boundary_no > 1:
    X = []
    for i in range(dim):
      ub_i=ub[i]
      lb_i=lb[i]
      X.append(np.random.rand(SearchAgents_no,1)*(np.array(ub_i)-np.array(lb_i))+np.array(lb_i))
    X = np.array(X).reshape(-1,dim)
  return X

def Ant_Random_walk_around_antlion(Dim,max_iter,lb, ub,antlion,current_iter):
  if np.array(lb).shape[0] == 1 and np.array(lb).shape[1] == 1: #Check if the bounds are scalar
    lb = np.ones((1,Dim))*np.array(lb)
    ub = np.ones((1,Dim))*np.array(ub)
  if np.array(lb).reshape(-1,1).shape[0] > np.array(lb).reshape(-1,1).shape[1]: # Check if boundary vectors are horizontal or vertical
    lb = np.array(lb).T
    ub = np.array(ub).T

  I=1 # I is the ratio in Equations (2.10) and (2.11)

  if current_iter > max_iter/10:
    I = 1+100*(current_iter/max_iter)

  if current_iter > max_iter/2:
    I = 1+1000*(current_iter/max_iter)

  if current_iter > max_iter*(3/4):
    I = 1+10000*(current_iter/max_iter)

  if current_iter>max_iter*(0.9):
    I = 1+100000*(current_iter/max_iter)
  
  if current_iter>max_iter*(0.95):
    I = 1+1000000*(current_iter/max_iter)

  # Dicrease boundaries to converge towards antlion
  lb=lb/(I) # Equation (2.10) in the paper
  ub=ub/(I) # Equation (2.11) in the paper

  # Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
  if np.random.rand() < 0.5:
      lb = lb + antlion # Equation (2.8) in the paper
  else:
      lb = -lb+antlion

  if np.random.rand() >= 0.5:
      ub = ub+antlion # Equation (2.9) in the paper
  else:
      ub = -ub+antlion

  # This function creates n random walks and normalize accroding to lb and ub
  # vectors 
  RWs = []
  for i in range(Dim):
    X = np.cumsum(2*(np.random.rand(max_iter,1)>0.5)-1).T # Equation (2.1) in the paper
    # [a b]--->[c d]
    a = min(X)
    b = max(X)
    c = lb[i]
    d = ub[i]      
    X_norm =((X-a)*(d-c))/(b-a)+c # Equation (2.7) in the paper
    RWs.append(X_norm)

  return np.array(RWs)

def Ant_RouletteWheelSelection(weights):
  accumulation = np.cumsum(weights)
  p = np.random.rand() * accumulation[-1]
  chosen_index = -1
  for index in range(len(accumulation)):
    if (accumulation[index] > p):
      chosen_index = index
      break

  choice = chosen_index
  return choice
