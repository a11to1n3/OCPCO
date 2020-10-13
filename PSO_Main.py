#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 19:48:49 2020

@author: duypham

Inspired by https://gist.github.com/tstreamDOTh/4af1d6b5a641deda16641181aa1e9ee8
"""
from PSO import *

#--- COST FUNCTION 
# function we are attempting to optimize (minimize)
def PSO_top_functions(P,num_relay,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
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

  return fittemp.min()

def PSO_Main(input_file,func):
  #GSA_MAIN Summary of this function goes here
  #   Detailed explanation goes here
  output_file=[]
  chart = []
  # Size of the swarm " no of objects "
  N = 50                        
  max_iter  = 1000         # Maximum number of "iterations"
  num_relay = len(input_file[:,0])

  for i in range(num_relay):
    N = input_file[i,1]
    CT = input_file[i,2]
    Ipick = input_file[i,3]
    Ifault = input_file[i,4]
    
    if i!=0:
      j=i-1
      TOP_actual = output_file[j][4]
    else:
      TOP_desired = input_file[i][5]
    
    tolerance_dn = input_file[i][6]
    discrimination_time = input_file[i][7]
    print(i)
    
    min_top = input_file[i][8]
    max_top = input_file[i][9]
    min_C1 = input_file[i][10]
    max_C1 = input_file[i][11]
    min_C2 = input_file[i][12]
    max_C2 = input_file[i][13]
    min_C3 = input_file[i][14]
    max_C3 = input_file[i][15]
    min_TDS = input_file[i][16]
    max_TDS = input_file[i][17]
    min_TMS = input_file[i][18]
    
    low = [min_TDS, min_C1, min_C2, min_C3]
    up = [max_TDS, max_C1, max_C2, max_C3]
    dim = 4

    initial = []
    bounds = []
    for j in range(dim):
      up_i = up[j]
      low_i = low[j]
      initial.append(np.random.rand(int(N))*(up_i-low_i)+low_i)
      bounds.append((low_i,up_i))
    
    if (i==0):
      pso = PSO(func,initial,bounds,N,max_iter,i,CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS)
      gBest,best_chart, operating_time = pso.output()
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, 0])
    else:
      pso = PSO(func,initial,bounds,N,max_iter,i,CT, Ipick, Ifault, TOP_actual, tolerance_dn, \
                                 discrimination_time, min_TDS, max_TDS, min_TMS)
      gBest,best_chart, operating_time = pso.output()
      DTime=abs(output_file[i-1][4]-operating_time)
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, DTime])
    chart.append(best_chart) 
  return output_file, chart

def operatingTimeFuction(Ifn,Ip,A,B,C,I50n,TDS,plotI50n):
  if (plotI50n):
    if (Ifn<=I50n):
      operatingTime = ((A/((Ifn/Ip)**B-1))+C)*TDS
    else:
      operatingTime = 0
  else:
      operatingTime = ((A/((Ifn/Ip)**B-1))+C)*TDS
  return operatingTime
