#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 19:48:49 2020

@author: duypham
"""
from Firefly import *

def FireFly_FitnessFunc(x, Ifault, Ipick,nRelay,TOP_Desired,Discrimination_Time, MinTMS,IECType):
  TDS = x[0]
  C1 = x[1]
  C2 = x[2]

  M= Ifault/Ipick
  
  if IECType:
      K=(C1/((M**C2)-1)) # IEC
  else:
      C3 = x[3];
      K=(C1/((M**C2)-1))+C3 # IEEE

  OperatingTime = K*TDS
  
  if nRelay == 1:
      TOP_min = TOP_Desired
  else:
      TOP_min = TOP_Desired + MinTMS

  if TOP_min <= OperatingTime:
      WE = OperatingTime
  else:
      WE = 1000

  return WE

def Firefly_Main(input_file,plot_fitness,IECType,FixABC):
  #Firefly_Main Summary of this function goes here
  #   Detailed explanation goes here
  output_file=[]
  chart = []
  # Size of the swarm " no of objects "
  # parameters [n N_iteration alpha betamin gamma]
  para=[50, 1000, 0.5, 0.2, 1];
  N = 50                        
  max_iter  = 1000         # Maximum number of "iterations"
  num_relay = len(input_file[:,0])
  dim = 4

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

    FixA = input_file[i][19]
    FixB = input_file[i][20]
    FixC = input_file[i][21]

    if FixABC:
      min_C1 = FixA
      min_C2 = FixB
      min_C3 = FixC
      max_C1 = FixA
      max_C2 = FixB
      max_C3 = FixC

    nvars = 5                                         # Number of design variables 
    if IECType:
      dim=3
      LB = np.array([min_TDS, min_C1, min_C2])                 # Lower bound
      UB = np.array([max_TDS, max_C1, max_C2])                 # Upper bound
    else:
      dim=4
      LB = np.array([min_TDS, min_C1, min_C2, max_C3])                 # Lower bound
      UB = np.array([max_TDS, max_C1, max_C2, max_C3])                 # Upper bound

    position0 = LB+(UB-LB)*np.random.rand(1,dim)
    
    if (i==0):
      ObjectiveFunction = lambda x: FireFly_FitnessFunc(x,Ifault, Ipick,i,TOP_desired,discrimination_time,min_TMS,IECType)
      position,fval,NumEval=FireFly_ffa_mincon(ObjectiveFunction,position0,LB,UB,para)
      operating_time = fval
      set_TDS = position[0]
      set_C1 = position[1]
      set_C2 = position[2]
      if IECType:
        set_C3 =0
      else:
        set_C3 = position[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, 0])
    else:
      ObjectiveFunction = lambda x: FireFly_FitnessFunc(x,Ifault, Ipick,i,TOP_desired,discrimination_time,min_TMS,IECType)
      position,fval,NumEval=FireFly_ffa_mincon(ObjectiveFunction,position0,LB,UB,para)
      operating_time = fval
      DTime=abs(output_file[i-1][4]-operating_time)
      set_TDS = position[0]
      set_C1 = position[1]
      set_C2 = position[2]
      if IECType:
        set_C3 =0
      else:
        set_C3 = position[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, DTime])
    chart.append(operating_time) 
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

