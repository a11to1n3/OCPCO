#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:58:50 2020

@author: duypham

Inspired by and inherited from Matlab version of Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb Programmed by Xin-She Yang at Cambridge University Programming dates: Nov 2008 to June 2009 Last revised: Dec 2009 (simplified version for demo only)
"""

from CS import *

def Cuckoo_Main(input_file,plot_fitness):
  #Cuckoo_Main Summary of this function goes here
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
    
    if (i==0):
      gBest,best_chart, operating_time = cuckoo_search(N, max_iter, i, CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_top, max_top, min_C1, max_C1, min_C2, max_C2, min_C3, max_C3, min_TDS, max_TDS,min_TMS)
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, 0])
    else:
      gBest,best_chart, operating_time = cuckoo_search(N, max_iter, i, CT, Ipick, Ifault, TOP_actual, tolerance_dn, \
                              discrimination_time, min_top, max_top, min_C1, max_C1, min_C2, max_C2, min_C3, max_C3, min_TDS, max_TDS,min_TMS)
      DTime=abs(output_file[i-1][4]-operating_time)
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, DTime])
    chart.append(best_chart) 
  return output_file, chart
