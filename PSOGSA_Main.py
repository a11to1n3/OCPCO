"""
@author: duypham

Inspired by and inherited from Matlab version of PSOGSA source code v3.0, Generated by SeyedAli Mirjalili, 2011. Adopted from: S. Mirjalili, S.Z. Mohd Hashim, New Hybrid PSOGSA Algorithm for Function Optimization, in IEEE International Conference on Computer and Information Application?ICCIA 2010), China, 2010, pp.374-377.
"""
from PSOGSA import *

def PSOGSA_Main(input_file,plot_fitness):
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
    
    if (i==0):
      gBestScore,gBest,best_chart, operating_time = PSOGSA(N, max_iter, i, CT, Ipick, Ifault, TOP_desired, tolerance_dn, \
                              discrimination_time, min_top, max_top, min_C1, max_C1, min_C2, max_C2, min_C3, max_C3, min_TDS, max_TDS,min_TMS)
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, 0])
    else:
      gBestScore,gBest,best_chart, operating_time = PSOGSA(N, max_iter, i, CT, Ipick, Ifault, TOP_actual, tolerance_dn, \
                              discrimination_time, min_top, max_top, min_C1, max_C1, min_C2, max_C2, min_C3, max_C3, min_TDS, max_TDS,min_TMS)
      DTime=abs(output_file[i-1][4]-operating_time)
      set_TDS = gBest[0]
      set_C1 = gBest[1]
      set_C2 = gBest[2]
      set_C3 = gBest[3]
      output_file.append([set_TDS, set_C1, set_C2, set_C3, operating_time, DTime])
    chart.append(best_chart) 
  return output_file, chart
