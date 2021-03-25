#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 23:41:50 2021

@author: ngzicong
"""
from utilities import *

###################################################################################################
#Testing using lecture example on slides 102++
#Edit soil params here
###################################################################################################
surcharge = 0 #kPa
duration = 0.5 #years of which surcharge is applied
length = 0.1 #m, PVD strip long side
width = 5e-3 #m, PVD strip short side
spacing = 1.5 #m c/c spacing
pattern = 'square' #arrangement
smear = 0.8 #smearing factor

#Strata thickness, initial and final effective stresses
#and double drain condition validity
H = 12 #m
Po, Pf = 36,93 #kPa 
double_drain = True #if mudstone/stiff clay base, set to False

#Soil consolidation properties
Pc = Po #kPa, Assumed NC in this case
e0 = 2.16
Cc = 0.9
Cr = 0.06
cv = 1 #m^2/year
ch = 2 #m^2/year

#Graph display #toggle True if needed
Uv_graph = False #Annotated Terzaghi plot
Uh_graph = False #Annotated Barron plot
graph = True #for settlement time plot
add = 2.5 #extra add*duration to observe residual settlement path 


Settlement_Time_1 = PVD_given_duration(length, width, 
                                       spacing, surcharge, 
                                       Po, Pf, Pc, H, 
                                       pattern = pattern, e0 = e0, Cc = Cc, Cr = Cr, cv = cv, ch = ch,
                                       smear = smear, duration = duration, double_drain = double_drain, add = add,
                                       Uv_graph = Uv_graph, Uh_graph = Uh_graph, 
                                       graph = graph, main_color = 'k', legend = False)

Settlement_Time_2 = PVD_given_duration(length, width, 
                                       spacing, 50, 
                                       Po, Pf, Pc, H, 
                                       pattern = pattern, e0 = e0, Cc = Cc, Cr = Cr, cv = cv, ch = ch,
                                       smear = smear, duration = duration, double_drain = double_drain, add = add,
                                       Uv_graph = Uv_graph, Uh_graph = Uh_graph, 
                                       graph = graph, new = False, main_color = 'r', legend = False)

#Settlement_Time_n = [(S_design, S_expected, S_residual, S_surcharge),surcharge_line[0], no_surcharge_line[0], design_line]

S_surcharge = max(Settlement_Time_1[0][-1],Settlement_Time_2[0][-1])
plt.axis([0,(duration+add*duration)*12, S_surcharge*1000,0.])
plt.legend(handles=[Settlement_Time_1[1],Settlement_Time_1[2], 
                    Settlement_Time_2[1], Settlement_Time_2[2],Settlement_Time_2[3], Settlement_Time_2[-1]], 
           loc = 'upper right')
plt.show()