# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:01:21 2024
@author: Takanori Harashima 
email : harashima@ims.ac.jp
"""

import itertools
import tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import cv2
from PIL import Image
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
import matplotlib.style as mplstyle
mplstyle.use('fast')
import math
from math import log10, floor
import pandas as pd
import glob
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D   

from tqdm import tqdm
import subprocess
import sys
import random
import matplotlib.animation as animation
import matplotlib.patches as patches
import cv2
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap
import datetime
now = datetime.datetime.now().strftime('%Y%m%d')
import time
import gc
import psutil
start = time.time()
def round_sig(x, sig):
    try:
        return np.around(x, sig-int(floor(log10(abs(x))))-1)
    except:
        return x
def MEMORY_RELEASE():    
    gc.collect()
    mem_using = psutil.virtual_memory().percent
  #  print('Using '+ str(mem_using) + '% memory')
    if mem_using > 90:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('!!!!!!!!!!!!!       Memory error     !!!!!!!!!!!!!!!!!')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        gc.collect()
        print(memory_error_memory_error_memory_error_memory_error_memory_error)
def folder_maker(folpath):
    if os.path.isdir(folpath) == False:
        os.mkdir(folpath)
from matplotlib.collections import LineCollection
def COLOR_LINE_PLOT(x,y,c,AX):
  #  x    = np.linspace(0,1, 100)
 #   y    = np.linspace(0,1, 100)
  #  cols = np.linspace(0,1,len(x))    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='jet')
    lc.set_array(c)
    lc.set_linewidth(1.0)
    line = AX.add_collection(lc)
    return line


plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["image.interpolation"] = "none"
plt.rcParams['font.family'] = 'Arial' 
plt.rcParams['font.size'] = 16 
plt.rcParams['axes.titlesize'] = 16 
plt.rcParams['lines.linewidth'] = 1 
#plt.rcParams['axes.xmargin'] = 0
#plt.rcParams['axes.ymargin'] = 0


###############################################################################

####################################################################################
# Geometric parameters (Fixed parameters)
particle_size = 100                    # Diameter of AuNP [nm]
r_particle = particle_size / 2         # Radius of AuNP [nm]
RNA_density =0.022                     # Surface density of RNAs on substrate[molecules/nm2]
DNA_density =0.10                      # Surface density of DNAs on AuNP [molecules/nm2]
##############################################################################
def END_TO_END_DISTANCE(lc,lp):
    return np.sqrt(2*lp**2 * (np.exp(-lc/lp) -1 + lc/lp))
# DNA leg end-to-end length [nm]
Lc_base = 0.34 # nm/base, Stacking displacement of DNAs and RNAs
N_base = 45 # bases
Lp_ssDNA = 2.0 # nm
Lc_ssDNA = Lc_base * N_base
L_DNA = END_TO_END_DISTANCE(Lc_ssDNA,Lp_ssDNA)
##############################################################################
# Accessible radius [nm]
path_radius = np.sqrt( (L_DNA+r_particle)**2 - r_particle**2  )
##############################################################################
# Mobile region of DNA/RNA duplex [nm]                
# Assuming only the rotation of AuNP
# Mobile region of the particle centroid when DNA/RNA duplex can be elongated to its contour length.# DNA leg + RNADNA length calculation ########################################
def cycroid(r,theta):
    x = r*(theta - np.sin(theta))
    y = r*(1 - np.cos(theta))
    R = np.sqrt(x**2 + y**2)
    return x,y,R
Lc_base     = 0.34 # nm/base, Stacking displacement of DNAs and RNAs
Lp_basepair = 50  # nm, Persistence length of dsDNA and DNA/RNA hybrid
Lp_base     = 2.0 # nm, Persistence length of ssDNA and ssRNA
# ssDNA15 base + dsDNA 15 base +ssRNA 5 base + dsRNA/DNA15 base + ssDNA 30 base
L_ssDNA15base  = END_TO_END_DISTANCE(Lc_base*15,Lp_base)
L_dsDNA15base  = END_TO_END_DISTANCE(Lc_base*15,Lp_basepair)
L_ssRNA5base   = END_TO_END_DISTANCE(Lc_base*5,Lp_base)
L_RNADNA15base = END_TO_END_DISTANCE(Lc_base*15,Lp_basepair)
L_ssDNA30base  = END_TO_END_DISTANCE(Lc_base*30,Lp_base)
L0_total = sum([L_ssDNA15base,L_dsDNA15base,L_ssRNA5base,L_RNADNA15base,L_ssDNA30base])
Lc_total = Lc_base*(15+15+5+15+30)
dLmax = Lc_total - L0_total
theta = np.linspace(0,2*np.pi,10000)
x,y,R = cycroid(r_particle,theta)
cx = r_particle*theta
R_mobile_nm = cx[np.where(R<dLmax)[0][-1]]
##############################################################################
##############################################################################                                       
####################################################################################
####################################################################################
# simulation parameter calculation
path_radius = ((L_DNA + r_particle)**2 - r_particle**2 )**0.5 
# Averaged number of accessible RNA molecules  
path_RNA_molecules = RNA_density * (np.pi*path_radius **2)
# Accessible AuNP surface area
path_DNA_radius_nm2 = (4*np.pi*r_particle**2) * (np.degrees(np.arcsin(path_radius/(L_DNA + r_particle)))*2/360)
# Averaged number of accessible DNA molecules   = Accessible AuNP surface area x DNA surface density
path_DNA_molecules = DNA_density * path_DNA_radius_nm2
# Pixel size = 1/DNA surface density
pixel_to_nm = np.sqrt(1/DNA_density) # nm
# mobile region [pixel]
R_mobile_pixel = R_mobile_nm / pixel_to_nm
# Accessible radius [pixel]
r_path = path_radius / pixel_to_nm  
DNA_size = round(r_path) 

print ("path_DNA_molecules",path_DNA_molecules)
print ("path_RNA_molecules",path_RNA_molecules)
print ('R_mobile_pixel', R_mobile_pixel)
print ('R_mobile nm', R_mobile_nm)

print('path_radius (nm)',path_radius)
print ('r_path nm', r_path *pixel_to_nm)
print ('r_path px', r_path )
print ("DNA_size",DNA_size)

###############################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################# variable parameters ############################# 
############################# variable parameters ############################# 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################
# Basic parameters
globalfol = r'Directory\to\perform\DNAmotor\simulation'
date = 'DATE'
RNA_size = 2000                        # Full range of RNA substrate (px)  (default:2000)
tmax_simu = 100000                     # Uplimit for simulation (sec)  (default:100000)
N_simu = 3                             # Nunber of trajectory generate (default:5)
frame_per_event = 1000                 # Span frame to check the progress of simulation (default:1000)
foli=0                                 # ID of the condition of kinetic parameters (default:0)
# Kinetic parameters
khyb_list = np.array([0.2])                  # DNA/RNA Hybridization rate [s-1]
konE_list = np.array([1.0]) *10**6           # RNase H binding rate [M-1 s-1]
kcat_list = np.array([3])                    # RNA hydrolysis rate [s-1] 
RNaseH_list = np.array([36])                 # RNase H condition [nM]
#######################
SIMULATION = True # Run the new simulation
SUMMERIZE = True # Output rough summary of the simulation
#!!!!!!!!!!!!!!!!
#########################################################################################################################################################################################
#########################################################################################################################################################################################
colortable_RNA = ListedColormap(['0.7','g','r','w']) # ['Unbound', 'Bound', 'Hydrolyzed', 'Empty']


def CENTER_TRANSLOCATE(cx, cy, DNA_state, cx_list, cy_list ,i ):
    
    xmin_disp, xmax_disp= min([cx-display_size,min(cx_list) ]), max([cx+display_size,max(cx_list) ])
    ymin_disp, ymax_disp= min([cy-display_size,min(cy_list) ]), max([cy+display_size,max(cy_list) ])
    global min_disp, max_disp
    min_disp, max_disp = min([xmin_disp,ymin_disp, min_disp]), max([xmax_disp,ymax_disp, max_disp])
    #################################
    
    #################################
    global min_disp_nm, max_disp_nm
    min_disp_nm, max_disp_nm = (min_disp-RNA_size//2)*pixel_to_nm + RNA_size//2, (max_disp-RNA_size//2)*pixel_to_nm + RNA_size//2
    #################################
    
    #   condition_move = DNA_state[cx-size:cx+size,cy-size:cy+size] == 1
    condition_move = DNA_state[min_disp:max_disp,min_disp:max_disp] == 1
    condition_move_cut = DNA_state[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ] == 1
    x_cut, y_cut = np.where(condition_move_cut)
    x, y = np.where(condition_move)
    
    img_sum = np.zeros((max_disp - min_disp,max_disp - min_disp)) 
    """
    plt.imshow(DNA_state[min_disp:max_disp,min_disp:max_disp])
    plt.show()
    plt.close()
    """    
    if len(x) == 0:
        print ('!!!!!!!!!!!!!!!!!!  DNA_state all zero  !!!!!!!!!!!!!!!!!!!!!')
        
    img_calc = np.zeros((DNA_size*6, DNA_size*6)) 
    for xi,yi in zip(x_cut, y_cut):   
        img = np.zeros((DNA_size*6, DNA_size*6), np.uint8)      
        cv2.circle(img, (xi ,yi) , round(R_mobile_pixel) , (1,1,1), thickness=-1)
        img_calc += (img).astype(np.float32) 
    img_calc = img_calc/len(x_cut)
    # y0:y1,x0:x1
    img_sum[(cy-min_disp) -DNA_size*3 :(cy-min_disp) + DNA_size*3 ,(cx-min_disp) -DNA_size*3 :(cx-min_disp)+ DNA_size*3 ] = img_calc
  #  img_sum = img_sum/len(x)

    """
    plt.close()
    plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.imshow(img_calc/len(x), clim=(0,1),cmap=cm.gist_ncar)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(img_sum, clim=(0,1),cmap=cm.gist_ncar)
    plt.colorbar()
    plt.scatter(x,y,color = 'k')
    plt.show(), plt.close()
    """

    def PLOT_MOBILE(AX6,FIX_OR_MOVE):
        
        AX6.set_title("Blue=before, Orange=after movement\n")
        cmap = cm.gray
        cmap.set_over('r')
        m6 = AX6.imshow(img_sum, extent=(min_disp,max_disp,min_disp,max_disp),clim=(0,0.999),cmap = cmap ,origin = 'lower')  
        cbar = fig.colorbar(m6, ax=AX6,label='Fraction_overlap\nred:Mobile_region', shrink=0.5, orientation="horizontal")
        c = patches.Circle(xy=(cx, cy), radius=r_path, ec='b',fill=False,lw=2)
        AX6.add_patch(c)
        
        AX6.scatter(x+ min_disp+0.5,y+ min_disp+0.5,marker = '+',color = 'c',s= 30)
        AX6.plot([cx+0.5],[cy+0.5],marker = '+',color = 'b',ms= 30)
        if FIX_OR_MOVE == 'move':
            AX6.plot([cx_new+0.5],[cy_new+0.5],marker = '+',color = 'orange',ms= 30)
            c_new = patches.Circle(xy=(cx_new, cy_new), radius=r_path, ec='orange',fill=False,lw=2)
            AX6.add_patch(c_new)
        AX6.set(xlabel='X (pixel)',ylabel='Y (pixel)',xlim = (min_disp,max_disp), ylim = (min_disp,max_disp) )
  
    
    # No movement occurs unless the mobile regions of all bound sites overlap.
    if len(np.where(img_calc == 1)[0]) == 0:
        if (i % frame_per_event == 0) & (i != 0):
            global fig, AX1, AX2, AX3, AX4, AX5, AX6
            fig = plt.figure(figsize=(18, 14))
            AX1, AX2, AX3, AX4, AX5, AX6 = fig.add_subplot(231), fig.add_subplot(232), fig.add_subplot(234), fig.add_subplot(235), fig.add_subplot(233), fig.add_subplot(236)
            PLOT_MOBILE(AX6,'fix')
        return cx, cy, 0 , 0


    y_id, x_id = np.where(img_calc == 1) 
    index = random.randint(0, len(x_id)-1)
    cx_new, cy_new =  x_id[index] + cx - DNA_size*3, y_id[index] + cy - DNA_size*3
  #  y_id = np.where(img_sum == 1)[0] 
    
    if (i % frame_per_event == 0) & (i != 0):
        fig = plt.figure(figsize=(18, 14))
        AX1, AX2, AX3, AX4, AX5, AX6 = fig.add_subplot(231), fig.add_subplot(232), fig.add_subplot(234), fig.add_subplot(235), fig.add_subplot(233), fig.add_subplot(236)
        PLOT_MOBILE(AX6,'move')
    
  #  print (cx_new, cy_new, cx_new -cx , cy_new -cy)
    return cx_new, cy_new, cx_new -cx , cy_new -cy





t01_sum,t12_sum,t23_sum,t34_sum,t45_sum,t56_sum,t67_sum = 0,0,0,0,0,0,0
t2_01_sum,t2_12_sum,t2_23_sum = 0,0,0

def SIMULATE_TRAJ(i_traj,trace_fol, khyb,konE,kcat):
    # Initialization
    RNA_state = np.zeros(RNA_size*RNA_size, dtype = 'int32')
    DNA_state = np.zeros((RNA_size, RNA_size))
    event_prob = np.zeros((RNA_size, RNA_size))
    event_type = np.zeros((RNA_size, RNA_size))

    # distribute RNAs
    n_vacant = int(RNA_size *RNA_size* (1-RNA_density/DNA_density)) 
    print(RNA_density/DNA_density)
    n_vacant_i = random.sample(range(RNA_size *RNA_size), k=n_vacant)
    RNA_state[np.array(n_vacant_i)] = 3
    RNA_state = RNA_state.reshape(RNA_size,RNA_size)

    ########################################################################################333
    Xi = np.ravel([np.ones(DNA_size*6+1)*i for i in range(DNA_size*6+1)]) 
    Yi = np.ravel([np.arange(DNA_size*6+1) for i in range(DNA_size*6+1)])
    df_xiyi = pd.DataFrame({'Xi':np.array(Xi, dtype= 'int32'),'Yi':Yi})
    
    Ri = ( (np.array(df_xiyi['Xi']) - DNA_size*3  )**2 +(np.array(df_xiyi['Yi']) - DNA_size*3  )**2 ) **0.5
    Ri_matrix = Ri.reshape([DNA_size*6+1,DNA_size*6+1])
    
 #   condition_out_circle_cut = Ri_matrix >= DNA_size 
#    condition_in_circle_cut = Ri_matrix < DNA_size
    condition_out_circle_cut = Ri_matrix > r_path
    condition_in_circle_cut = Ri_matrix <= r_path
    condition_in_circle_cut_01 = np.where(condition_in_circle_cut, 1, 0)
    
    """
    plt.imshow(Ri_matrix*condition_in_circle_cut_01,cmap = 'gist_ncar_r')
    plt.colorbar(), plt.show(), plt.close()
    """
    ########################################################################################333
    
    
    cx, cy = int(RNA_size/2), int(RNA_size/2)

    # v3.09 Simulation starts from single bound site
    event_x,event_y = 0,0
    event = 1
    lx = np.arange(cx - DNA_size, cx + DNA_size)
    ly = np.arange(cy - DNA_size, cy + DNA_size)
    istart = True
    for xi in lx:
        for yi in ly:
            d = ((xi-cx)**2+(yi-cy)**2)**0.5
            if (d <= r_path) & (RNA_state[xi, yi] != 3):
                if istart == True:
                    RNA_state[xi, yi] = 1
                    DNA_state[xi, yi] = 1
                    event_type[xi, yi] = 1.5
                    print(xi,yi)
                    event_x,event_y = xi,yi
                    ##############
                    ##############
                    istart = False
                    ##############
                    ##############
    RNA_state = np.array(RNA_state, dtype = 'int32')
    DNA_state = np.array(DNA_state, dtype = 'int32')
    i = 0
    t = 0
    dx,dy = 0,0

    cx_list, cy_list, dt_list, dx_list, dy_list = [], [], [], [], []
    event_popu = {}
    
    global min_disp, max_disp,display_size
    display_size = DNA_size*3
    
    min_disp, max_disp =  cx-display_size, cx+display_size+1
    
    while True:
        t0 = time.time()
        
        def STATE_SHOW(DNA_state, RNA_state, event_type, event_prob):
            print ("sample "+str(i_traj))
            t_now = time.time()- start
            print ("run_time "+str(t_now) +' sec')
            print ('nokori_time_to_complete_this_sample :' + str(round(t_now/(t+ 0.01)*tmax_simu/60, 1)) + ' min' )
            print ('i='+str(i), ' t='+str(t))            
         #   print ('enzyme = '+str(enzyme_count), ' / '+str(enzyme_max))
            # plt.subplot(221)
            t6_0  = time.time() 

            for AXi in [AX6] :                      
                AXi.text(min_disp+(max_disp-min_disp)*0.65 ,max_disp-(max_disp-min_disp)*0.85, 'Event = '+str(i)+'\n'+'t = '+ str(round_sig(t,3)) + ' s' ,color = 'w')

            for AXi in [AX2, AX5]:                      
                AXi.text(min_disp_nm+(max_disp_nm-min_disp_nm)*0.65 ,max_disp_nm-(max_disp_nm-min_disp_nm)*0.85, 'Event = '+str(i)+'\n'+'t = '+ str(round_sig(t,3)) + ' s' ,color = 'k')


            """
            AX1.set_title( "Sample_" +('00000'+str(i_traj))[-4:])
            m1 = AX1.imshow(DNA_state[min_disp:max_disp,min_disp:max_disp].T, clim=(0,1),cmap = ListedColormap(['navy','w'])\
                       ,extent = (min_disp,max_disp,min_disp,max_disp), origin='lower')   
            AX1.set_xlabel('x / pixel'),AX1.set_ylabel('y / pixel')            
            cbar = fig.colorbar(m1, ax=AX1,label='DNA state',ticks=[0.25, 0.75], shrink=0.7, orientation="horizontal")
            cbar.ax.set_xticklabels(['Unbound', 'Bound'])
            """
            #######################################################################################################################
            #######################################################################################################################
            cx_nm = (cx-RNA_size//2)*pixel_to_nm + RNA_size//2
            cy_nm = (cy-RNA_size//2)*pixel_to_nm + RNA_size//2
            AX1.set_title("p" +('00000'+str(i_traj))[-3:]+ ', [RNaseH]='+ str(RNaseH) + 'nM'+'\n')
            zoom_size = 20
            zoom_size_nm = zoom_size * pixel_to_nm
            AX1.text(cx_nm+zoom_size_nm*0.30 ,cy_nm-zoom_size_nm*0.7, 'Event = '+str(i)+'\n'+'t = '+ str(round_sig(t,3)) + ' s' ,color = 'k', ha='left')
            scalebar = 20 # nm
            AX1.plot( [cx_nm+zoom_size_nm*0.5,cx_nm+zoom_size_nm*0.5 + scalebar],[cy_nm+zoom_size_nm*0.7,cy_nm+zoom_size_nm*0.7], lw = 8,c='k')
            AX1.text(cx_nm+zoom_size_nm*0.5+ scalebar/2 ,cy_nm+zoom_size_nm*0.81, str(scalebar)+' nm', va='center', ha='center',weight='bold',fontsize= 20,zorder=10)

        #    m1 = AX1.imshow(RNA_state[cy-zoom_size:cy+zoom_size,cx-zoom_size:cx+zoom_size].T, clim=(0,4),cmap=colortable_RNA\
        #               ,extent = (cy-zoom_size,cy+zoom_size,cx-zoom_size,cx+zoom_size), origin='lower')    
    #        m1 = AX1.imshow(RNA_state[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, clim=(0,4),cmap=colortable_RNA\
            m1 = AX1.imshow(RNA_state[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, clim=(0,4),cmap=colortable_RNA\
                       ,extent = ( cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm, cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm), origin='lower')    
     #                  ,extent = (cx-zoom_size*pixel_to_nm,cx+zoom_size*pixel_to_nm,cy-zoom_size*pixel_to_nm,cy+zoom_size*pixel_to_nm), origin='lower')    
            AX1.set(xlabel= 'X (nm)',ylabel= 'Y (nm)',xlim=(cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm),ylim=(cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm) )
            cbar = fig.colorbar(m1, ax=AX1,label='RNA state',ticks=[0.5,1.5,2.5,3.5], shrink=0.7, orientation="horizontal")
            cbar.ax.set_xticklabels(['Unbound', 'Bound', 'Hydrolyzed', 'Empty'])
            #######################################################################################################################
            #######################################################################################################################

            AX2.set_title('Particle_size:' + str(particle_size) +'\n')
            m2 = AX2.imshow(RNA_state[min_disp:max_disp,min_disp:max_disp].T, clim=(0,4),cmap=colortable_RNA\
                       ,extent = ( min_disp_nm, max_disp_nm, min_disp_nm, max_disp_nm), origin='lower')    
#                       ,extent = (min_disp,max_disp,min_disp,max_disp), origin='lower')    
            AX2.set(xlabel= 'X (nm)',ylabel= 'Y (nm)')
            cbar = fig.colorbar(m2, ax=AX2,label='RNA state',ticks=[0.5,1.5,2.5,3.5], shrink=0.7, orientation="horizontal")
            cbar.ax.set_xticklabels(['Unbound', 'Bound', 'Hydrolyzed', 'Empty'])
            
            AX3.set_title('event_type')
     #       m3 = AX3.imshow(event_type[min_disp:max_disp,min_disp:max_disp].T, clim=(0,5),cmap=ListedColormap(['darkred','navy','g','m','darkorange'])\
            m3 = AX3.imshow(event_type[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, clim=(0,5),cmap=ListedColormap(['darkred','navy','g','m','darkorange'])\
                       ,extent = ( cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm, cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm), origin='lower')    
            AX3.set(xlabel= 'X (nm)',ylabel= 'Y (nm)',xlim=(cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm),ylim=(cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm) )
            cbar = fig.colorbar(m3, ax=AX3,ticks=[0.5,1.5,2.5,3.5,4.5], shrink=0.7, orientation="horizontal")
            cbar.ax.set_xticklabels(['Empty\nor\nInnaccessible','k$_{on}^{DNA}$','k$_{on}^E$', 'k$_{cat}^E$', 'Hydrolyzed'])
            AX4.set_title('event_prob')
            
         #   m4 = AX4.imshow(event_prob[min_disp:max_disp,min_disp:max_disp].T, norm=colors.LogNorm(vmin=10**-3,vmax=10),cmap = 'jet'\
            m4 = AX4.imshow(event_prob[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, norm=colors.LogNorm(vmin=10**-3,vmax=10),cmap = 'jet'\
                       ,extent = ( cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm, cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm), origin='lower')    
            AX4.set(xlabel= 'X (nm)',ylabel= 'Y (nm)',xlim=(cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm),ylim=(cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm) )
            fig.colorbar(m4, ax=AX4,label='Time / s' , shrink=0.7, orientation="horizontal")
            
            """
            m = AX5.imshow()    
            fig.colorbar(m ,ax=AX5,shrink=0.7, orientation="horizontal")
            """
            cbar = fig.colorbar(cm.ScalarMappable(norm=None, cmap=ListedColormap(['w','w'])), ax=AX5,shrink=0.7, orientation="horizontal", drawedges=True)
            cbar.outline.set_color('w')            
            cbar.dividers.set_color('w')
            cbar.ax.axis('off')
            
            AX5.set_title("(cx,cy)="+str((cx-RNA_size//2, cy-RNA_size//2))  + ' px'+ '\n' )
         #   AX5.plot([int(RNA_size/2), int(RNA_size/2) ],[min_disp, max_disp ], lw = 1,color = 'gray')
         #   AX5.plot([min_disp, max_disp ], [int(RNA_size/2), int(RNA_size/2)], lw = 1,color = 'gray')
            for AX in [AX1,AX2,AX5]: 
                AX.axhline(RNA_size//2, lw = 1,color = 'gray',zorder=1)
                AX.axvline(RNA_size//2, lw = 1,color = 'gray',zorder=1)
            cx_nm_list = (np.array(cx_list) - RNA_size//2 ) * pixel_to_nm + RNA_size//2
            cy_nm_list = (np.array(cy_list) - RNA_size//2 ) * pixel_to_nm + RNA_size//2
            AX5.plot(cx_nm_list,cy_nm_list, lw = 1,color = 'k',alpha = 0.2,marker= 'o',ms = 3)
            AX5.plot(cx_nm,cy_nm, lw = 1,color = 'y',marker= '+',ms = 10)
            AX5.set_xlim( min_disp_nm, max_disp_nm)
            AX5.set_ylim( min_disp_nm, max_disp_nm)
            AX5.set_xlabel('X (nm)'),AX5.set_ylabel('Y (nm)')
            AX5.set_aspect('equal')
            
            
            t6_1  = time.time() 
            t6_01 = t6_1 - t6_0

            
            c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
   #         c = patches.Circle(xy=(cx, cy), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
            AX1.add_patch(c)
            c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
            AX2.add_patch(c)
            c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
            AX3.add_patch(c)
            c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
            AX4.add_patch(c)
            c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
            AX5.add_patch(c)
            
            CENTERING=True
            plt.tight_layout()
            fig.canvas.draw()
            
            if CENTERING:
                for AXi in [AX1,AX2, AX3, AX4, AX5, AX6]:
                   # print(AXi.get_xticklabels())
                    xlabels = [ str( int(item.get_text().replace('−','-')) - int(RNA_size/2) ).replace('-','−') for item in AXi.get_xticklabels() \
                               if item.get_text() != '']
                    AXi.set_xticklabels(xlabels)
                    ylabels = [ str( int(item.get_text().replace('−','-')) - int(RNA_size/2) ).replace('-','−') for item in AXi.get_yticklabels() \
                               if item.get_text() != '']
                    AXi.set_yticklabels(ylabels)
                    AXi.invert_yaxis()
            t6_2  = time.time() 
            t6_12 = t6_2 - t6_1

            plt.savefig(trace_fol + os.sep + ('0000'+str(i_traj))[-3:]+'_'+ ('00000'+str(i))[-6:] +'.png' ,dpi = 50)

            if i%1000 == 0:
                print(i)
                #plt.show()
                MEMORY_RELEASE()
            plt.clf()
            plt.close()
            
            t6_3  = time.time() 
            t6_23 = t6_3 - t6_2     
            
            print ("t6_01 "+str(t6_01) +' sec')
            print ("t6_12 "+str(t6_12) +' sec')
            print ("t6_23 "+str(t6_23) +' sec')

            
        global t01_sum,t12_sum,t23_sum,t34_sum,t45_sum,t56_sum,t67_sum
        t1  = time.time() 
        t01 = t1 - t0
        t01_sum += t01
        
      #  print(i)

        # no event happen at the first iteration
        if i != 0:
            global t2_01_sum,t2_12_sum,t2_23_sum
            t2_0  = time.time() 

            event_prob_cut = event_prob[cx-DNA_size*2:cx+DNA_size*2+1: ,cy-DNA_size*2:cy+DNA_size*2+1 ]
            event_type_cut = event_type[cx-DNA_size*2:cx+DNA_size*2+1: ,cy-DNA_size*2:cy+DNA_size*2+1 ]
            
            #########################################
            ########## Gillespie Algorithm ##########
            # Calculate time step
            r_tot = khyb * len(np.where(event_type_cut==1)[0]) \
                  + konE * len(np.where(event_type_cut==2)[0])\
                  + kcat * len(np.where(event_type_cut==3)[0])
            delta_t = random.expovariate(r_tot)
            
            # Make reaction probability list
            w = np.ravel(event_prob_cut)
            w = w / sum(w)
            dice = np.arange(len(w))
            # Choice reaction site
            react_choice = random.choices(dice, weights = w)[0]
            x_choice = react_choice//event_prob_cut.shape[1]
            y_choice = react_choice%event_prob_cut.shape[1]
            event = event_type_cut[x_choice,y_choice]
            #########################################
            #########################################
                        
            event_x, event_y = cx -DNA_size*2 + x_choice, cy-DNA_size*2 + y_choice
            if event == 1:
                # DNA/RNA hybridization
                RNA_state[event_x, event_y] = 1
                DNA_state[event_x, event_y] = 1
                event_type[event_x, event_y] = 1.5 # Waiting RNase H binding 

            elif event == 2:
                # RnaseH binding
#                print ('RnaseH bind!')
            #    print (event_time[event_position])
                event_type[event_x, event_y] = 2.5
            elif event == 3:
                # RNA hydrolysis
             #   print ('hydrolyzed ')
            #    print (event_time[event_position])
                RNA_state[event_x, event_y] = 2
                DNA_state[event_x, event_y] = 0
                event_prob[event_x, event_y] = 0
                event_type[event_x, event_y] = 4
            
          

            """
            print(event)            
            event_position = np.where(event_prob == delta_t)
            event = event_type[event_position]
            print(event_position)
            print(event)
            print(eventdsafdsafdsaf)
            """
            
            # event 0 ： inert
            # event 1 : Waiting RNA/DNA hybridization
            # event 2 : Waiting RNase H binding 
            # event 3 : Waiting RNA hydrolysis
    
            t2_1  = time.time() 
            t2_01 = t2_1 - t2_0
            t2_01_sum += t2_01
        
    
            cx_list.append(cx)
            cy_list.append(cy)
            dt_list.append(delta_t)
            dx_list.append(dx)
            dy_list.append(dy)

            t2_2  = time.time() 
            t2_12 = t2_2 - t2_1
            t2_12_sum += t2_12

            cx, cy, dx, dy = CENTER_TRANSLOCATE(cx, cy, DNA_state,cx_list, cy_list,i )               

            t2_3  = time.time() 
            t2_23 = t2_3 - t2_2
            t2_23_sum += t2_23

            
            t += delta_t

        """
        plt.imshow((Ri_matrix/DNA_size)**1)
        plt.show()
        plt.close()
        print (afdaf)
        """
        
        t2  = time.time() 
        t12 = t2 - t1
        t12_sum += t12

        t3  = time.time() 
        t23 = t3 - t2
        t23_sum += t23
        

      #  condition_out_circle_cut = condition_out_circle[cx-DNA_size*3:cx+DNA_size*3 ,cy-DNA_size*3:cy+DNA_size*3 ]
    #    print(DNA_state[cx - DNA_size:cx + DNA_size+1, cy - DNA_size:cy + DNA_size+1]==0)
        
        event_prob[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ][condition_out_circle_cut] = 0
        event_prob[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ][condition_out_circle_cut] = 0
        event_type[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ][condition_out_circle_cut] = 0
        DNA_state[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ][condition_out_circle_cut]  = 0
        RNA_state[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ][(condition_out_circle_cut) & (RNA_state[cx-DNA_size*3:cx+DNA_size*3+1 ,cy-DNA_size*3:cy+DNA_size*3+1 ] != 2) & (RNA_state[cx-DNA_size*3:cx+DNA_size*3+1 ,cy-DNA_size*3:cy+DNA_size*3+1 ] != 3) ] = 0
        
        """
        plt.close()
        plt.title('2')
        plt.imshow(event_prob[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ], norm=colors.LogNorm(vmin=10**-3,vmax=1000),cmap = 'jet')
    #    plt.colorbar()
        plt.show()
        plt.close()
        """

        """
        event_prob[condition_out_circle] = 10000
        event_type[condition_out_circle] = 0
        DNA_state[condition_out_circle] = 0
        RNA_state[(condition_out_circle) & (RNA_state != 2) & (RNA_state != 3) ] = 0
        """
        ##################################################################################################################
        t4  = time.time() 
        t34 = t4 - t3
        t34_sum += t34
        ###########################################################################################################
        ###########################################################################################################
        ###########################################################################################################
        # Limit the area when we consider the reaction
        event_type_cut = event_type[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ]
        event_prob_cut = event_prob[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ]
        DNA_state_cut = DNA_state[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ]
        RNA_state_cut = RNA_state[cx-DNA_size*3:cx+DNA_size*3+1: ,cy-DNA_size*3:cy+DNA_size*3+1 ]
        
        
   #     condition_in_circle_cut = condition_in_circle[cx-DNA_size*3:cx+DNA_size*3: ,cy-DNA_size*3:cy+DNA_size*3 ]
        condition_reactive_cut = np.logical_not( ((DNA_state_cut == 0) & (RNA_state_cut == 2)) | ((DNA_state_cut == 0) & (RNA_state_cut == 3)) )
        condition_in_circle_reactive_cut = np.where( (condition_in_circle_cut) & (condition_reactive_cut) )

        """
        plt.close()
        plt.title('(DNA_state_cut == 0) & (RNA_state_cut == 0)')
        plt.imshow( (DNA_state_cut == 0) & (RNA_state_cut == 0),cmap = 'jet')
     #   plt.colorbar()
        plt.show()
        plt.close()
        """

        t5  = time.time() 
        t45 = t5 - t4
        t45_sum += t45

        #######################################################################
        #######################################################################
        # shuffle ID of reaction site
        C = np.copy(condition_in_circle_reactive_cut).T
        np.random.shuffle(C)        
        for xy in C:
            xi,yi= xy[0],xy[1]
            xi_org, yi_org =  xy[0]+cx-DNA_size*3 ,xy[1]+cy-DNA_size*3
        
            # {unbound RNA present} but {khyb did not start} → start khyb 
            if (event_type_cut[xi, yi] == 0) & (RNA_state_cut[xi, yi] == 0):
                event_prob[xi_org, yi_org] = khyb
                event_type[xi_org, yi_org] = 1
            #before enzyme binding → start konE
            if(event_type_cut[xi, yi] == 1.5):
                event_prob[xi_org, yi_org] = konE
                event_type[xi_org, yi_org] = 2
            #Enzyme binding → start kcat
            if (event_type_cut[xi, yi] == 2.5):
                event_prob[xi_org, yi_org] = kcat
                event_type[xi_org, yi_org] = 3
            #Hydrolysis completed → RNA state = 4
            if (event_type_cut[xi, yi] == 4):
                event_prob[xi_org, yi_org] = 0
                event_type[xi_org, yi_org] = 4
        #######################################################################
        #######################################################################
        
        t6  = time.time() 
        t56 = t6 - t5
        t56_sum += t56
        
        
        
        if (i % frame_per_event == 0) & (i != 0):            
            STATE_SHOW(DNA_state, RNA_state, event_type, event_prob)
        
                
        if (np.all(DNA_state[cx - DNA_size:cx + DNA_size+1, cy - DNA_size:cy + DNA_size+1] == 0)):
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   detachment   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(i)
            _ = CENTER_TRANSLOCATE(cx, cy, DNA_state,cx_list, cy_list,frame_per_event )   
            STATE_SHOW(DNA_state, RNA_state, event_type, event_prob)
            break

        if t > tmax_simu:
            STATE_SHOW(DNA_state, RNA_state, event_type, event_prob)
            print("Simulation time is over preset tmax")
            break
        print(event)

        event_type_cut = event_type[cx - DNA_size:cx + DNA_size+1, cy - DNA_size:cy + DNA_size+1]
        event_popu[i] = {'N1':len(np.where(event_type_cut==1)[0]),\
                         'N2':len(np.where(event_type_cut==2)[0]),\
                         'N3':len(np.where(event_type_cut==3)[0]),\
                         'N4':len(np.where(event_type_cut==4)[0]),\
                         'event':int(event),
                         'ex':event_x,
                         'ey':event_y
                         }
        print('#############################')
        print(i)
        print(event_popu[i])
        print('#############################')
        t7  = time.time() 
        t67 = t7 - t6
        t67_sum += t67        

        if (i % frame_per_event == 0) & (i != 0):            
            print ("t01_sum "+str(t01_sum) +' sec')
            print ("t12_sum "+str(t12_sum) +' sec')
            print ("t2_01_sum "+str(t2_01_sum) +' sec')
            print ("t2_12_sum "+str(t2_12_sum) +' sec')
            print ("t2_23_sum "+str(t2_23_sum) +' sec')
            print ("t23_sum "+str(t23_sum) +' sec')
            print ("t34_sum "+str(t34_sum) +' sec')
            print ("t45_sum "+str(t45_sum) +' sec')
            print ("t56_sum "+str(t56_sum) +' sec')
            print ("t67_sum "+str(t67_sum) +' sec')

        i += 1
        
        ########################################################################################################################
        ########################################################################################################################

    t_list = [sum(dt_list[:ti]) for ti in range(len(dt_list))]
    dr_list = (np.array(dx_list)**2+np.array(dy_list)**2)**0.5
    cx_list, cy_list = np.array(cx_list),np.array(cy_list)
    
    SD_list = (np.array(cx_list) - int(RNA_size/2))**2 + \
        (np.array(cy_list) - int(RNA_size/2))**2
 
    X0,Y0 = RNA_size/2,RNA_size/2
    X0_nm,Y0_nm = X0*pixel_to_nm,Y0*pixel_to_nm
    
    X_nm, Y_nm = (cx_list - X0)*pixel_to_nm, (cy_list - Y0)*pixel_to_nm
    
    ###########################################################################
    ###########################################################################
    FIG = plt.figure(figsize=(14, 12))
    ax = FIG.add_subplot(111)
    xmin, xmax = np.nanmin(X_nm),np.nanmax(X_nm)
    ymin, ymax = np.nanmin(Y_nm),np.nanmax(Y_nm)
    tmin, tmax = np.nanmin(t_list),np.nanmax(t_list)
    T = tmax - tmin    
    DX, DY = xmax - xmin , ymax - ymin 
    DXY = max(DX,DY)/2 * 1.05
    X0, Y0 = (xmin + xmax)/2, (ymin + ymax)/2    
    Xmin, Xmax = X0-DXY, X0+DXY
    Ymin, Ymax = Y0-DXY, Y0+DXY
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    def PLOT_TRAJECTORY(X,Y,t,AX):
        AX.plot(X,Y,lw = 0)
        l = COLOR_LINE_PLOT(X,Y,t,AX)
        if DX < DY: # For vertical trajectories, the scale bar is on the right side
            AX.plot([xmax+20,xmax+120],[ymax-50,ymax-50],c='k', lw= 10, solid_capstyle='butt')
        if DX >= DY:# For horizontal trajectories, scale bar is lower side
            AX.plot([xmax-100,xmax],[ymax+50,ymax+50],c='k', lw= 10, solid_capstyle='butt')   
        FIG.colorbar(l, ax=AX,shrink=0.6,label= 'Time (s)',orientation = 'vertical')
        AX.set_aspect('equal')
        AX.set(xlabel='X (nm)',ylabel='Y (nm)')
        AX.invert_yaxis()
      #  ax1.set(xlim=(min(df['X (nm)']),max1(df['X (nm)']+120)), ylim=(max1(df['Y (nm)'])*1.05,min(df['Y (nm)'])*1.05) )
        AX.axis('off')
        
    ax.set_title( trajfol.replace(os.sep,'\n') +'\n'+ ('00000'+str(i_traj))[-4:]+".png" \
                 + '\n' + 'Net displacement:' + str( round((X_nm[-1]**2 + Y_nm[-1]**2)**0.5 )) + ' nm'\
                 + '\n' + 'Time: ' + str(round(tmax,1)) + ' sec' ,loc ='left')
        
    PLOT_TRAJECTORY(X_nm,Y_nm,t_list,ax)
    plt.savefig(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+"_trace.png",dpi = 50)    
  #  plt.show()
    plt.close()
    ###########################################################################
    ###########################################################################
    
    plt.figure(figsize=(12, 18))
    gs = GridSpec(5, 4)
    
    ax = plt.subplot(gs[0:2,0:2])

    
    """
    ax.plot(cx_list, cy_list, lw = 0.5, color = 'k',zorder = 0)
    s = ax.scatter((cx_list - X0)*pixel_to_nm, (cy_list - Y0)*pixel_to_nm, c=t_list, cmap='jet',s= 5, vmin= 0, vmax=tmax)
    fig.colorbar(s ,ax = ax, label='Time (s)',shrink = 0.5)
    """
    l = COLOR_LINE_PLOT(X_nm, Y_nm ,t_list,ax )
    fig.colorbar(l, ax = ax, label = 'Time (s)',shrink = 0.5)           
    ax.set(xlabel='X (nm)',ylabel='Y (nm)',xlim=(-X0_nm,X0_nm),ylim=(-Y0_nm,Y0_nm))
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.grid()

    ax = plt.subplot(gs[0:2,2:4])
    l = COLOR_LINE_PLOT(X_nm, Y_nm ,t_list,ax )
    fig.colorbar(l, ax = ax, label = 'Time (s)',shrink = 0.5)           
    ax.set(xlabel='X (nm)',ylabel='Y (nm)',xlim=(-500,500),ylim=(500,-500))
    ax.set_aspect('equal')
    ax.grid()



    
    

    plt.subplot(gs[2,0])
    plt.title(('00000'+str(i_traj))[-4:])
    plt.step(t_list, np.array(SD_list)**0.5, lw = 0.5, color = 'k',zorder = 0)
    plt.scatter(t_list, np.array(SD_list)**0.5,  c=t_list, cmap='jet',s= 5, vmin= 0, vmax=tmax ,zorder = 2)
    plt.xlim(1, 1000), plt.ylim(1, 1000)
    plt.xlabel('Time (s)'), plt.ylabel('Net displacement / pixel')
    plt.colorbar(label='Time (s)')
    plt.xscale('log'), plt.yscale('log')
    plt.grid()
    
    #################################################################################################################
    # save data #####################################################################################################
    df_event = pd.DataFrame.from_dict(event_popu, orient='index')
    df_XY = pd.DataFrame({'t':t_list,'X (pixel)':cx_list,'Y (pixel)':cy_list,\
                  'X (nm)':np.array(cx_list)*pixel_to_nm,'Y (nm)':np.array(cy_list)*pixel_to_nm,\
                  'dr (nm)':dr_list*pixel_to_nm,})

    df_plot = pd.concat([df_XY,df_event],axis = 1)
    df_plot['N_bound'] =  df_plot['N2'] + df_plot['N3']
    
    df_plot.to_csv(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".csv")

    
    pd.DataFrame(RNA_state).to_pickle(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+"_RNAstate.gz")

    #################################################################################################################
    #################################################################################################################

    # Show time course    
    ax_state = plt.subplot(gs[3,:])
    ax_state.step(df_plot['t'],df_plot['N1'],where='post',color = 'navy',label ='k$_{on}^{DNA}$')
    ax_state.step(df_plot['t'],df_plot['N2'],where='post',color = 'g',label ='k$_{on}^{E}$')
    ax_state.step(df_plot['t'],df_plot['N3'],where='post',color = 'm',label ='k$_{cat}^{E}$')
    ax_state.step(df_plot['t'],df_plot['N4'],where='post',color = 'darkorange',label ='Hydrolyzed')
    ax_state.set(ylabel='Number of sites',xlabel='Time (s)',yscale='log',ylim=(0.9,300) )            
    ax_state.set_yticks([2,3,4,5,6,7,8,9,20,30,40,50,60,70,80,90,200,300], minor=True)            
    ax_state.set_yticklabels(["2","3","4","","","","","","","","","","","","","","",""], minor=True)            
    ax_state.set_yticks([1,10,100])            
    ax_state.set_yticklabels(["1","10","100"])            
    ax_state.grid(which = 'both',axis='y')
    ax_state.legend(loc ='upper left',fontsize = 15)
    
    ax_dr = plt.subplot(gs[4,:])
    ax_dr.scatter(df_plot['t'],df_plot['dr (nm)'],)
    ax_dr.step(df_plot['t'],df_plot['dr (nm)'],where='post')
    ax_dr.set(ylabel='Displacement (nm)',xlabel='Time (s)',ylim =(0,50))
    ax_dr.grid(which = 'both',axis='y')
    
    ti_dr_nonzero = np.where(df_plot['dr (nm)']!=0)[0]
    ti_dr_zero = np.append( np.where(df_plot['dr (nm)']==0)[0] , np.array([len(df_plot)]) )
    ti_ini_step = np.array(sorted(list(set(ti_dr_nonzero) & set(ti_dr_zero+1))) )
    ti_end_step = np.array(sorted(list(set(ti_dr_nonzero) & set(ti_dr_zero-1))) ) 
    for ti_ini, ti_end in zip(ti_ini_step,ti_end_step):
        ax_dr.axvspan(list(df_plot['t'])[ti_ini], list(df_plot['t'])[ti_end], facecolor='lightpink')
        ax_state.axvspan(list(df_plot['t'])[ti_ini], list(df_plot['t'])[ti_end], facecolor='lightpink')
       
    plt.tight_layout()
    plt.savefig(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".png", dpi = 75)    
   # plt.show()
    plt.close()    
    """
    dt = 1 # sec
    for i in range(len(df_traj)):
        t_for_MSD = np.linspace(0, max(t_list), int(max(t_list)/dt))
    """    

    
    return df_plot





#################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Main program
global khyb, konE, kcat
for khyb in khyb_list:
    for konE in konE_list:
        for kcat in kcat_list:
            foli+=1
            foli_str =  ('000' + str(foli))[-3:]
            khyb_str = f"khyb={khyb:.2f}"
            kcat_str = f"kcatE={kcat:.1f}"
            konE_scaled = konE / 10**6
            konE_str = f"konE={konE_scaled:.1f}x106"
            # Combine string
            dirname = foli_str + '_' + f"{khyb_str}_{kcat_str}_{konE_str}"
            print('############################################')
            print('############################################')
            print(dirname)
            print('############################################')
            print('############################################')
            # make working folder
            workfol = globalfol + os.sep +  dirname
            folder_maker(workfol)            

            #start simulation for preset RNase H
            for RNaseH in RNaseH_list:
                print('RNaseH ' + str(RNaseH) + ' nM')
            
                keys = {'RNA_size':RNA_size,'N_simu':N_simu,'tmax':tmax_simu,'khyb':khyb,'kcat':kcat,'konE':konE,'RNaseH':RNaseH\
                        ,'particle_size':particle_size,'RNA_density':RNA_density,'DNA_density':DNA_density,'L_DNA':L_DNA\
                        ,'R_mobile_nm':R_mobile_nm,'frame_per_event':frame_per_event}                
                title_params = ['RNaseH','N_simu','tmax','frame_per_event']
                title_keys = {k: v for k, v in keys.items() if k in title_params}
                                
                savefol = workfol + os.sep + date + '_'+ str(title_keys).replace(':','=')
                folder_maker(savefol)                
                trajfol = savefol + os.sep + 'progress'
                folder_maker(trajfol)
                # save the parameter used.
                ####################################################################################
                df_params = pd.DataFrame(keys, index=[0])
                df_params.to_csv(savefol + os.sep + 'parameters.csv')
                ###############################################################################
                
                if SIMULATION:
                    # Run simulations for "N_simu" times
                    for i_traj in range(N_simu):
                        start = time.time()
                        trace_i_str = ('0000'+str(i_traj))[-3:]
                        trace_fol = trajfol + os.sep + trace_i_str
                        folder_maker(trace_fol)
                        print("\n################################################################################")
                        print("################################################################################")
                        print('simulation file is \n')
                        print(trajfol + os.sep + ('0000'+str(i_traj))[-4:] +'.png')
                        print('\n')    

                
                        if os.path.isfile(trajfol + os.sep + ('0000'+str(i_traj))[-4:] +'.png'):
                            print('file exist. Skip this')
                            continue        
                        #################################################################################
                        #################################################################################
                        df_traj = SIMULATE_TRAJ(i_traj,trace_fol, khyb,konE* RNaseH *10**-9, kcat )
                        df_traj.to_csv(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".csv")
                        #################################################################################
                        #################################################################################
                
                        trace_i_str = ('0000'+str(i_traj))[-3:]
                        if frame_per_event==1:
                          #  command = r'ffmpeg -y -framerate 100 -i "'+trace_i_str+os.sep+trace_i_str+ '_%06d.png" -c:v libx264 -preset slow -crf 22 '+trace_i_str+ '.avi'
                            command = r'ffmpeg -y -framerate 100 -i "'+trace_i_str+os.sep+trace_i_str+ '_%06d.png" -vf crop=iw:680:0:0 -pix_fmt yuv420p -c:v libx264 -preset slow  -crf 18 '+trace_i_str+ '.mp4'            
                            print(command)
                            p = subprocess.run(command, stdout=subprocess.PIPE, shell=True, cwd=trajfol,universal_newlines=True)




                if SUMMERIZE:
                    df_traj_sum = pd.DataFrame()
                    for i_traj in tqdm(range(N_simu)):
                        print(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".csv")
                        df = pd.read_csv(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".csv", index_col = None)
                        df_traj_sum = pd.concat([df_traj_sum,df])
                    df_traj_sum.to_pickle(savefol + os.sep + 'df_traj.gs')
                   
                    
                    # plot trajectories
                    def PLOT_ALL_TRAJ(df_traj_sum, MIN, MAX,BIN):
                        TMAX = max(df_traj_sum['t'])
                        fig = plt.figure(figsize=(7*2,7),tight_layout=True)
                        
                        ax = fig.add_subplot(121) 
                        ax.set_title(trajfol.replace(os.sep , '\n') +'\n\n',loc = 'left',fontsize= 10)
                        net_list = []
                        life_list = []
                        for i_traj in tqdm(range(N_simu)):
                            df = pd.read_csv(trajfol + os.sep + ('00000'+str(i_traj))[-4:]+".csv", index_col = None)
                            X,Y = np.array(df['X (nm)'] - df['X (nm)'].iloc[0]), np.array(df['Y (nm)'] - df['Y (nm)'].iloc[0])
                            ax.plot(X,Y,lw = 0)
                    
                            points = np.array([X, Y]).T.reshape(-1, 1, 2)
                            segments = np.concatenate([points[:-1], points[1:]], axis=1)
                            lc = LineCollection(segments, cmap='jet')
                            lc.set_array(df['t'])
                            lc.set_linewidth(0.8)
                            lc.set_clim(vmin = 0, vmax = TMAX)
                            line = ax.add_collection(lc)
                            
                            net = (X[-1]**2 + Y[-1]**2) ** 0.5
                            net_list.append(net)
                            life_list.append(max(df['t']))
                            
                        ax.set_aspect('equal')
                        ax.set(xlabel='X (nm)',ylabel='Y (nm)')
                        ax.set(xlim=(MIN, MAX),ylim=(MIN, MAX),xticks = np.arange(MIN, MAX+1,BIN),yticks = np.arange(MIN, MAX+1,BIN))
                        fig.colorbar(line, ax = ax,label ='Time (s)',shrink=0.7)
                        ax.grid()
                        
                        
                        ax = fig.add_subplot(222)  
                        ax.axvline(x=np.mean(net_list),color = 'r', ls = 'dashed',lw = 2)
                        ax.hist(net_list, bins = np.logspace(-1,4,60))
                        ax.set(xscale='log')
                        ax.set(xlabel='Net displacement (nm)',ylabel='Counts',xlim =(10**-1,10**4))
                        ax.legend(title= 'Mean = '+ str(round(np.mean(net_list))) + ' nm'\
                                         '\nMin = '+ str(round(np.min(net_list))) + ' nm'\
                                         '\nMax = '+ str(round(np.max(net_list))) + ' nm')
                        ax.grid()
                        
                        ax = fig.add_subplot(224)  
                        t_surv = np.array( sorted(life_list) )
                        N_surv = np.array( [ len(np.where( t_surv >= ti)[0])-1 for ti in t_surv]) / len(life_list)
                        
                        t_surv = np.append(np.array([0]),t_surv)
                        N_surv = np.append(np.array([1]),N_surv)
                        
                        ax.step(t_surv, N_surv,where='post')
                        ax.scatter(t_surv, N_surv)
                        ax.set_ylim(0,1.1)
                        ax.set_xlim(0.1,10000)
                        ax.set_xscale('log')
                        ax.set_xlabel('Time (s)')
                        ax.set_ylabel('Fraction bound')
                        ax.grid()
                        
                        plt.savefig(savefol + os.sep + 'All_trajectories_'+str(MIN) + '_to_' +str(MAX) + '_nm.png')
                     #   plt.show()
                        plt.clf()
                        plt.close()
                    
                    if not os.path.isfile(savefol + os.sep + 'All_trajectories_'+str(-3000) + '_to_' +str(3000) + '_nm.png'):                        
                        PLOT_ALL_TRAJ(df_traj_sum, -3000, 3000,1000)
                        PLOT_ALL_TRAJ(df_traj_sum, -2000, 2000,1000)
                        PLOT_ALL_TRAJ(df_traj_sum, -1000, 1000,500)
                        PLOT_ALL_TRAJ(df_traj_sum, -800, 800,400)
                        PLOT_ALL_TRAJ(df_traj_sum, -400, 400,200)
                        PLOT_ALL_TRAJ(df_traj_sum, -200, 200,100)
                    
                #################################################################################
                


