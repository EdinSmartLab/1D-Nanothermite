# -*- coding: utf-8 -*-
"""
######################################################
#             1D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the post-processing script for 1D conduction:
    -Called from command line by:
        python Post-processing.py [Data directory relative to current directory]
    -Reads input file to get necessary parameters
    -Reads x meshgrid array (.npy) for graph output
    -Reads variable arrays (.npy files) and outputs graphs (.png) for 
    each time step in directory

Features:
    -Graphs of Temperature, reaction progress, reaction rate

Desired:
    -display variables available OR pass in desired graphs as argument for
    executing script
    
"""

import numpy as np
#import CoolProp.CoolProp as CP
import os
import sys
import string as st
from matplotlib import pyplot, cm
#from mpl_toolkits.mplot3d import Axes3D

pyplot.ioff()

print('######################################################')
print('#            1D Conduction Post-processing           #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')

inputargs=sys.argv
if len(inputargs)>1:
    inp_file=inputargs[1]
else:
    print 'Usage is: python Post.py [Input File]\n'
    print 'where\n'
    print '[Input File] is the Post-processing input file'
    print '***********************************'
    sys.exit('Post-processing halted')

##############################################################
#               Read post-processing file
##############################################################
try:
    fin=open(inp_file, 'r')
except:
    sys.exit('Cannot find post-processing input file')
    
for line in fin:
    if st.find(line, ':')>0 and st.find(line, '#')!=0:
        line=st.split(line, ':')
        if line[0]=='Directory':
            dir_files=st.split(line[1], '\n')[0]
        elif line[0]=='Times':
            if st.find(line[1], ',')>0:
                times=st.split(line[1], ',')
                times[-1]=st.split(times[-1], '\n')[0]
            else:
                times=st.split(line[1], '\n')[0]
        elif line[0]=='x_min':
            xmin=float(line[1])
        elif line[0]=='x_max':
            try:
                xmax=float(line[1])
            except:
                xmax=line[1]
        elif line[0]=='1D_Plots':
            OneD_graphs=line[1]
        elif line[0]=='Temp_min':
            temp_min=float(line[1])
        elif line[0]=='Temp_max':
            temp_max=float(line[1])
        elif line[0]=='Temp_pts':
            temp_pts=int(line[1])
        elif line[0]=='Phi_Plots':
            Phi_graphs=line[1]

fin.close()

try:
    os.chdir(dir_files)
except:
    sys.exit('Directory "'+dir_files+'" not found')

##############################################################
#               Read Solver file
##############################################################
A0=-1.0
Ea=-1.0
source='False'
try:
    input_file=open('Input_file.txt')
except:
    try:
        input_file=open('Input_file_stats.txt')
    except:
        sys.exit('Input file missing')

titles=[]
while A0<0 or Ea<0 or source=='False':
    line=input_file.readline()
    if st.find(line, 'Ea')==0:
        Ea=float(st.split(line, ':')[1])
    elif st.find(line, 'A0')==0:
        A0=float(st.split(line, ':')[1])
    elif st.find(line, 'Source_Kim')==0:
        source=st.split(line, ':')[1]
    elif st.find(line, 'Species')==0:
        titles=st.split(st.split(st.split(line, ':')[1], '\n')[0], ',')
    elif st.find(line, 'Length')==0 and type(xmax) is str:
        xmax=float(st.split(line, ':')[1])*1000
input_file.close()

##############################################################
#               Times to process (if ALL is selected)
##############################################################
if type(times) is str:
    times=os.listdir('.')
    i=len(times)
    j=0
    while i>j:
        if st.find(times[j],'T')==0 and st.find(times[j],'.npy')>0:
            times[j]=st.split(st.split(times[j],'_')[1],'.npy')[0]
            j+=1
        else:
            del times[j]
            i-=1

##############################################################
#               Figure details (NOT USED)
##############################################################
# NO OPTIONS YET
##############################################################
#               Generate graphs
##############################################################
X=np.load('X.npy', False)
for time in times:
    T=np.load('T_'+time+'.npy', False)
    if st.find(source,'True')>=0:
        eta=np.load('eta_'+time+'.npy', False)
        Y_tot=np.zeros_like(X)
    
    # 1D temperature profile at centreline
     
    fig=pyplot.figure(figsize=(6, 6))
    pyplot.plot(X*1000, T)
    pyplot.xlabel('$x$ (mm)')
    pyplot.ylabel('T (K)')
    pyplot.xlim([xmin,xmax])
    pyplot.title('Temperature distribution t='+time+' ms')
    fig.savefig('T_'+time+'.png',dpi=300)
    pyplot.close(fig)
    
    if st.find(source,'True')>=0:
        # Progress contour
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.plot(X*1000, eta)
        pyplot.xlabel('$x$ (mm)')
        pyplot.ylabel('$\eta$ (-)')
        pyplot.xlim([xmin,xmax])
        pyplot.title('Progress distribution t='+time+' ms');
        fig.savefig('eta_'+time+'.png',dpi=300)
        pyplot.close(fig)
        
        # Reaction rate contour
        if st.find(Phi_graphs,'True')>=0:
            phi=A0*(1-eta)*np.exp(-Ea/8.314/T)
            fig=pyplot.figure(figsize=(6, 6))
            pyplot.plot(X*1000, phi)
            pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            pyplot.xlabel('$x$ (mm)')
            pyplot.ylabel('$d\eta/dt$ ($s^{-1}$)')
            pyplot.xlim([xmin,xmax])
            pyplot.title('Reaction rate t='+time+' ms')
            fig.savefig('Phi_'+time+'.png',dpi=300)
            pyplot.close(fig)
    try:
        P=np.load('P_'+time+'.npy', False)
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.plot(X*1000,P)
        pyplot.xlabel('$x$ (mm)')
        pyplot.ylabel('Pressure (Pa)')
        pyplot.xlim([xmin,xmax])
        pyplot.title('Pressure t='+time+' ms');
        fig.savefig('P_'+time+'.png',dpi=300)
        pyplot.close(fig)
    except:
        print 'Processed '+time
        continue
    
        # Mass fraction contours
    for i in range(len(titles)):
        Y_0=np.load('rho_'+titles[i]+'_'+time+'.npy', False)
        fig=pyplot.figure(figsize=(6, 6))
        pyplot.plot(X*1000,Y_0)
        pyplot.xlabel('$x$ (mm)')
        pyplot.ylabel('$m$ ($kg/m^3$)')
        pyplot.xlim([xmin,xmax])
        pyplot.title('Density; $'+titles[i]+'$, t='+time+' ms');
        fig.savefig('rho_'+titles[i]+'_'+time+'.png',dpi=300)
        pyplot.close(fig)
        Y_tot+=Y_0
            
        
    print 'Processed '+time
    print '     Mass balance residual: %8f'%(np.amin(Y_tot))

print '\nPost-processing complete'