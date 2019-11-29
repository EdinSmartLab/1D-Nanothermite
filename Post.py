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
import os
import sys
import string as st
from matplotlib import pyplot
from FileClasses import FileIn

pyplot.ioff()

def interpolate(k1, k2, func):
    if func=='Linear':
        return 0.5*k1+0.5*k2
    else:
        return 2*k1*k2/(k1+k2)

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
try:
    input_file=FileIn('Input_file.txt',False)
#    open('Input_file.txt')
except:
    try:
        input_file=FileIn('Input_file_stats.txt',False)
    except:
        sys.exit('Input file missing')

titles=['g','s']

settings={}
sources={}
Species={}
BCs={}
input_file.Read_Input(settings, sources, Species, BCs)
xmax=float(settings['Length'])*1000
try:
    settings['rho_IC']=st.split(settings['rho_IC'], ',')
except:
    settings['rho_IC']=float(settings['rho_IC'])

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
    if st.find(sources['Source_Kim'],'True')>=0:
        eta=np.load('eta_'+time+'.npy', False)
        Y_tot=0.0
    
    # 1D temperature profile at centreline
     
    fig=pyplot.figure(figsize=(6, 6))
    pyplot.plot(X*1000, T)
    pyplot.xlabel('$x$ (mm)')
    pyplot.ylabel('T (K)')
    pyplot.xlim([xmin,xmax])
    pyplot.title('Temperature distribution t='+time+' ms')
    fig.savefig('T_'+time+'.png',dpi=300)
    pyplot.close(fig)
    
    if st.find(sources['Source_Kim'],'True')>=0:
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
            phi=sources['A0']*(1-eta)*np.exp(-sources['Ea']/8.314/T)
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
        # Pressure plot
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
        Y_tot+=np.sum(Y_0)/len(Y_0)
    
    # Velocity plot
    por=settings['Porosity']+\
        (1-Y_0/(float(settings['rho_IC'][1])*settings['Porosity']))\
        *(1-settings['Porosity'])
    perm=por**3*settings['Carmen_diam']**2/(settings['Kozeny_const']*(1-por)**2)
    u=-interpolate(perm[1:], perm[:-1], settings['diff_interpolation'])\
        /settings['Darcy_mu']*(P[1:]-P[:-1])/(X[1:]-X[:-1])
    fig=pyplot.figure(figsize=(6, 6))
    pyplot.plot(X[1:]*1000,u)
    pyplot.xlabel('$x$ (mm)')
    pyplot.ylabel('Velocity (m/s)')
    pyplot.xlim([xmin,xmax])
    pyplot.title('Darcy velocity t='+time+' ms');
    fig.savefig('u_'+time+'.png',dpi=300)
    pyplot.close(fig)
        
    print 'Processed '+time
    print '     Mass balance residual: %8f'%(float(Y_tot))

print '\nPost-processing complete'