# -*- coding: utf-8 -*-
"""
######################################################
#             1D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains functions to do boundary conditions:
    -
    
Features:
    -

Desired:
    -
    -
    
"""

import numpy as np
#import string as st

class BCs():
    def __init__(self, BC_dict, dx):
        self.BCs=BC_dict
        self.dx=dx
        
    # Energy BCs
    def Energy(self, E, T_prev, dt, rho, Cv, hx):
        # Left face
        if self.BCs['bc_left_E'][0]=='T':
            E[0]=self.BCs['bc_left_E'][1]*rho[0]*Cv[0]
                
        else:
            if self.BCs['bc_left_E'][0]=='F':
                q=self.BCs['bc_left_E'][1]
                Bi=0
                
            else:
                q=self.BCs['bc_left_E'][1][0]*self.BCs['bc_left_E'][1][1] # h*Tinf
                Bi=-self.BCs['bc_left_E'][1][0]*T_prev[0] # h*Tij
            
            E[0]+=(Bi+q)*dt/hx[0]
            
        # Right face
        if self.BCs['bc_right_E'][0]=='T':
            E[-1]=self.BCs['bc_right_E'][1]*rho[-1]*Cv[-1]
                
        else:
            if self.BCs['bc_right_E'][0]=='F':
                q=self.BCs['bc_right_E'][1]
                Bi=0
                
            else:
                q=self.BCs['bc_right_E'][1][0]*self.BCs['bc_right_E'][1][1] # h*Tinf
                Bi=-self.BCs['bc_right_E'][1][0]*T_prev[-1] # h*Tij
            
            E[-1]+=(Bi+q)*dt/hx[-1]
            
        # Apply radiation BCs
        if self.BCs['bc_left_rad']!='None':
            E[0]+=dt/hx[0]*\
                self.BCs['bc_left_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_left_rad'][1]**4-T_prev[0]**4)
        if self.BCs['bc_right_rad']!='None':
            E[-1]+=dt/hx[-1]*\
                self.BCs['bc_right_rad'][0]*5.67*10**(-8)*\
                (self.BCs['bc_right_rad'][1]**4-T_prev[-1]**4)
            
    # Conservation of mass BCs
    def mass(self, m, P, Ax, Ay):
        # Left face
        for i in range(len(self.BCs['bc_left_mass'])/3):
            st=self.BCs['bc_left_mass'][2+3*i][0]
            en=self.BCs['bc_left_mass'][2+3*i][1]
            # Gradient
            if self.BCs['bc_left_mass'][3*i]=='grad':
                m[st:en,0]=m[st:en,1]-self.BCs['bc_left_mass'][1+3*i]*self.dx[st:en,0]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=m[-1,1]-self.BCs['bc_left_mass'][-2]*self.dx[-1,0]
            # Pressure flux
            elif self.BCs['bc_left_mass'][3*i]=='grad_P':
                m[st:en,0]=m[st:en,1]-self.BCs['bc_left_mass'][1+3*i]*self.dx[st:en,0]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=m[-1,1]-self.BCs['bc_left_mass'][-2]*self.dx[-1,0]
            # Constant
            else:
                m[st:en,0]=self.BCs['bc_left_mass'][1+3*i]
                if len(self.BCs['bc_left_mass'])/3-i==1:
                    m[-1,0]=self.BCs['bc_left_mass'][-2]
        # Right face
        for i in range(len(self.BCs['bc_right_mass'])/3):
            st=self.BCs['bc_right_mass'][2+3*i][0]
            en=self.BCs['bc_right_mass'][2+3*i][1]
            if self.BCs['bc_right_mass'][3*i]=='grad':
                m[st:en,-1]=self.BCs['bc_right_mass'][1+3*i]*self.dx[st:en,-1]+m[st:en,-2]
                if len(self.BCs['bc_right_mass'])/3-i==1:
                    m[-1,-1]=self.BCs['bc_right_mass'][-2]*self.dx[-1,-1]+m[-1,-2]
        
        # South face
        for i in range(len(self.BCs['bc_south_mass'])/3):
            st=self.BCs['bc_south_mass'][2+3*i][0]
            en=self.BCs['bc_south_mass'][2+3*i][1]
            if self.BCs['bc_south_mass'][3*i]=='grad':
                m[0,st:en]=m[1,st:en]-self.BCs['bc_south_mass'][1+3*i]*self.dy[0,st:en]
                if len(self.BCs['bc_south_mass'])/3-i==1:
                    m[0,-1]=m[1,-1]-self.BCs['bc_south_mass'][-2]*self.dy[0,-1]
                    
        # North face
        for i in range(len(self.BCs['bc_north_mass'])/3):
            st=self.BCs['bc_north_mass'][2+3*i][0]
            en=self.BCs['bc_north_mass'][2+3*i][1]
            if self.BCs['bc_north_mass'][3*i]=='grad':
                m[-1,st:en]=self.BCs['bc_north_mass'][1+3*i]*self.dy[-1,st:en]+m[-2,st:en]
                if len(self.BCs['bc_north_mass'])/3-i==1:
                    m[-1,-1]=self.BCs['bc_north_mass'][-2]*self.dy[-1,-1]+m[-2,-1]
        return 0
    
    # Pressure BCs (eventually lead to momentum)
    def P(self, P):
        # Left face
        if self.BCs['bc_left_P'][0]=='grad':
            P[0]=P[1]-self.BCs['bc_left_P'][1]*self.dx[0]
        
        # Right face
        if self.BCs['bc_right_P'][0]=='grad':
            P[-1]=self.BCs['bc_right_P'][1]*self.dx[-1]+P[-2]
            
        return 0