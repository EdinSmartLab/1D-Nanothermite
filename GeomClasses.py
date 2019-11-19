# -*- coding: utf-8 -*-
"""
######################################################
#             1D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the 1D domain classes:
    -holds conservative variable (energy) at each node
    -holds thermal properties at each node
    -holds x coordinate arrays
    -holds dx discretization array
    -calculates thermal properties
    -meshing function (biasing feature not functional in solver)
    -function to return temperature given conservative variable (energy)
    -calculate CV 'volume' at each node

Requires:
    -length of domain
    -number of nodes across length
    -values for thermal properties (depends on what method of calculation)
    

"""
import numpy as np
import string as st
from MatClasses import Cp, therm_cond

class OneDimLine():
    def __init__(self, settings, Species, solver, rank):
        
        self.L=settings['Length']
        self.Nx=settings['Nodes_x']
        self.model=settings['Model']
        self.x=np.zeros(self.Nx)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.rank=rank
        self.porosity_0=settings['Porosity']
        
        # Variables for conservation equations
        self.E=np.zeros(self.Nx) # Lumped energy
        self.max_iter=settings['Max_iterations']
        self.conv=settings['Convergence']
        
        # Thermal properties (solid model)
        self.k=settings['k_s']
        self.k_mode=settings['k_model']
        self.rho=settings['rho_IC']
        self.Cv=settings['Cv_s']
        
        # Process the thermal properties options
        if (type(self.k) is str):
            self.k=st.split(self.k, ',')
        if (type(self.Cv) is str):
            self.Cv=st.split(self.Cv, ',')
                
        # Species model options
        if self.model=='Species':
            self.species_keys=['g','s']
            self.rho=st.split(self.rho, ',')
            self.mu=settings['Darcy_mu']
            self.part_diam=settings['Carmen_diam']
            self.R=settings['gas_constant']
            self.Cv_g=Species['Cv_g']
            self.Cp_g=Species['Cp_g']
            self.k_g=Species['k_g']
            # Process thermal properties options
            if (type(self.Cv_g) is str):
                self.Cv_g=st.split(self.Cv_g, ',')
            if (type(self.Cp_g) is str):
                self.Cp_g=st.split(self.Cp_g, ',')
            if (type(self.k_g) is str):
                self.k_g=st.split(self.k_g, ',')
        
        
        # Object declarations for property calculations
        self.Cp_calc=Cp()
        self.k_calc=therm_cond()
        
        # Biasing options       
        self.xbias=[settings['bias_type_x'], settings['bias_size_x']]
        self.isMeshed=False
        
        # MPI information (will be set by another function)
        self.proc_left=-1
        self.proc_right=-1
        
    # Discretize domain and save dx and dy
    def mesh(self):
        # Discretize x
        if self.xbias[0]=='OneWayUp':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
#            print 'One way biasing in x: smallest element at x=%2f'%self.L
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='OneWayDown':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
#            print 'One way biasing in x: smallest element at x=0'
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayEnd':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
#            print 'Two way biasing in x: smallest elements at x=0 and %2f'%self.L
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayMid':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
#            print 'Two way biasing in x: smallest elements around x=%2f'%(self.L/2)
#            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        else:
            self.dx[:]=self.L/(self.Nx-1)
#            print 'No biasing schemes specified in x'
        
        for i in range(self.Nx-1):
            self.x[i+1]=self.x[i]+self.dx[i]
        
        self.X=self.x
        
        self.isMeshed=True
    
    # Define other variables for calculations after MPI
    def create_var(self, Species):
        self.eta=np.zeros_like(self.E) # extent of reaction
        self.P=np.zeros_like(self.E) # pressure
        self.T_guess=np.zeros_like(self.E)
        self.porosity=np.ones_like(self.E)*self.porosity_0
        
        # Species
        self.rho_species={}
        try:
            self.rho_0=np.ones_like(self.E)*self.rho*(1-self.porosity_0)
        except:
            self.rho_0=np.zeros_like(self.E)
        por=[self.porosity,(1-self.porosity)]
        if self.model=='Species':
            for i in range(len(self.species_keys)):
                self.rho_species[self.species_keys[i]]=np.ones_like(self.E)\
                    *float(self.rho[i])*por[i]
            self.rho_0=self.rho_species[self.species_keys[1]]
            self.perm=self.porosity**3*self.part_diam**2\
                /(180*(1-self.porosity)**2)
        
    # Calculate and return dimensions of CV
    def CV_dim(self):
        hx=np.zeros_like(self.E)
        
        hx[0]      =0.5*(self.dx[0])
        hx[1:-1]   =0.5*(self.dx[1:-1]+self.dx[:-2])
        hx[-1]     =0.5*(self.dx[-1])
        
        return hx
    
    # Calculate temperature and thermodynamic properties
    # Have it return only rho_s; where is it used?
    def calcProp(self, T_guess=300, init=False):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        rhoC=np.zeros_like(self.eta)
        Cv=np.zeros_like(self.eta)
        Cp=np.zeros_like(self.eta)
        
        ############ Specific heat of solid phase (either model)
        # Solid phase (eta dependent)
        if (type(self.Cv) is list) and (self.Cv[0]=='eta'):
            Cv=self.eta*float(self.Cv[2])+(1-self.eta)*float(self.Cv[1])
        
        # Solid phase (temperature dependent for given element)
        elif (type(self.Cv) is list) and (self.Cv[1]=='Temp'):
            # Constant temperature value
            if len(self.Cv)>2:
                Cv=self.Cp_calc.get_Cv(np.ones_like(self.E)*float(self.Cv[2]), self.Cv[0])
            # Temperature dependent
            else:
                Cv=self.Cp_calc.get_Cv(T_guess, self.Cv[0])
        
        # Solid phase (constant)
        else:
            Cv[:]=self.Cv
        
        ############ Thermal conductivity of solid phase (either model)
        # Solid phase (eta dependent)
        if (type(self.k) is list) and (self.k[0]=='eta'):
            k=self.eta*float(self.k[2])+(1-self.eta)*float(self.k[1])
        
        # Solid phase (temperature dependent for given element)
        elif (type(self.k) is list) and (self.k[1]=='Temp'):
            # Constant temperature value
            if len(self.k)>2:
                k=self.k_calc.get_k(np.ones_like(self.E)*float(self.k[2]), self.k[0])
            # Temperature dependent
            else:
                k=self.k_calc.get_k(T_guess, self.k[0])
        
        # Solid phase (constant)
        else:
            k[:]=self.k
        
        ############ When species model is active
        if self.model=='Species':
            k_g=np.zeros_like(self.eta)
            # Changing porosity/permeability
            self.porosity=self.porosity_0+\
                self.rho_species[self.species_keys[0]]/self.rho_0*(1-self.porosity_0)
            self.perm=self.porosity**3*self.part_diam**2\
                /(72*(1-self.porosity)**2)
            
            # Heat capacity of Solid phase
            rhoC=self.rho_species[self.species_keys[1]]*Cv
#            rhoC=self.rho*(1-self.porosity)*Cv # REPLICATE CASE 10 (CASE 10d,e)
            
            # Heat capacity of Gas phase
            if (type(self.Cv_g) is list) and (self.Cv_g[0]=='eta'):
                Cv=self.eta*float(self.Cv_g[2])+(1-self.eta)*float(self.Cv_g[1])
            
            # Gas phase (temperature dependent for given element)
            elif (type(self.Cv_g) is list) and (self.Cv_g[1]=='Temp'):
                # Constant temperature value
                if len(self.Cv_g)>2:
                    Cv=self.Cp_calc.get_Cv(np.ones_like(self.E)*float(self.Cv_g[2]), self.Cv_g[0])
                # Temperature dependent
                else:
                    Cv=self.Cp_calc.get_Cv(T_guess, self.Cv_g[0])
            
            # Gas phase (constant)
            else:
                Cv[:]=self.Cv_g
            rhoC+=self.rho_species[self.species_keys[0]]*Cv
            
            # Temperature calculation
            T=self.E/rhoC
            self.T_guess=T
            # Iteratively solve temperature (temperature dependent properties)
#            T_0=np.ones_like(self.eta)
#            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
#            i=0
#            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
#                T_0=T.copy()
#                rhoC=(1-self.porosity)*self.rho_species[self.species_keys[1]]*Cv
#                rhoC+=self.porosity*self.rho_species[self.species_keys[0]]*self.Cp_calc.get_Cv(T_guess, self.pore_gas)
#                T=self.E/rhoC
#                i+=1
#                if init:
#                    break
            # Specific heat of Gas phase (Cp)
#                Cv_Al2O3=self.Cp_calc.get_Cp(np.ones_like(T)*2327,'Al2O3')
#                Cv_Cu=self.Cp_calc.get_Cp(np.ones_like(T)*2843,'Cu')
#                Cp=(0.351*Cv_Al2O3+0.649*Cv_Cu)
            
            if (type(self.Cp_g) is list) and (self.Cp_g[0]=='eta'):
                Cp=self.eta*float(self.Cp_g[2])+(1-self.eta)*float(self.Cp_g[1])
            
            # (temperature dependent for given element)
            elif (type(self.Cp_g) is list) and (self.Cp_g[1]=='Temp'):
                # Constant temperature value
                if len(self.Cp_g)>2:
                    Cp=self.Cp_calc.get_Cp(np.ones_like(self.E)*float(self.Cp_g[2]), self.Cp_g[0])
                # Temperature dependent
                else:
                    Cp=self.Cp_calc.get_Cp(T_guess, self.Cp_g[0])
            
            # Gas phase (constant)
            else:
                Cp[:]=self.Cp_g
            
            # Thermal conductivity of gas phase
            # eta dependent
            if (type(self.k_g) is list) and (self.k_g[0]=='eta'):
                k_g=self.eta*float(self.k_g[2])+(1-self.eta)*float(self.k_g[1])
            
            # temperature dependent for given element
            elif (type(self.k_g) is list) and (self.k_g[1]=='Temp'):
                # Constant temperature value
                if len(self.k_g)>2:
                    k_g=self.k_calc.get_k(np.ones_like(self.E)*float(self.k_g[2]), self.k_g[0])
                # Temperature dependent
                else:
                    k_g=self.k_calc.get_k(T_guess, self.k_g[0])
            
            # constant
            else:
                k_g[:]=self.k_g
            
            # Thermal conductivity models
            if self.k_mode=='Parallel':
                k=self.porosity*k_g+(1-self.porosity)*k
            elif self.k_mode=='Geometric':
                k=k*(k_g/k)**(self.porosity)
            elif self.k_mode=='Series':
                k=(self.porosity/k_g+(1-self.porosity)/k)**(-1)
            
        ############ Plain heat transfer model
        else:
            rho[:]=self.rho*(1-self.porosity_0)
            
            rhoC=rho*Cv
            T=self.E/rhoC
            self.T_guess=T
            # Iteratively solve temperature (temperature dependent properties)
#            T_0=np.ones_like(self.eta)
#            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
#            i=0
#            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
#                T_0=T.copy()
#                rhoC=rho*Cv
#                T=self.E/rhoC
#                i+=1
#                if init:
#                    break 
        if init:
            return rhoC
        else:
            return T, k, rhoC, Cp
