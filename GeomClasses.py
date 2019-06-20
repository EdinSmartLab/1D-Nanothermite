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
import copy
from MatClasses import Diff_Coef, Cp, therm_cond

class OneDimLine():
    def __init__(self, settings, Species, solver, rank):
        
        self.L=settings['Length']
        self.Nx=settings['Nodes_x']
        self.x=np.zeros(self.Nx)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.rank=rank
        self.porosity=settings['Porosity']
        self.pore_gas=settings['pore_gas']
        self.species_keys=[]
        if bool(Species):
            self.species_keys=Species['keys']
        
        # Variables for conservation equations
        self.E=np.zeros(self.Nx) # Lumped energy
        self.max_iter=settings['Max_iterations']
        self.conv=settings['Convergence']
        
        # Thermal properties
        self.k=settings['k']
        self.rho=settings['rho']
        self.Cv=settings['Cp']
#        if type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
#            line=st.split(self.rho, ',')
#            self.rho0=float(line[1])
#            self.rho1=float(line[2])
        if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
            line=st.split(self.Cv, ',')
            self.Cv0=float(line[1])
            self.Cv1=float(line[2])
        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
            line=st.split(self.k, ',')
            self.k0=float(line[1])
            self.k1=float(line[2])
        self.mu=settings['Darcy_mu']
        self.perm=self.porosity**3*settings['Particle_diam']**2\
            /(72*(1-self.porosity)**2)
        self.R=settings['gas_constant']
        
        self.Diff=Diff_Coef()
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
        
        # Species
        self.rho_species={}
        self.rho_0=np.zeros_like(self.E)
        por=[self.porosity,(1-self.porosity)]
        if bool(self.species_keys):
            for i in range(len(self.species_keys)):
                self.rho_species[self.species_keys[i]]=np.ones_like(self.E)\
                    *Species['Specie_IC'][i]
                self.rho_0+=por[i]*self.rho_species[self.species_keys[i]]
        
          
    # Calculate and return dimensions of CV
    def CV_dim(self):
        hx=np.zeros_like(self.E)
        
        hx[0]      =0.5*(self.dx[0])
        hx[1:-1]   =0.5*(self.dx[1:-1]+self.dx[:-2])
        hx[-1]     =0.5*(self.dx[-1])
        
        return hx
    
    # Calculate temperature dependent properties
#    def calcProp(self):
#        k=np.zeros_like(self.eta)
#        rho=np.zeros_like(self.eta)
#        Cv=np.zeros_like(self.eta)
#        Cp=np.zeros_like(self.eta)
#        D=copy.deepcopy(self.rho_species)
#        
#        # Species densities and specific heat
#        if bool(self.rho_species):
##            m_tot=np.zeros_like(self.E) # Use rho to be bulk rho
#            for i in range(len(self.species_keys)):
##                self.rho_species[self.species_keys[i]]=\
##                    self.m_species[self.species_keys[i]]/(por[i]*self.CV_vol())
#                Cv+=self.rho_species[self.species_keys[i]]*self.Cv_species[self.species_keys[i]]
#                Cp+=self.rho_species[self.species_keys[i]]*self.Cp_species[self.species_keys[i]]
#                rho+=self.rho_species[self.species_keys[i]]
#            Cv/=rho
#            Cp/=rho
#        
#        # Calculate properties based on eta or constant
#        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
#            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
#        elif type(self.k) is float:
#            k[:]=self.k
#        if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
#            Cv=self.eta*self.Cv1+(1-self.eta)*self.Cv0
#        elif type(self.Cv) is float:
#            Cv[:]=self.Cv
#        if type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
#            rho=self.eta*self.rho1+(1-self.eta)*self.rho0
#        elif type(self.rho) is float:
#            rho[:]=self.rho
#        
#        # Mass diffusion coefficient; g, s
#        if bool(D):
#            for i in self.species_keys:
#                D[i][:]=self.Diff.get_Diff(300,i)
#        
#        
#        return k, rho, Cv, Cp, D
    
    # Calculate temperature and thermodynamic properties
    # Have it return only rho_s; where is it used?
    def calcProp(self, T_guess=300, init=False):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        Cv=np.zeros_like(self.eta)
        Cp=np.zeros_like(self.eta)
        D=copy.deepcopy(self.rho_species)
        
#        # EXPERIMENTAL: Change porosity as function of eta
#        self.porosity
        por=[self.porosity,(1-self.porosity)]
        
        # Density
        if (type(self.rho) is str) and (st.find(self.rho, 'spec')>=0):
            for i in range(len(self.species_keys)):
                rho+=por[i]*self.rho_species[self.species_keys[i]]
#        elif type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
#            rho=self.eta*self.rho1+(1-self.eta)*self.rho0
        else:
            rho[:]=self.rho*(1-self.porosity)+self.Cp_calc.rho[self.pore_gas]*self.porosity
        
        # Specific heat (Cv)
        if (type(self.Cv) is str) and (st.find(self.Cv, 'spec')>=0):
            T_0=np.ones_like(self.eta)
            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
            i=0
            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
                T_0=T.copy()
                # Reactants
                Cv_Al=self.Cp_calc.get_Cv(T_0,'Al')
                Cv_CuO=self.Cp_calc.get_Cv(T_0,'CuO')
                # Products (gaseous phase, need to account for Cp vs Cv)
                Cv_Al2O3=self.Cp_calc.get_Cv(T_0,'Al2O3')
                Cv_Cu=self.Cp_calc.get_Cv(T_0,'Cu')
                
#                Cv=self.eta*(0.351*Cv_Al2O3+0.649*Cv_Cu)\
#                    +(1-self.eta)*(0.186*Cv_Al+0.814*Cv_CuO)
                
                Cv=(self.rho_species[self.species_keys[0]]*por[0]*(0.351*Cv_Al2O3+0.649*Cv_Cu)\
                    +self.rho_species[self.species_keys[1]]*por[1]*(0.186*Cv_Al+0.814*Cv_CuO))/rho
                
                T=self.E/Cv/rho
                i+=1
                if init:
                    break
        elif (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
#            Cv=self.eta*self.Cv1+(1-self.eta)*(self.Cv0)
#            Cv=(self.eta*self.Cv1+(1-self.eta)*self.Cv0)*(1-self.porosity)\
#                +self.Cp_calc.get_Cv(300, self.pore_gas)*self.porosity
            T_0=np.ones_like(self.eta)
            T=np.ones_like(self.eta)*T_guess # Initial guess for temperature
            i=0
            while np.amax(np.abs(T_0-T)/T)>self.conv and i<self.max_iter:
                T_0=T.copy()
                try:
                    Cv=(self.rho_species[self.species_keys[0]]*(self.eta*self.Cv1+(1-self.eta)*self.Cv0)*(1-self.porosity)\
                    +self.rho_species[self.species_keys[1]]*self.Cp_calc.get_Cv(T_0, self.pore_gas)*self.porosity)\
                        /rho
                except:
                    Cv=(self.rho*(self.eta*self.Cv1+(1-self.eta)*self.Cv0)*(1-self.porosity)\
                    +self.Cp_calc.rho[self.pore_gas]*self.Cp_calc.get_Cv(T_0, self.pore_gas)*self.porosity)\
                        /rho
                T=self.E/Cv/rho
                i+=1
                if init:
                    break
        else:
            Cv[:]=self.Cv*(1-self.porosity)\
                +self.Cp_calc.get_Cv(T_guess, self.pore_gas)*self.porosity
            T=self.E/Cv/rho
        
        # Specific heat (Cp) and diffusion coefficients (Dij)
        if bool(self.rho_species):
#            for i in range(len(self.species_keys)):
##                Cp+=self.rho_species[self.species_keys[i]]*por[i]*self.Cp_species[self.species_keys[i]]/rho
#                D[self.species_keys[i]][:]=self.Diff.get_Diff(T,self.species_keys[i])
            # Products (only these have gas phases)
            if self.species_keys[0]=='Ar':
                # Argon as only gas specie
                Cp=self.Cp_calc.get_Cp(T, 'Ar')
            else:
                # Special mix of products, no argon or air present
#                Cv_Al2O3=self.Cp_calc.get_Cp(T,'Al2O3')
                Cv_Al2O3=self.Cp_calc.get_Cp(np.ones_like(T)*2327,'Al2O3')
#                Cv_Cu=self.Cp_calc.get_Cp(T,'Cu')
                Cv_Cu=self.Cp_calc.get_Cp(np.ones_like(T)*2843,'Cu')
                
#                Cp=self.rho_species[self.species_keys[0]]*por[0]*(0.351*Cv_Al2O3+0.649*Cv_Cu)/rho
                Cp=(0.351*Cv_Al2O3+0.649*Cv_Cu)
#                Cp=Cv
                
        
        # Thermal conductivity
        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
#            ks=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
#            kf=self.k_calc.get_k(T, self.pore_gas)
#            k[:]=ks*(kf/ks)**(self.porosity)
        elif type(self.k) is float:
            k[:]=self.k
#            kf=self.k_calc.get_k(T, self.pore_gas)
#            k[:]=self.k*(kf/self.k)**(self.porosity)
        
        if init:
            return rho, Cv
        else:
            return T, k, rho, Cv, Cp, D
