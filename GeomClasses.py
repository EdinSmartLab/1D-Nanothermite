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
from MatClasses import Diff_Coef

class OneDimLine():
    def __init__(self, settings, Species, solver, rank):
        
        self.L=settings['Length']
        self.Nx=settings['Nodes_x']
        self.x=np.zeros(self.Nx)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        self.rank=rank
        self.porosity=settings['Porosity']
        self.species_keys=[]
        if bool(Species):
            self.species_keys=Species['keys']
        
        # Variables for conservation equations
        self.E=np.zeros(self.Nx) # Lumped energy
        
        # Thermal properties
        self.k=settings['k']
        self.rho=settings['rho']
        self.Cv=settings['Cp']
        if type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
            line=st.split(self.rho, ',')
            self.rho0=float(line[1])
            self.rho1=float(line[2])
        if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
            line=st.split(self.Cv, ',')
            self.Cv0=float(line[1])
            self.Cv1=float(line[2])
        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
            line=st.split(self.k, ',')
            self.k0=float(line[1])
            self.k1=float(line[2])
        self.mu=settings['Darcy_mu']
        self.perm=settings['Darcy_perm']
        self.R=settings['gas_constant']
        
        self.Diff=Diff_Coef()
        
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
        
        # Species
#        self.m_species={}
#        self.mu_species={}
#        self.mv_species={}
        self.rho_species={}
        self.Cp_species={}
        self.rho_0=np.zeros_like(self.E)
        if bool(self.species_keys):
            i=0
            for key in self.species_keys:
#                self.m_species[key]=np.zeros_like(self.E)
#                self.mu_species[key]=np.zeros_like(self.E)
#                self.mv_species[key]=np.zeros_like(self.E)
                self.rho_species[key]=np.ones_like(self.E)*Species['Specie_IC'][i]
                self.Cp_species[key]=np.ones_like(self.E)*Species['Specie_Cp'][i]
                self.rho_0+=self.rho_species[key]
                i+=1
        
          
    # Calculate and return dimensions of CV
    def CV_dim(self):
        hx=np.zeros_like(self.E)
        
        hx[0]      =0.5*(self.dx[0])
        hx[1:-1]   =0.5*(self.dx[1:-1]+self.dx[:-2])
        hx[-1]     =0.5*(self.dx[-1])
        
        return hx
    
    # Calculate temperature dependent properties
    def calcProp(self):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        Cv=np.zeros_like(self.eta)
        D=copy.deepcopy(self.m_species)
        
        # Species densities and specific heat
        por=[self.porosity,(1-self.porosity)]
        if bool(self.m_species):
#            m_tot=np.zeros_like(self.E) # Use rho to be bulk rho
            for i in range(len(self.species_keys)):
#                self.rho_species[self.species_keys[i]]=\
#                    self.m_species[self.species_keys[i]]/(por[i]*self.CV_vol())
                Cv+=self.rho_species[self.species_keys[i]]*self.Cp_species[self.species_keys[i]]
                rho+=self.rho_species[self.species_keys[i]]
            Cv/=rho
        
        # Calculate properties based on eta or constant
        if type(self.k) is str and (st.find(self.k, 'eta')>=0):
            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
        elif type(self.k) is float:
            k[:]=self.k
        if (type(self.Cv) is str) and (st.find(self.Cv, 'eta')>=0):
            Cv=self.eta*self.Cv1+(1-self.eta)*self.Cv0
        elif type(self.Cv) is float:
            Cv[:]=self.Cv
        if type(self.rho) is str and (st.find(self.rho, 'eta')>=0):
            rho=self.eta*self.rho1+(1-self.eta)*self.rho0
        elif type(self.rho) is float:
            rho[:]=self.rho
        
        # Mass diffusion coefficient; Al, CuO, Al2O3, Cu
        if bool(D):
            for i in self.species_keys:
#                D[i][:,;]=self.Diff.get_Diff(300,i)
                D[i][:]=0
        
        
        return k, rho, Cv, D
    
    # Calculate temperature from energy
    def TempFromConserv(self):
        k,rho,Cv,D=self.calcProp()
        return self.E/Cv/rho
