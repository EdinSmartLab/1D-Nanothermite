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
    def __init__(self, settings, Species, solver):
        
        self.L=settings['Length']
        self.Nx=settings['Nodes_x']
        self.x=np.zeros(self.Nx)
        self.dx=np.zeros(self.Nx) # NOTE: SIZE MADE TO MATCH REST OF ARRAYS (FOR NOW)
        
        # Variables for conservation equations
        self.E=np.zeros(self.Nx) # Lumped energy
        self.eta=np.zeros_like(self.E) # extent of reaction
        self.P=np.zeros_like(self.E) # pressure
        
        # Species
        self.m_species={}
        self.mu_species={}
        self.mv_species={}
        self.rho_species={}
        self.Cp_species={}
        if bool(Species):
            self.species_keys=Species['keys']
            i=0
            for key in self.species_keys:
                self.m_species[key]=np.zeros_like(self.E)
                self.mu_species[key]=np.zeros_like(self.E)
                self.mv_species[key]=np.zeros_like(self.E)
                self.rho_species[key]=np.ones_like(self.E)*Species['Specie_rho'][i]
                self.Cp_species[key]=np.ones_like(self.E)*Species['Specie_Cp'][i]
                i+=1
        self.m_0=np.zeros_like(self.E)
        
        # Thermal properties
        self.k=settings['k']
        self.rho=settings['rho']
        self.Cv=settings['Cp']
        if type(self.rho) is str:
            line=st.split(self.rho, ',')
            self.rho0=float(line[1])
            self.rho1=float(line[2])
        if type(self.Cv) is str:
            line=st.split(self.Cv, ',')
            self.Cv0=float(line[1])
            self.Cv1=float(line[2])
        if type(self.k) is str:
            line=st.split(self.k, ',')
            self.k0=float(line[1])
            self.k1=float(line[2])
        
        self.Diff=Diff_Coef()
        
        # Biasing options       
        self.xbias=[settings['bias_type_x'], settings['bias_size_x']]
        self.isMeshed=False
        # Other useful calculations (put elsewhere??)
        
        
    # Discretize domain and save dx and dy
    def mesh(self):
        # Discretize x
        if self.xbias[0]=='OneWayUp':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
            print 'One way biasing in x: smallest element at x=%2f'%self.L
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='OneWayDown':
            smallest=self.xbias[1]
            self.dx[:-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,self.Nx-1)
            self.dx[-1]=self.dx[-2]
            print 'One way biasing in x: smallest element at x=0'
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayEnd':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
            print 'Two way biasing in x: smallest elements at x=0 and %2f'%self.L
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        elif self.xbias[0]=='TwoWayMid':
            smallest=self.xbias[1]
            self.dx[:int(self.Nx/2)]=np.linspace(2*self.L/(self.Nx-1)-smallest,smallest,(self.Nx-1)/2)
            self.dx[int(self.Nx/2):-1]=np.linspace(smallest,2*self.L/(self.Nx-1)-smallest,(self.Nx-1)/2)
            self.dx[-1]=self.dx[-2]
            print 'Two way biasing in x: smallest elements around x=%2f'%(self.L/2)
            print 'Element size range: %2f, %2f'%(smallest, 2*self.L/(self.Nx-1)-smallest)
        else:
            self.dx[:]=self.L/(self.Nx-1)
            print 'No biasing schemes specified in x'
        
        for i in range(self.Nx-1):
            self.x[i+1]=self.x[i]+self.dx[i]
        
        self.X=self.x
        
        self.isMeshed=True
    
    # Calculate and return volume of each node
    def CV_vol(self):
        v=np.zeros_like(self.E)
        dx=self.dx
        v[0]      =0.5*(dx[0])
        v[1:-1]   =0.5*(dx[1:-1]+dx[:-2])
        v[-1]     =0.25*(dx[-1])
        
        return v
    
    # Calculate temperature dependent properties
    def calcProp(self):
        k=np.zeros_like(self.eta)
        rho=np.zeros_like(self.eta)
        Cv=np.zeros_like(self.eta)
        D=copy.deepcopy(self.m_species)
        
        # Calculate properties based on eta or constant
        if type(self.k) is str:
            k=(self.eta/self.k1+(1-self.eta)/self.k0)**(-1)
        else:
            k[:]=self.k
        if type(self.Cv) is str:
            Cv=self.eta*self.Cv1+(1-self.eta)*self.Cv0
        else:
            Cv[:]=self.Cv
        if type(self.rho) is str:
            rho=self.eta*self.rho1+(1-self.eta)*self.rho0
        else:
            rho[:]=self.rho
        
        # Mass diffusion coefficient; Al, CuO, Al2O3, Cu
        if bool(D):
            for i in self.species_keys:
#                D[i][:,;]=self.Diff.get_Diff(300,i)
                D[i][:]=0
        
        # Species densities
        por=[0.6,0.4]
        if bool(self.m_species):
            for i in range(len(self.species_keys)):
                self.rho_species[self.species_keys[i]]=\
                    self.m_species[self.species_keys[i]]/(por[i]*self.CV_vol())
                
        return k, rho, Cv, D
    
    # Calculate temperature from energy
    def TempFromConserv(self):
        k,rho,Cv,D=self.calcProp()
        return self.E/Cv/rho/self.CV_vol()
