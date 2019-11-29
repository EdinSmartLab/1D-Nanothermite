# -*- coding: utf-8 -*-
"""
######################################################
#             1D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the solver classes for 2D planar and axisymmetric 
Heat Conduction:
    -Modifies geometry class variables directly
    -Calculate time step based on Fourier number
    -Compute conduction equations
    -Add source terms as needed
    -Applies boundary conditions (can vary along a side)

Features/assumptions:
    -time step based on Fourrier number and local discretizations in x
    -equal node spacing in x
    -thermal properties can vary in space (call from geometry object)
    -Radiation boundary conditions

"""

import numpy as np
import copy
import string as st
import Source_Comb
import BCClasses
from mpi4py import MPI

# 2D solver (Cartesian coordinates)
class OneDimLineSolve():
    def __init__(self, geom_obj, settings, Sources, BCs, solver, size, comm):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx=geom_obj.dx
        self.Fo=settings['Fo']
        self.dt=settings['dt']
        self.conv=settings['Convergence']
        self.countmax=settings['Max_iterations']
        self.rank=geom_obj.rank # MPI info
        self.size=size # MPI info
        self.comm=comm # MPI comms
        self.diff_inter=settings['diff_interpolation']
        self.conv_inter=settings['conv_interpolation']
        
        # Define source terms and pointer to source object here
        self.get_source=Source_Comb.Source_terms(Sources['Ea'], Sources['A0'], Sources['dH'], Sources['gas_gen'])
        self.source_unif=Sources['Source_Uniform']
        self.source_Kim=Sources['Source_Kim']
        self.ign=st.split(Sources['Ignition'], ',')
        self.ign[1]=float(self.ign[1])
        
        # BC class
        self.BCs=BCClasses.BCs(BCs, self.dx)
        # Modify BCs if no process is next to current one
        if self.Domain.proc_left>=0:
            self.BCs.BCs['bc_left_E']=['F', 0.0, (0, -1)]
            self.BCs.BCs['bc_left_P']=['none', 0.0, (0, -1)]
        if self.Domain.proc_right>=0:
            self.BCs.BCs['bc_right_E']=['F', 0.0, (0, -1)]
            self.BCs.BCs['bc_right_P']=['none', 0.0, (0, -1)]
            
    # Time step check with dx, dy, Fo number
    def getdt(self, k, rhoC, h):
        # Stability check for Fourrier number
        if self.time_scheme=='Explicit':
            self.Fo=min(self.Fo, 1.0)
        elif self.Fo=='None':
            self.Fo=1.0
        
        dt=self.Fo*rhoC/k*(h)**2
        return np.amin(dt)
    
    # Interpolation function
    def interpolate(self, k1, k2, func):
        if func=='Linear':
            return 0.5*k1+0.5*k2
        else:
            return 2*k1*k2/(k1+k2)
        
    # Main solver (1 time step)
    def Advance_Soln_Cond(self, nt, t, hx, ign):
#    def Advance_Soln_Cond(self, nt, t, hx):
        max_Y,min_Y=0,1
        # Calculate properties
        T_0, k, rhoC, Cp=self.Domain.calcProp(self.Domain.T_guess)
        if self.dt=='None':
            dt=self.getdt(k, rhoC, hx)
            # Collect all dt from other processes and send minimum
            dt=self.comm.reduce(dt, op=MPI.MIN, root=0)
            dt=self.comm.bcast(dt, root=0)
        else:
            dt=min(self.dt,self.getdt(k, rhoC, hx))
            # Collect all dt from other processes and send minimum
            dt=self.comm.reduce(dt, op=MPI.MIN, root=0)
            dt=self.comm.bcast(dt, root=0)
        
        if (np.isnan(dt)) or (dt<=0):
            return 1, dt, ign
        if self.Domain.rank==0:
            print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        # If strang splitting
        if self.time_scheme=='Strang_split':
            dt_strang=[0.5*dt, dt]
        else:
            dt_strang=[dt]
        
        # Copy needed variables and set pointers to other variables
        E_0=self.Domain.E.copy()
        T_c=T_0.copy()
        if self.Domain.model=='Species':
            rho_0=copy.deepcopy(self.Domain.rho_species)
            rho_spec=copy.deepcopy(self.Domain.rho_species)
            species=self.Domain.species_keys
            mu=self.Domain.mu
            perm=self.Domain.perm
        
        # Beginning of strang splitting routine (2 step process)
        for i in range(len(dt_strang)):
            ###################################################################
            # Calculate source terms
            ###################################################################
            # Source terms
            E_unif,E_kim=0,0
            if i==0:
                if self.source_unif!='None':
                    E_unif      = self.source_unif
                if self.source_Kim=='True' or self.Domain.model=='Species':
                    E_kim, deta =self.get_source.Source_Comb_Kim(self.Domain.rho_0, T_c, self.Domain.eta, dt_strang[i])
            
            ###################################################################
            # Conservation of Mass
            ###################################################################
            if self.Domain.model=='Species':
                
                # Calculate pressure
                self.Domain.P=rho_spec[species[0]]/self.Domain.porosity*self.Domain.R*T_c
                
                # Use Darcy's law to directly calculate the velocities at the faces
                flx=np.zeros_like(self.Domain.P)
                
                # Left face
                flx[1:]+=dt_strang[i]/hx[1:]\
                    *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],self.conv_inter)\
                    *(-self.interpolate(perm[1:], perm[:-1],self.diff_inter)/mu\
                    *(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])
                    
                # Right face
                flx[:-1]-=dt_strang[i]/hx[:-1]\
                    *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1], self.conv_inter)\
                    *(-self.interpolate(perm[1:], perm[:-1], self.diff_inter)/mu\
                    *(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])
                    
                self.Domain.rho_species[species[0]]=rho_0[species[0]]+flx
                self.Domain.rho_species[species[1]]=rho_0[species[1]].copy()
                
                # Source terms
                dm0,dm1=self.get_source.Source_mass(deta, self.Domain.porosity, self.Domain.rho_0)
                self.Domain.rho_species[species[0]]+=dm0*dt_strang[i]
                self.Domain.rho_species[species[1]]-=dm1*dt_strang[i]
                
                # Apply pressure BCs
                self.Domain.P=rho_spec[species[0]]/self.Domain.porosity*self.Domain.R*T_c
                eflx=self.BCs.P(self.Domain.P, self.Domain.R, T_c)
                self.Domain.rho_species[species[0]]+=eflx*self.Domain.porosity
                
                # Check max and min for divergence
                max_Y=max(np.amax(self.Domain.rho_species[species[0]]),\
                          np.amax(self.Domain.rho_species[species[1]]))
                min_Y=min(np.amin(self.Domain.rho_species[species[0]]),\
                          np.amin(self.Domain.rho_species[species[1]]))
                
                # Apply BCs
#                self.BCs.mass(self.Domain.m_species[species[0]], self.Domain.P, Ax, Ay)
            ###################################################################
            # Conservation of Energy
            ###################################################################
            self.Domain.E=E_0.copy()
            # Heat diffusion
                #left faces
            self.Domain.E[1:]   -= dt_strang[i]/hx[1:]\
                        *self.interpolate(k[:-1],k[1:], self.diff_inter)\
                        *(T_c[1:]-T_c[:-1])/self.dx[:-1]
            
                # Right face
            self.Domain.E[:-1] += dt_strang[i]/hx[:-1]\
                        *self.interpolate(k[1:],k[:-1], self.diff_inter)\
                        *(T_c[1:]-T_c[:-1])/self.dx[:-1]
            
            # Source terms
            self.Domain.E +=E_unif*dt_strang[i]
            self.Domain.E +=E_kim *dt_strang[i]
            
            if self.Domain.model=='Species':
                # Porous medium advection
                eflx*=Cp*T_c # Convert mass flux to energy flux
                    # Incoming fluxes (Darcy and diffusion)
                eflx[1:]+=dt_strang[i]/hx[1:]\
                    *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],self.conv_inter)*\
                    (-self.interpolate(perm[1:], perm[:-1],self.diff_inter)/mu\
                    *(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])\
                    *self.interpolate(Cp[1:],Cp[:-1],self.conv_inter)\
                    *self.interpolate(T_c[1:],T_c[:-1],self.conv_inter)
                    
                    # Outgoing fluxes (Darcy and diffusion)
                eflx[:-1]-=dt_strang[i]/hx[:-1]\
                    *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],self.conv_inter)*\
                    (-self.interpolate(perm[1:], perm[:-1], self.diff_inter)/mu\
                    *(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])\
                    *self.interpolate(Cp[1:],Cp[:-1],self.conv_inter)\
                    *self.interpolate(T_c[1:],T_c[:-1],self.conv_inter)
                
                self.Domain.E +=eflx
            
            # Apply boundary conditions
            self.BCs.Energy(self.Domain.E, T_0, dt_strang[i], rhoC, hx)
        
        # Check for ignition
        if ign==0 and self.source_Kim=='True':
            if ((self.ign[0]=='eta' and np.amax(self.Domain.eta)>=self.ign[1])\
                or (self.ign[0]=='Temp' and np.amax(T_c)>=self.ign[1])):
                ign=1
        
        # Save previous temp as initial guess for next time step
        self.Domain.T_guess=T_0.copy()
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amin(self.Domain.E)<=0):
            return 2, dt, ign
        elif (np.amax(self.Domain.eta)>1.0) or (np.amin(self.Domain.eta)<-10**(-9)):
            return 3, dt, ign
        elif self.Domain.model=='Species' and ((min_Y<-10)\
                  or np.isnan(max_Y)):
            return 4, dt, ign
        else:
            return 0, dt, ign
