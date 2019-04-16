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
#import CoolProp.CoolProp as CP
#import temporal_schemes
import Source_Comb
import BCClasses

# 2D solver (Cartesian coordinates)
class OneDimLineSolve():
    def __init__(self, geom_obj, settings, Sources, BCs, solver):
        self.Domain=geom_obj # Geometry object
        self.time_scheme=settings['Time_Scheme']
        self.dx=geom_obj.dx
        self.Fo=settings['Fo']
        self.dt=settings['dt']
        self.conv=settings['Convergence']
        self.countmax=settings['Max_iterations']
        
        # Define source terms and pointer to source object here
        self.get_source=Source_Comb.Source_terms(Sources['Ea'], Sources['A0'], Sources['dH'])
        self.source_unif=Sources['Source_Uniform']
        self.source_Kim=Sources['Source_Kim']
        
        # BC class
        self.BCs=BCClasses.BCs(BCs, self.dx)
        
    # Time step check with dx, dy, Fo number
    def getdt(self, k, rho, Cv, vol):
        # Stability check for Fourrier number
        if self.time_scheme=='Explicit':
            self.Fo=min(self.Fo, 1.0)
        elif self.Fo=='None':
            self.Fo=1.0
        
        dt=self.Fo*rho*Cv/k*(vol)**2
        return np.amin(dt)
    
    # Interpolation function
    def interpolate(self, k1, k2, func):
        if func=='Linear':
            return 0.5*k1+0.5*k2
        else:
            return 2*k1*k2/(k1+k2)
        
    # Main solver (1 time step)
    def Advance_Soln_Cond(self, nt, t, vol, Ax):
        max_Y,min_Y=0,1
        # Calculate properties
        k, rho, Cv, D=self.Domain.calcProp()
        mu=10**(-5)
        perm=10**(-11)
        
        if self.dt=='None':
            dt=self.getdt(k, rho, Cv, vol)
        else:
            dt=min(self.dt,self.getdt(k, rho, Cv, vol))
            
        if (np.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1, dt
        print 'Time step %i, Step size=%.7f, Time elapsed=%f;'%(nt+1,dt, t+dt)
        
        # Copy needed variables and set pointers to other variables
        T_c=self.Domain.TempFromConserv()
        if bool(self.Domain.m_species):
            m_c=copy.deepcopy(self.Domain.m_species)
            rho_spec=self.Domain.rho_species
            species=self.Domain.species_keys
            Cp_spec=self.Domain.Cp_species
            mu_c=copy.deepcopy(self.Domain.mu_species)
            mv_c=copy.deepcopy(self.Domain.mv_species)
            
            # Velocity
#            u=mu_c[species[0]]/m_c[species[0]]
#            v=mv_c[species[0]]/m_c[species[0]]
        
        ###################################################################
        # Calculate source and Porous medium terms
        ###################################################################
        # Source terms
        E_unif,E_kim=0,0
        if self.source_unif!='None':
            E_unif      = self.get_source.Source_Uniform(self.source_unif, vol)
        if self.source_Kim=='True':
#            self.Domain.eta=self.Domain.m_species[:,:,2]/0.25
            E_kim, deta =self.get_source.Source_Comb_Kim(rho, T_c, self.Domain.eta, vol, dt)
#            E_kim, deta =self.get_source.Source_Comb_Umbrajkar(rho, T_c, self.Domain.eta, self.Domain.CV_vol(), dt)
        
        
        ###################################################################
        # Conservation of Mass
        ###################################################################
        if bool(self.Domain.m_species):
            # Adjust pressure
            print '     Gas mass: %f, %f'%(np.amax(self.Domain.m_species['g'])*10**6,np.amin(self.Domain.m_species['g'])*10**6)
            print '     Gas density: %f, %f'%(np.amax(rho_spec['g']),np.amin(rho_spec['g']))
            self.Domain.P=self.Domain.m_species['g']/102*1000*8.314*300/(0.6*vol)
#            self.BCs.P(self.Domain.P)
            print '     Pressure: %f, %f'%(np.amax(self.Domain.P),np.amin(self.Domain.P))
            
            # Use Darcy's law to directly calculate the velocities at the faces
            # Ingoing fluxes
            flx=np.zeros_like(self.Domain.P)
            
            # Left face
            flx[1:]+=Ax[1:]*dt\
                *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])
    #            self.interpolate(u[:,1:], u[:,:-1], 'Linear')
            
            # Right face
            flx[:-1]-=Ax[:-1]*dt\
                *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1], 'Linear')*\
                (-perm/mu*(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])
    #            self.interpolate(u[:,1:], u[:,:-1], 'Linear')
            
            
            print '    Gas fluxes in x: %f, %f'%(np.amax(flx)*10**(9),np.amin(flx)*10**(9))
            
            self.Domain.m_species[species[0]]+=flx
            
            # Source terms
    #        dm=deta*dt*(m_c[species[0]]+m_c[species[1]])
#            dm=np.zeros_like(deta)
            dm=deta*dt*(self.Domain.m_0)
#            dm[dm<10**(-9)]=0
#            print '     Mass generated: %f, %f'%(np.amax(dm)*10**(9),np.amin(dm)*10**(9))
    #        (m_c[species[0]]+m_c[species[1]])
            self.Domain.m_species[species[0]]+=dm
            self.Domain.m_species[species[1]]-=dm
                    
            max_Y=max(np.amax(self.Domain.m_species[species[0]]),\
                      np.amax(self.Domain.m_species[species[1]]))
            min_Y=min(np.amin(self.Domain.m_species[species[0]]),\
                      np.amin(self.Domain.m_species[species[1]]))
            
            # Apply BCs
#            self.BCs.mass(self.Domain.m_species[species[0]], self.Domain.P, Ax, Ay, vol)
        
        ###################################################################
        # Conservation of Momentum (x direction; gas)
        ###################################################################
        # Fluxes
#        self.Domain.mu_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(mu_c[species[0]][:,1:]+mu_c[species[0]][:,:-1])*self.interpolate(u[:,1:], u[:,:-1], 'Linear')
#        self.Domain.mu_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(mu_c[species[0]][1:,:]+mu_c[species[0]][:-1,:])*self.interpolate(v[1:,:], v[:-1,:], 'Linear')
#        self.Domain.mu_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(mu_c[species[0]][:,1:]+mu_c[species[0]][:,:-1])*self.interpolate(u[:,1:], u[:,:-1], 'Linear')
#        self.Domain.mu_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(mu_c[species[0]][1:,:]+mu_c[species[0]][:-1,:])*self.interpolate(v[1:,:], v[:-1,:], 'Linear')
#        
#        # Pressure
#        self.Domain.mu_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(self.Domain.P[:,1:]+self.Domain.P[:,:-1])
#        self.Domain.mu_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(self.Domain.P[:,1:]+self.Domain.P[:,:-1])
#                
#        # Porous medium losses
#        self.Domain.mu_species[species[0]]-=mu/perm*u*vol*dt
#        
#        ###################################################################
#        # Conservation of Momentum (y direction; gas)
#        ###################################################################
#        # Fluxes
#        self.Domain.mv_species[species[0]][:,1:]+=Ax[:,1:]*dt\
#            *0.5*(mv_c[species[0]][:,1:]+mv_c[species[0]][:,:-1])*0.5*(u[:,1:]+u[:,:-1])
#        self.Domain.mv_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(mv_c[species[0]][1:,:]+mv_c[species[0]][:-1,:])*0.5*(v[1:,:]+v[:-1,:])
#        self.Domain.mv_species[species[0]][:,:-1]-=Ax[:,:-1]*dt\
#            *0.5*(mv_c[species[0]][:,1:]+mv_c[species[0]][:,:-1])*0.5*(u[:,1:]+u[:,:-1])
#        self.Domain.mv_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(mv_c[species[0]][1:,:]+mv_c[species[0]][:-1,:])*0.5*(v[1:,:]+v[:-1,:])
#        
#        # Pressure
#        self.Domain.mv_species[species[0]][1:,:]+=Ay[1:,:]*dt\
#            *0.5*(self.Domain.P[1:,:]+self.Domain.P[:-1,:])
#        self.Domain.mv_species[species[0]][:-1,:]-=Ay[:-1,:]*dt\
#            *0.5*(self.Domain.P[1:,:]+self.Domain.P[:-1,:])
#                
#        # Porous medium losses
#        self.Domain.mu_species[species[0]]-=mu/perm*v*vol*dt
        
        
        ###################################################################
        # Conservation of species
        ###################################################################
#        if bool(self.Domain.m_species):
#            # Mole ratios
#            mole_ratio={}
#            mole_ratio[species[0]]=-2.0/5 # Al
#            mole_ratio[species[1]]=-3.0/5 # CuO
#            mole_ratio[species[2]]=1.0/4  # Al2O3
#            mole_ratio[species[3]]=3.0/4  # Cu
#            
#            for i in species:
#                # Calculate flux coefficients
#                aW,aE,aS,aN=self.get_Coeff(self.dx,self.dy, dt, rho*D[i], 'Linear')
#                
#                # Diffusion contribution (2nd order central schemes)
#                self.Domain.m_species[i][:,1:]    = aW[:,1:]    * m_c[i][:,:-1]
#                self.Domain.m_species[i][:,0]     = aE[:,0]     * m_c[i][:,1]
#                
#                self.Domain.m_species[i][:,1:-1] += aE[:,1:-1]  * m_c[i][:,2:]
#                self.Domain.m_species[i][1:,:]   += aS[1:,:]    * m_c[i][:-1,:]
#                self.Domain.m_species[i][:-1,:]  += aN[:-1,:]   * m_c[i][1:,:]
#                self.Domain.m_species[i]         -= (aW+aE+aS+aN)*m_c[i]
#            
#                # Species generated/destroyed during reaction
#                self.Domain.m_species[i]+=mole_ratio[i]*deta
#                
#                # Species advected from Porous medium equations [TO BE CONTINUED]
#                
#                
#                # Apply data from previous time step
#                self.Domain.m_species[i]*= dt
#                self.Domain.m_species[i]+= m_c[i]
##                print(self.Domain.m_species[i])
#                # IMPLICITLY MAKING SPECIES FLUX 0 AT BOUNDARIES
#                max_Y=max(np.amax(self.Domain.m_species[i]), max_Y)
#                min_Y=min(np.amin(self.Domain.m_species[i]), min_Y)
#        print(self.Domain.m_species)
        ###################################################################
        # Conservation of Energy
        ###################################################################
        # Heat diffusion
            #left faces
        self.Domain.E[1:]   -= dt*self.interpolate(k[:-1],k[1:], 'Harmonic')\
                    *(T_c[1:]-T_c[:-1])/self.dx[:-1]*Ax[1:]
            # Right face
        self.Domain.E[:-1] += dt*self.interpolate(k[:-1],k[1:], 'Harmonic')\
                    *(T_c[1:]-T_c[:-1])/self.dx[:-1]*Ax[:-1]
        
        # Source terms
        self.Domain.E +=E_unif*dt
        self.Domain.E +=E_kim *dt
        
        if bool(self.Domain.m_species):
            # Porous medium advection
            eflx=np.zeros_like(self.Domain.P)
                # Incoming fluxes
            eflx[1:]+=dt\
                *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])\
                *0.5*(T_c[1:]+T_c[:-1])*0.5*(Cp_spec[species[0]][1:]+Cp_spec[species[0]][:-1])
                # Outgoing fluxes
            eflx[:-1]-=dt\
                *self.interpolate(rho_spec[species[0]][1:],rho_spec[species[0]][:-1],'Linear')*\
                (-perm/mu*(self.Domain.P[1:]-self.Domain.P[:-1])/self.dx[:-1])\
                *0.5*(T_c[1:]+T_c[:-1])*0.5*(Cp_spec[species[0]][1:]+Cp_spec[species[0]][:-1])
            
            print '    Gas energy flux in x: %f, %f'%(np.amax(eflx)*10**(9),np.amin(eflx)*10**(9))
            self.Domain.E +=eflx
#        # Radiation effects
#        self.Domain.T[1:-1,1:-1]+=0.8*5.67*10**(-8)*(T_c[:-2,1:-1]**4+T_c[2:,1:-1]**4+T_c[1:-1,:-2]**4+T_c[1:-1,2:]**4)
        
        # Apply boundary conditions
        self.BCs.Energy(self.Domain.E, T_c, dt, rho, Cv, vol)
        
        ###################################################################
        # Divergence/Convergence checks
        ###################################################################
        if (np.isnan(np.amax(self.Domain.E))) \
        or (np.amin(self.Domain.E)<=0):
            print '***********Divergence detected - energy************'
            return 2, dt
        elif (np.amax(self.Domain.eta)>1.0) or (np.amin(self.Domain.eta)<-10**(-9)):
            print '***********Divergence detected - reaction progress************'
            return 3, dt
        elif bool(self.Domain.m_species) and ((min_Y<-10**(-9))\
                  or np.isnan(max_Y)):
            print '***********Divergence detected - species mass ****************'
            return 4, dt
        else:
            return 0, dt
