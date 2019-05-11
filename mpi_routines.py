# -*- coding: utf-8 -*-
"""
######################################################
#             1D Heat Conduction Solver              #
#              Created by J. Mark Epps               #
#          Part of Masters Thesis at UW 2018-2020    #
######################################################

This file contains the MPI routines:
    -
    
Features:
    -Ignition condition met, will change north BC to that of right BC
    -Saves temperature and reaction data (.npy) depending on input file 
    settings

"""

import numpy as np
import string as st

class MPI_comms():
    def __init__(self, comm, rank, size, Sources, Species):
        self.comm=comm
        self.rank=rank
        self.size=size
        self.Sources=Sources
        self.Species=Species
        
    # Function to split global array to processes
    # Use for MPI_discretize and restart
    def split_var(self, var_global, domain):
        var_local=np.zeros(2)
        # Far left domain
        if self.rank==0:
            var_local=var_global[:domain.Nx+1]
        # Far right domain
        elif self.rank==(self.size-1):
            var_local=var_global[self.rank*domain.Nx-1:]
        # Interior domain
        else:
            var_local=var_global[self.rank*domain.Nx-1:(self.rank+1)*domain.Nx+1]
        
        return var_local
    
    # MPI discretization routine
    def MPI_discretize(self, domain):
        if domain.Nx%self.size!=0:
            return 1
        domain.Nx/=self.size
        
        # Divide global variables
        domain.X=self.split_var(domain.X, domain)
        domain.dx=self.split_var(domain.dx, domain)
        domain.E=self.split_var(domain.E, domain)
        
        # Identify neighboring processes
        domain.proc_left=self.rank-1
        domain.proc_right=self.rank+1
        if self.rank==(self.size-1):
            domain.proc_right=-1
        
        return 0
    
    # Update ghost nodes for processes
    def update_ghosts(self, domain):
        # Send to the left, receive from the right
        a=np.ones(1)*domain.E[-1]
        self.comm.Send(domain.E[1], dest=domain.proc_left)
        self.comm.Recv(a, source=domain.proc_right)
        domain.E[-1]=a
        # Send to the right, receive from the left
        a=np.ones(1)*domain.E[0]
        self.comm.Send(domain.E[-2], dest=domain.proc_right)
        self.comm.Recv(a, source=domain.proc_left)
        domain.E[0]=a
        
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            # Send to the left, receive from the right
            a=np.ones(1)*domain.eta[-1]
            self.comm.Send(domain.eta[1], dest=domain.proc_left)
            self.comm.Recv(a, source=domain.proc_right)
            domain.eta[-1]=a
            # Send to the right, receive from the left
            a=np.ones(1)*domain.eta[0]
            self.comm.Send(domain.eta[-2], dest=domain.proc_right)
            self.comm.Recv(a, source=domain.proc_left)
            domain.eta[0]=a
        if bool(self.Species):
            # Send to the left, receive from the right
            self.comm.Send(domain.P[1], dest=domain.proc_left)
            a=np.ones(1)*domain.P[-1]
            self.comm.Recv(a, source=domain.proc_right)
            domain.P[-1]=a
            # Send to the right, receive from the left
            self.comm.Send(domain.P[-2], dest=domain.proc_right)
            a=np.ones(1)*domain.P[0]
            self.comm.Recv(a, source=domain.proc_left)
            domain.P[0]=a
            for i in self.Species['keys']:
                # Send to the left, receive from the right
                self.comm.Send(domain.m_species[i][1], dest=domain.proc_left)
                a=np.ones(1)*domain.m_species[i][-1]
                self.comm.Recv(a, source=domain.proc_right)
                domain.m_species[i][-1]=a
                # Send to the right, receive from the left
                self.comm.Send(domain.m_species[i][-2], dest=domain.proc_right)
                a=np.ones(1)*domain.m_species[i][0]
                self.comm.Recv(a, source=domain.proc_left)
                domain.m_species[i][0]=a
                
    # General function to compile a variable from all processes
    def compile_var(self, var, Domain):
        var_global=var[:-1].copy()
        if self.rank==0:
            for i in range(self.size-1):
                len_arr=self.comm.recv(source=i+1)
                dat=np.empty(len_arr)
                self.comm.Recv(dat, source=i+1)
                var_global=np.block([var_global, dat])
        elif (Domain.proc_left>=0) and (Domain.proc_right>=0):
            len_arr=len(var)-2
            self.comm.send(len_arr, dest=0)
            self.comm.Send(var[1:-1], dest=0)
        else:
            len_arr=len(var)-1
            self.comm.send(len_arr, dest=0)
            self.comm.Send(var[1:], dest=0)
        len_arr=self.comm.bcast(len(var_global), root=0)
        if self.rank!=0:
            var_global=np.empty(len_arr)
        self.comm.Bcast(var_global, root=0)
        return var_global
        
    # Function to save data to npy files
    def save_data(self, Domain, time, vol):
        T=self.compile_var(Domain.TempFromConserv(vol), Domain)
        np.save('T_'+time, T, False)
        # Kim source term
        if st.find(self.Sources['Source_Kim'],'True')>=0:
            eta=self.compile_var(Domain.eta, Domain)
            np.save('eta_'+time, eta, False)
        if bool(self.Species):
            P=self.compile_var(Domain.P, Domain)
            np.save('P_'+time, P, False)
            for i in self.Species['keys']:
                m_i=self.compile_var(Domain.m_species[i], Domain)
                np.save('m_'+i+'_'+time, m_i, False)