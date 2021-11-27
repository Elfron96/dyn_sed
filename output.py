#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 18:12:49 2021

@author: fdufois
"""

# -*- coding: utf-8 -*-

import numpy as np
import netCDF4


class ModelOutput:
    """
    Main class for model output
    """
    
    def create_outputfile(self,fileout):
            """ Create a netCDF4 file for output time step
            """
            self.nc = netCDF4.Dataset(fileout, 'w')
            # Create dimensions
            self.nc.createDimension('time', None)
            self.nc.createDimension('z',self.z.shape[0])
            self.nc.createDimension('zi',self.zi.shape[0])
    
            # Create dimension variables
            self.nc.createVariable('time', 'd', ('time',))
            self.nc.variables['time'].units = 'seconds'
            self.nc.variables['time'].long_name = 'time in seconds'

           
            self.nc.createVariable('z', 'f', ('z',))
            self.nc.variables['z'].units = ''
            self.nc.variables['z'].long_name = 'z level at u,C points'
            self.nc.variables['z'][:] = self.z[:]
    
            self.nc.createVariable('zi', 'f', ('zi',))
            self.nc.variables['zi'].units = ''
            self.nc.variables['zi'].long_name = 'z level at kz points'
            self.nc.variables['zi'][:] = self.zi[:]
            

    def create_variables(self,l_sed=False):
            """ Create variables in the netCDF file
            """
            print('CREATE NETCDF VARIABLES')
            
            self.nc.createVariable('dt', 'f', ('time',))
            self.nc.variables['dt'].units = 'seconds'
            self.nc.variables['dt'].long_name = 'time step seconds'    
            
            self.nc.createVariable('tenfon', 'f', ('time',))
            self.nc.variables['tenfon'].units = 'N/m2'
            self.nc.variables['tenfon'].long_name = 'Bottom shear stress'
            
            self.nc.createVariable('kz', 'f', ('time', 'zi'))
            self.nc.variables['kz'].units = 'm2/s'
            self.nc.variables['kz'].long_name = 'turbulent diffusivity'
               
            self.nc.createVariable('u', 'f', ('time', 'z'))
            self.nc.variables['u'].units = 'm/s'
            self.nc.variables['u'].long_name = 'horizontal velocity'
            
            if l_sed:
                self.nc.createVariable('c', 'f', ('time', 'z'))
                self.nc.variables['c'].units = 'kg/m3'
                self.nc.variables['c'].long_name = 'sediment concentration'
     
                self.nc.createVariable('mass', 'f', ('time',))
                self.nc.variables['mass'].units = 'kg'
                self.nc.variables['mass'].long_name = 'Total sediment mass'
    
                self.nc.createVariable('hsed', 'f', ('time',))
                self.nc.variables['hsed'].units = 'm'
                self.nc.variables['hsed'].long_name = 'Bed height'   

            
        
    def write_variables(self,ind_t,time,l_sed=False):
            """ Write variables at each output time step
            """
            print('WRITE VARIABLES in Netcdf file',ind_t,time)   
            self.nc.variables['time'][ind_t] = time
            self.nc.variables['dt'][ind_t] = self.dt_s[ind_t]   
            self.nc.variables['tenfon'][ind_t] = self.tenfon_s[ind_t] 
            self.nc.variables['u'][ind_t,:] = self.u_s[ind_t,:]
            self.nc.variables['kz'][ind_t,:] = self.kz_s[ind_t,:]
            
            if l_sed:           
                self.nc.variables['c'][ind_t,:] = self.c_s[ind_t,:]
                self.nc.variables['mass'][ind_t] = self.mass_s[ind_t] 
                self.nc.variables['hsed'][ind_t] = self.hsed_s[ind_t]  
                
            self.nc.sync()
            
