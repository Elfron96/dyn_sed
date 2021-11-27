# -*- coding: utf-8 -*-
"""
## 1DV MODEL for hydrodynamics and sediment transport
#---------------------------------------------------- 
F. Dufois & H. Muller - October 2018 
Converted from a Matlab code created by R Verney et P. Le Hir ( october 2014 - november 2016)
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import output

plt.close("all")

start_time = time.process_time()

#Generic Parameters

# physical parameters
# -------------------
rhow=1025      # volumic mass (density) (kg/m3)
muw=0.001      # dynamic viscosity (kg/m/s)
nuw=muw/rhow   # kinematic viscosity (m2/s)
vk=0.41        # Von Karman constant
grav=9.81       # gravity (m2/s)

# numerical parameters : time management
# --------------------------------------
dtinit=.01      # initial time step (s)
dtmax=10.       # maximum time step (s)
dts=60.         # output time step (s)
ts=0           # output time counter (initialized to 0)(s)
tfin=50000.    # total running time (s)

nts=np.int(np.ceil(tfin/dts))   # number of output time steps

# hydrodynamic parameters
# -----------------------
slope=0.00005   # free surface slope (0.00005 or 0.00001, or 0.0001 /0.004 for waves)
forcing=slope*grav   # if *rhow, becomes a pressure gradient
z0_val=0.001  # bed roughness (m)


# sediment parameters
wsc=0.01      # Settling velocity if constant (m/s) 
csed=500       # Bed sediment concentration (kg/m3)
cinit=0.1        # initial sediment concentration in suspension (kg/m3)


## grid definition (vertical / "z")
# ---------------------------------
nz=40       # number of layers (nz U&C layers, nz+1 kz interfaces)
# z is located at the center of the layer
# layer height (for linear grid spacing)
dz=.05      # vertical grid spacing, constant (m)

z=np.zeros(nz)
zi=np.zeros(nz+1)

zi[0]=0
for k in range(1,nz+1):  
    zi[k]=zi[k-1]+dz 
       
z=(zi[0:nz]+zi[1:nz+1])/2 # interface (kz) height


## tables declaration and initialization
# -------------------------------------

# output results
dt_s=np.zeros(nts)          # output time step (s)
s=0                        # output counter (initialized to 0)

s=0                     # output counter (initialized to 0)
# netcdf output
# Creation of Model Output object
out = output.ModelOutput()
l_save=False

if l_save:
    out.zi=zi
    out.z=z
    out.create_outputfile('out.nc')
    out.create_variables()


out.dt_s=np.zeros(nts)          # output time step (s)
out.temps_s=np.zeros(nts)       # output elpased time from t0 (s)
out.u_s=np.zeros((nts,nz))          # output current (m/s, >0 ou <0)
out.tenfon_s=np.zeros(nts)       # output bed shear stress (N/m2)
out.kz_s=np.zeros((nts,nz+1))       # diffusion turbulente





# computed variables (hydrodynamics) at the middle of the layer
u=np.zeros(nz)              # current (m/s)
ua=np.zeros(nz)             # current (m/s) at t-1

# computed variables (hydrodynamics) at the interfaces
kz=np.zeros(nz+1)          # vertical turbulent viscosity (m2/s), initialized to 0
kdudz=np.zeros(nz+1)
lgm=np.zeros(nz+1)            # mixing length

ustar=0
kzmax=0                   # maximum turbulent viscosity (m2/s), initialized to 0

tt=0                       # elapsed time, incrementing during simulation (s)
dt=dtinit                  # time step (s), initialised as dtinit



 
# turbulent diffusivity / here we chose it equal to turbulent viscosity




## calculation loop
# --------------------------------------------------------
# --------------------------------------------------------

while tt+dt<tfin:

    ## computation of viscosity (kzu) and current velocity (u)
    # --------------------------------------------------------
    
    ua[:]=u[:]                  # copying u at t-1 to ua


    
## Hydrodynamics
# --------------------------------------------------------
    
    # mixing length and turbulente diffusivity (calculated at layer centers)
    # z(nz) is the total water depth
    
    lgm[1:nz] =  vk * zi[1:nz] * np.sqrt((1-(zi[1:nz]/(dz*nz))))
    
    kz[1:nz]= nuw + (lgm[1:nz]**2 *(abs(ua[1:nz]-ua[0:nz-1])/dz) )
    

    #   Time step update
    #   dynamic update of dt to account for stability criteria
    kzmax=max(kz)
    dt= dtinit# stability criterium when explicit formulation of diffusion process
    
    
    tt=tt+dt
    # Momentum equation
    forcingatt=forcing* 1*np.cos(np.pi*2/(12.4*3600) * tt)
 
    # computing diffusive fluxes, and adding boundary conditions 
    # -------------------------------------
    kdudz[1:nz] = kz[1:nz] * ((ua[1:nz]-ua[0:nz-1])/dz)
    
    # surface boundary condition    
    kdudz[nz]= 0 
    
    # bottom boundary condition   
    ustar= vk*ua[0]/np.log2(dz/(2*z0_val))
    kdudz[0]= ustar  # boundary condition (kz*du/dz=tau/rho=ustar**2)

     
    # solving momentum equation
    # -------------------------------------   

    u[0:nz] = ua[0:nz] + (forcingatt + ((kdudz[1:nz+1]-kdudz[0:nz])/dz))*dt
 
    # bed shear stress
    ustar= vk*ua[0]/np.log2(dz/(2*z0_val)) # boundary condition
    
    tenfon=rhow*ustar**2

    
    if tt>ts and tt<tfin:   
        out.u_s[s,:]=u
        out.dt_s[s]=dt
        out.temps_s[s]=tt
        out.tenfon_s[s]=tenfon
        out.kz_s[s,:]=kz      
         

      
        progress=100*np.float(s+1)/np.float(nts+1)
        print("calculation # completed:",np.floor(progress),"% dt=",dt," u_star=",ustar," m/s")


        if l_save:
            out.write_variables(s,tt)
        
        s=s+1
        ts=dts+ts

if l_save:      
    out.nc.close()        

#
#
print(time.process_time() - start_time, "seconds")
 
  # Figures

fig=plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(out.temps_s/86400,out.tenfon_s)

ax.set_xlabel('time (hours)')
ax.set_ylabel('BSS (Pa)')
plt.show()

fig=plt.figure(2)
ax = fig.add_subplot(111)
for kk in range(10):
   ax.plot(out.u_s[-kk-1,:],np.log10(z))
 
    
ax.set_xlabel('u (m/s)')
ax.set_ylabel('log10(z)')
plt.title("Last 10 velocity profiles")
plt.show()


fig=plt.figure(3)
ax = fig.add_subplot(311)
ax.plot(out.temps_s/86400,out.tenfon_s)
plt.title("BSS")

ax = fig.add_subplot(312)
cax=ax.pcolor(out.temps_s/86400,z,np.transpose(out.u_s))
cax.set_clim(-1,1)
fig.colorbar(cax)
plt.title("U(m/s)")


fig=plt.figure(5)
ax = fig.add_subplot(111)
cax=ax.pcolor(out.temps_s/86400,zi,np.transpose(out.kz_s))
cax.set_clim(0,0.005)
fig.colorbar(cax)
plt.title("Kz")
plt.show()
