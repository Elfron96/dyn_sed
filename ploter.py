from typing import List
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from windrose import WindroseAxes
import matplotlib.cm as cm

nc_1 = netCDF4.Dataset('1.nc','r')
nc_2 = netCDF4.Dataset('2.nc','r')
nc_3 = netCDF4.Dataset('3.nc','r')
nc_4 = netCDF4.Dataset('4.nc','r')

list_nc = [nc_1,nc_2,nc_3,nc_4]

u_1 = nc_1.variables['u'][:]
u_2 = nc_2.variables['u'][:]
u_3 = nc_3.variables['u'][:]
u_4 = nc_4.variables['u'][:]

plt.figure()
plt.plot(np.mean(u_1[50:],axis=0),label='z_0 = 0.001 ; GP = 0.0005 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
plt.plot(np.mean(u_2[50:],axis=0),label='z_0 = 0.0001 ; GP = 0.0001 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
plt.plot(np.mean(u_3[50:],axis=0),label='z_0 = 0.001 ; GP = 0.0001 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
plt.plot(np.mean(u_4[50:],axis=0),label='z_0 = 0.0001 ; GP = 0.0005 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
plt.title('Profils de U avec variations de z_0 et du gradient de pression')
plt.legend()
plt.show()




fig=plt.figure(1)
ax = fig.subplots(4,1)
print(ax.shape)

for i in range(len(list_nc)):
    ax[i].plot(list_nc[i]['time'][:]/86400,list_nc[i]['tenfon'][:])
    

    ax[i].set_xlabel('time (hours)')
    ax[i].set_ylabel('BSS (Pa)')
    ax[i].set_title(str(i+1))


fig=plt.figure(2)
ax = fig.subplots(4,1)
for i in range(len(list_nc)):
    for kk in range(10):
        ax[i].plot(list_nc[i]['u'][-kk-1,:],np.log10(list_nc[i]['z'][:]))

        
        ax[i].set_xlabel('u (m/s)')
        ax[i].set_ylabel('log10(z)')
        ax[i].set_title("Last 10 velocity profiles on simulation n° "+str(i+1))



fig=plt.figure(3)
ax = fig.subplots(len(list_nc),1)
for i in range(len(list_nc)):
    ax[i].plot(list_nc[i]['time'][:]/86400,list_nc[i]['tenfon'][:])
    ax[i].set_title("BSS")

fig=plt.figure(4)
ax = fig.subplots(len(list_nc),1)
for i in range(len(list_nc)):
    cax=ax[i].pcolor(list_nc[i]['time'][:]/86400,list_nc[i]['z'][:],np.transpose(list_nc[i]['u']))
    cax.set_clim(-1,1)
    fig.colorbar(cax)
    ax[i].set_title("U(m/s) simulation n° "+str(i+1))


# fig=plt.figure(5)
# ax = fig.add_subplot(111)
# cax=ax.pcolor(list_nc[i]['time'][:]/86400,list_nc[i]['zi'][:],np.transpose(list_nc[i]['kz'][:]))
# cax.set_clim(0,0.005)
# fig.colorbar(cax)
# plt.title("Kz")





plt.show()

nc_1.close()
nc_2.close()
nc_3.close()
nc_4.close()