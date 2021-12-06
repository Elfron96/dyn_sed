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

nc_5 = netCDF4.Dataset('5.nc','r')
nc_6 = netCDF4.Dataset('6.nc','r')
nc_7 = netCDF4.Dataset('7.nc','r')
nc_8 = netCDF4.Dataset('8.nc','r')
nc_9 = netCDF4.Dataset('8.nc','r')

nc_10 = netCDF4.Dataset('10.nc','r')
nc_11 = netCDF4.Dataset('11.nc','r')
nc_12 = netCDF4.Dataset('12.nc','r')
nc_13 = netCDF4.Dataset('13.nc','r')

nc_14 = netCDF4.Dataset('14.nc','r')
nc_15 = netCDF4.Dataset('15.nc','r')
nc_16 = netCDF4.Dataset('16.nc','r')

list_nc = [nc_1,nc_2,nc_3,nc_4]
list_nc = [nc_5,nc_6,nc_7,nc_8,nc_9]
list_nc = [nc_10,nc_11,nc_12,nc_13]
list_nc = [nc_14,nc_15,nc_16]

u_1 = nc_1.variables['u'][:]
u_2 = nc_2.variables['u'][:]
u_3 = nc_3.variables['u'][:]
u_4 = nc_4.variables['u'][:]

# plt.figure()
# plt.plot(np.mean(u_1[50:],axis=0),label='z_0 = 0.001 ; GP = 0.0005 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
# plt.plot(np.mean(u_2[50:],axis=0),label='z_0 = 0.0001 ; GP = 0.0001 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
# plt.plot(np.mean(u_3[50:],axis=0),label='z_0 = 0.001 ; GP = 0.0001 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
# plt.plot(np.mean(u_4[50:],axis=0),label='z_0 = 0.0001 ; GP = 0.0005 ; 5cm resolution, hauteur totale = 2m ; dt = 0.01')
# plt.title('Profils de U avec variations de z_0 et du gradient de pression')
# plt.legend()
# plt.show()

plt.figure()
for i in range(len(list_nc)):
    plt.plot(list_nc[0]['time'][:]/86400,list_nc[i]['tenfon'][:], label='Modèle Z_0 n°{}'.format(i+1))
    plt.ylabel('BSS (Pa)')
    plt.xlabel('Temps')
plt.legend()

plt.figure()
for i in range(len(list_nc)):
    # print(list_nc[i]['kz'].shape,list_nc[0]['time'].shape)
    plt.plot(list_nc[0]['time'][:]/86400,list_nc[i]['kz'][:,5], label='Modèle Z_0 n°{}'.format(i+1))
    plt.ylabel('Kz')
    plt.xlabel('Temps')
plt.legend()
plt.show()

# fig=plt.figure(1)
# ax = fig.subplots(len(list_nc),1)
# print(ax.shape)

# for i in range(len(list_nc)):
#     ax[i].plot(list_nc[i]['time'][:]/86400,list_nc[i]['tenfon'][:])
    

#     ax[i].set_xlabel('time (hours)')
#     ax[i].set_ylabel('BSS (Pa)')
#     ax[i].set_title(str(i+5))


# fig=plt.figure(2)
# ax = fig.subplots(len(list_nc),1)
# for i in range(len(list_nc)):
#     for kk in range(10):
#         ax[i].plot(list_nc[i]['u'][-kk-1,:],np.log10(list_nc[i]['z'][:]))

        
#         ax[i].set_xlabel('u (m/s)')
#         ax[i].set_ylabel('log10(z)')
#         ax[i].set_title("Last 10 velocity profiles on simulation n° "+str(i+5))



# fig=plt.figure(3)
# ax = fig.subplots(len(list_nc),1)
# for i in range(len(list_nc)):
#     ax[i].plot(list_nc[i]['time'][:]/86400,list_nc[i]['tenfon'][:])
#     ax[i].set_title("BSS")

# fig=plt.figure(4)
# ax = fig.subplots(len(list_nc),1)
# for i in range(len(list_nc)):
#     cax=ax[i].pcolor(list_nc[i]['time'][:]/86400,list_nc[i]['z'][:],np.transpose(list_nc[i]['u']))
#     cax.set_clim(-1,1)
#     fig.colorbar(cax)
#     ax[i].set_title("U(m/s) simulation n° "+str(i+5))


# # fig=plt.figure(5)
# # ax = fig.add_subplot(111)
# # cax=ax.pcolor(list_nc[i]['time'][:]/86400,list_nc[i]['zi'][:],np.transpose(list_nc[i]['kz'][:]))
# # cax.set_clim(0,0.005)
# # fig.colorbar(cax)
# # plt.title("Kz")





# plt.show()


nc_1.close()
nc_2.close()
nc_3.close()
nc_4.close()
nc_5.close()
nc_6.close()
nc_7.close()
nc_8.close()
nc_9.close()
nc_10.close()
nc_11.close()
nc_12.close()
nc_13.close()
nc_14.close()
nc_15.close()
nc_16.close()