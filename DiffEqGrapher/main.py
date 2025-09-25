import RungeKutta as rk
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sigfig
from sigfig import round
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({'font.size': 12})

#Physics related variables
iangle = 45. * np.pi / 180.
thrust0 = 933000.
thrust1 = 267700.
thrust2 = 152000.
thrustf = 0.
m0 = 35300.
m1 = 12223.
m2 = 5191.0
mf = 1591.0
f0 = 346.42
f1 = 103.95
f2 = 53.333
ff = 0.
k0 = 1.0062
k1 = 0.63821
k2 = 0.62864
kf = 0.27940
g= 9.8067

#time related variables
t0 = 0.
vx0 = 0.
sx0 = 0.
vy0 = 0.
sy0 = 0.
ax0 = thrust0 / m0 * np.cos(iangle)
ay0 = thrust0 / m0 * np.sin(iangle) - g

dt = 0.01
t_end = 907.08
n_steps = int(round((t_end-t0)/dt))    # number of timesteps
boost1 = 60.
boost2 = 120.
boost3 = 180.

#initiating arrays
vx_arr = np.zeros(n_steps + 1)   # create an array of zeros for Y
sx_arr = np.zeros(n_steps + 1) 
vy_arr = np.zeros(n_steps + 1)   # create an array of zeros for Y
sy_arr = np.zeros(n_steps + 1)   # create an array of zeros for P
ax_arr = np.zeros(n_steps + 1)   # create an array of zeros for P
ay_arr = np.zeros(n_steps + 1)   # create an array of zeros for P
z_arr = np.zeros(n_steps + 1)

t_arr = np.zeros(n_steps + 1)   # create an array of zeros for t
t_arr[0] = t0              # add starttime to array
vx_arr[0] = vx0              # add initial value of Y to array
sx_arr[0] = sx0 
vy_arr[0] = vy0              # add initial value of Y to array
sy_arr[0] = sy0 
ax_arr[0] = ax0 
ay_arr[0] = ay0 

#Differential equations
def dvxdt(t, vx, vy, sx, sy):
   #helper functions
   drag = 0
   mass = 0
   thrust = 0
   if(t < boost1 and t >= 0):
      mass = m0 - f0 * t
      thrust = thrust0
      drag = k0 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vx
   elif(t < boost2 and t >= boost1):
      mass = m1 - f1 * (t - boost1)
      thrust = thrust1
      drag = k1 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vx
   elif(t < boost3 and t >= boost2):
      mass = m2 - f2 * (t - boost2)
      thrust = thrust2
      drag = k2 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vx
   else:
      mass = mf
      thrust = thrustf
      drag = kf * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vx
   
   angle = iangle
   if(t > 0):
      angle = iangle

   #differential equation
   return (thrust * np.cos(angle) - drag) / (mass)

def dvydt(t, vx, vy, sx, sy):
   #helper functions
   gravity = g * (6378100 / (6378100 + sy))**2
   drag = 0
   mass = 0
   thrust = 0
   if(t < boost1 and t >= 0):
      mass = m0 - f0 * t
      thrust = thrust0
      drag = k0 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vy
   elif(t < boost2 and t >= boost1):
      mass = m1 - f1 * (t - boost1)
      thrust = thrust1
      drag = k1 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vy
   elif(t < boost3 and t >= boost2):
      mass = m2 - f2 * (t - boost2)
      thrust = thrust2
      drag = k2 * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vy
   else:
      mass = mf
      thrust = thrustf
      drag = kf * np.e**(-sy/10400) * np.sqrt((vx*vx)+(vy*vy)) * vy

   angle = iangle
   if(t > 0):
      angle = iangle

   #diff eq
   return (thrust * np.sin(angle) - drag) / (mass) - gravity

def dsxdt(t, vx, vy, sx, sy):
   return vx

def dsydt(t, vx, vy, sx, sy):
   return vy

#calculating values
tmax = 0
ymax = 0
for i in range (1, n_steps + 1): 
   t = t_arr[i-1]
   vx = vx_arr[i-1]
   sx = sx_arr[i-1]
   vy = vy_arr[i-1]
   sy = sy_arr[i-1]
   ax = ax_arr[i-1]
   ay = ay_arr[i-1]
   vx_arr[i], vy_arr[i], sx_arr[i], sy_arr[i], ax_arr[i], ay_arr[i] = rk.RK4_4(t, vx, vy, sx, sy, dt, dvxdt, dvydt, dsxdt, dsydt)
   t_arr[i] = t + dt
   if(sy_arr[i-1] > 0 and sy_arr[i] <= 0):
      print("Flight duration: {duration}".format(duration=t_arr[i]))
      print("Flight range: {range}".format(range=sx_arr[i]))
   if(vy_arr[i-1] > 0 and vy_arr[i] <= 0):
      print("Max altitude time: {time}".format(time=t_arr[i]))
      print("Max altitude: {max}".format(max=sy_arr[i]))
      tmax = t_arr[i]
      ymax = sy_arr[i]

print("End vert displacement: {end}".format(end=sy_arr[-1]))
print("End hor displacement: {end}".format(end=sx_arr[-1]))
print("End vert velocity: {end}".format(end=vy_arr[-1]))
print("End hor velocity: {end}".format(end=vx_arr[-1]))
print("End vert acceleration: {end}".format(end=ay_arr[-1]))
print("End hor acceleration: {end}".format(end=ax_arr[-1]))




#plotting values
fig = plt.figure()
t0 = int(t0 / dt)
t1 = int(boost1 / dt)
t2 = int(boost2 / dt)
t3 = int(boost3 / dt)
tf = int(t_end / dt)
tm = int(tmax / dt)
#arrays to graph
xarr = sx_arr
yarr = sy_arr
zarr = z_arr
laserx = [100000, 1658100]
lasery = [2000000, 623290]
laserz = [-300000, 0]
#relevant points
ypoints = [yarr[tm], lasery[0]]
xpoints = [xarr[tm], laserx[0]]
zpoints = [zarr[tm], laserz[0]]

sizes = [200, 200]
colors = ['red', 'black']

'''ax = plt.axes(projection ='3d')
ax.plot(xarr[t0:t1],zarr[t0:t1],yarr[t0:t1], linewidth = 4, label = 'Stage 1', color = 'g')
ax.plot(xarr[t1:t2],zarr[t1:t2],yarr[t1:t2], linewidth = 4, label = 'Stage 1', color = 'b')
ax.plot(xarr[t2:t3],zarr[t2:t3],yarr[t2:t3], linewidth = 4, label = 'Stage 1', color = 'y')
ax.plot(xarr[t3:tf],zarr[t3:tf],yarr[t3:tf], linewidth = 4, label = 'Stage 1', color = 'r')
ax.plot(laserx, laserz, lasery, linewidth=4, label = 'Laser Path', color='black')
ax.scatter(xpoints, zpoints, ypoints, sizes=sizes, c=colors)
ax.set_title('3D Flight Path', fontsize = 12)    # add some title to your plot
ax.set_xlabel('Hor Disp (in m)', fontsize = 12)
ax.set_ylabel('Depth Disp (in m)', fontsize = 12)
ax.set_zlabel('Vert Disp (in m)', fontsize = 12)
plt.show()'''

# create figure
#plt.scatter(xpoints, ypoints, sizes=sizes, c=colors)
plt.plot(xarr[t0 : t1], yarr[t0 : t1], linewidth = 4, label = 'Stage 1', color = 'g') 
plt.plot(xarr[t1 : t2], yarr[t1 : t2], linewidth = 4, label = 'Stage 2', color = 'b') 
plt.plot(xarr[t2 : t3], yarr[t2 : t3], linewidth = 4, label = 'Stage 3', color = 'y') 
plt.plot(xarr[t3 : ], yarr[t3 :], linewidth = 4, label = 'Free Fall', color = 'r') 


for i in range(len(xpoints)):
   x, y = round(xpoints[i], sigfigs=5), round(ypoints[i], sigfigs=5)
   plt.annotate(xy=(xpoints[i] - 3, ypoints[i] + 100), text="({x}, {y})".format(x=x, y=y))
plt.plot(t_arr[t3 :], vy_arr[t3 :], linewidth = 4, label = 's', color='r') 


   # plot Y to t 
plt.title('Vertical Displacement against time', fontsize = 12)    # add some title to your plot
plt.xlabel('Time (in seconds)', fontsize = 12)
plt.ylabel('Vert Disp (in m)', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)                        # show grid
#plt.axis([t0, t_end, 0, 50])     # show axes measures
plt.legend()
plt.show()