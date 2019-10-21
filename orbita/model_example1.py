from numpy import array, linspace
from math import sin, cos, pi
from pylab import plot, xlabel, ylabel, show
from scipy.integrate import odeint
from math import pi, sqrt, atan 
from numpy import linspace, shape, array
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from pylab import plot, show, savefig, ylim, xlim, grid, xlabel, ylabel,imshow

from vpython import sphere, scene, vector, color, arrow, text, sleep

arrow_size = 0.1

arrow_x = arrow(pos=vector(0,0,0), axis=vector(arrow_size,0,0), color=color.red)
arrow_y = arrow(pos=vector(0,0,0), axis=vector(0,arrow_size,0), color=color.green)
arrow_z = arrow(pos=vector(0,0,0), axis=vector(0,0,arrow_size))

Q= 1#UA
e= 0.8
a= Q/(1+e) #UA
b= sqrt((1-(e*e))*(a*a)) #AU
q= 2*a - Q #UA
R = 0.1 
thes = 45*pi/180.
omes = 0.

def orbit (init, t):
    dx = init[0]
    dy = init[1]
    dth = init[4]
	            
    GM = 4*pi*pi # UA^3 / yr^2
    r = sqrt(init[2]*init[2] + init[3]*init[3])
    ac = -GM / (r*r*r)

    dv_x = ac*init[2]
    dv_y = ac*init[3]
					    
    dthe = (sqrt(((dv_x*dv_x)/((2/r)-(1/a))) + (dv_y*dv_y)/((2/r)-(1/b)))) 
    #print(dx,dy)
    return array([dv_x, dv_y, dx, dy, dthe], float) 


init = [0., 2.*pi, a, 0, 0]
n_steps = 1000
t_start = 0.
t_final = 15.
t_delta = (t_final - t_start) / n_steps
t = linspace(t_start, t_final, n_steps)


sol,outodeint = odeint(orbit, init, t, full_output=True)

vxx, vyy, xx, yy, th = sol.T



# =====================

scene.range = 2 # m
xp=1
yp=0
zp =0
sleeptime = 0.0001

prtcl = sphere(pos=vector(xp,yp,zp), radius=R, color=color.cyan)

time_i = 0
t_run = 0

#for i in omega:
#    print(i)

while time_i<t_final:
    sleep(sleeptime)
    prtcl.pos = vector( xx[time_i], yy[time_i], time_i)
    t_run=time_i*delta
    time_i+=time_i
