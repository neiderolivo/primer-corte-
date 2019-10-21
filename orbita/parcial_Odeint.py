"""
simulación de dos pendulos ligados por un resorte
"""

from numpy import array, linspace, shape,matrix,transpose
from math import sin, cos, pi,sqrt,atan
from scipy.integrate import odeint
from vpython import  curve,box,cylinder,helix,sphere, scene, vector, color, arrow, text, sleep

def solucion (f,t,g,l,k,m):
    df1=f[1]
    df2=-(g/l)*f[0]-(k/m)*(f[0]-f[2])
    df3=f[3]
    df4=-(g/l)*f[2]+(k/m)*(f[0]-f[2])
    return array([df1,df2,df3,df4])
l=10
k=3
m=50
g=9.8
thetaI=5*pi/180
thetaIp=0
phiI=15*pi/180
phiIp=0
tI=0
tf=1000
n=1000
t=linspace(tI,tf,n)
condINI=array([thetaI,thetaIp,phiI,phiIp])
sol=odeint(solucion,condINI,t,args=(g,l,k,m))
xp=l*thetaI
yp=-l
zp=0
r=1



pendulo1=sphere(pos=vector(xp,yp,zp),radius=r,color=color.green)
pendulo2=sphere(pos=vector(xp+l,yp,zp),radius=r,color=color.green)
#cuerda1=curve(vector(0,0,0),pendulo1.pos)
#cuerda2=curve(vector(l,0,0),pendulo2.pos)

cuerdas1=cylinder(pos=pendulo2.pos,axis=vector(0,0,0),radius=0.1,color=color.white)
cuerdas2=cylinder(pos=pendulo2.pos,axis=vector(l,0,0),radius=0.1,color=color.white)
base=box(pos=vector(l/2,0,0),size=vector(l+2,0.1,0.1),color=color.orange)
spring=helix(pos=pendulo2.pos,axis=vector(0,0,0),radius=0.3,constant=k,thickness=0.1,coils=10,color=color.white)
ti=0
while ti<tf:
    sleep(0.01)
    pendulo1.pos=vector(l*sol[ti,0],yp,zp)
    pendulo2.pos=vector(l*sol[ti,2]+l,yp,zp)
    cuerdas2.pos=pendulo2.pos
    cuerdas2.axis=vector(l,0,0)-cuerdas2.pos
    cuerdas1.pos=pendulo1.pos
    cuerdas1.axis=vector(0,0,0)-cuerdas1.pos
    spring.pos=pendulo1.pos
    spring.axis=pendulo2.pos-spring.pos
    ti = ti + 1
