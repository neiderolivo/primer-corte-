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
k=2
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



pendulo1=sphere(pos=vector(xp,yp,zp),radius=r,color=color.yellow,make_trail=True)#,trail_type="points")
pendulo2=sphere(pos=vector(xp+l,yp,zp),radius=r,color=color.red,make_trail=True)#,trail_type="points")
cuerda1=curve(vector(0,0,0),pendulo1.pos)
cuerda2=curve(vector(l,0,0),pendulo2.pos)
base=box(pos=vector(l/2,0,0),size=vector(l,0.1,0.1),color=color.magenta)
spring=helix(pos=pendulo1.pos,axis=pendulo2.pos-pendulo1.pos,radius=0.3,constant=k,thickness=0.1,coils=6,color=color.white)
ti=0
while ti<tf:
    sleep(0.01)
    pendulo1.pos=vector(l*sol[ti,0],yp,zp)
    pendulo2.pos=vector(l*sol[ti,2]+l,yp,zp)
    cuerda1=curve(vector(0,0,0),pendulo1.pos)
    cuerda2=curve(vector(l,0,0),pendulo2.pos)
    spring.pos=pendulo1.pos
    spring.axis=pendulo2.pos-spring.pos
    ti = ti + 1
