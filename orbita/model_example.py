
####################################
#sistema soalr externo
#planetas jupyter, saturno, neptuno, y urano.
#######################################

from numpy import array, linspace, shape,matrix,transpose
from math import sin, cos, pi,sqrt,atan
from scipy.integrate import odeint
from vpython import sphere, scene, vector, color, arrow, text, sleep
scene.range = 70 # m

xp =1.
yp =0.
zp =0.

sleeptime = 0.0001
n_steps = 900
t_start = 0.
t_final = 200.
t_delta = (t_final - t_start) / n_steps
t = linspace(t_start, t_final, n_steps)

#arrow_size = 3

#arrow_x = arrow(pos=vector(0,0,0), axis=vector(arrow_size,0,0), color=color.red)
#arrow_y = arrow(pos=vector(0,0,0), axis=vector(0,arrow_size,0), color=color.green)
#arrow_z = arrow(pos=vector(0,0,0), axis=vector(0,0,arrow_size))
#____________________________________________________________--
QJu= 5.219#UA jupiter
eJu= 0.0484
aJu= QJu/(1+eJu) #UA
bJu=sqrt((1-(eJu*eJu))*(aJu*aJu)) #AU
qJu= 2*aJu - QJu #UA
RJu= 0.2

thes= 45*pi/180.
omes= 0.

def orbitJu (initJu, t):
    dxJu = initJu[0]
    dyJu = initJu[1]
    dthJu = initJu[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rJu=(aJu*(1-eJu**2))/(1+eJu*cos(thes))
    acJu = -GM / (rJu*rJu*rJu)

    dv_xJu = acJu*initJu[2]
    dv_yJu = acJu*initJu[3]

    dtheJu = (sqrt(((dv_xJu*dv_xJu)/((2/rJu)-(1/aJu))) + (dv_yJu*dv_yJu)/((2/rJu)-(1/bJu))))
    return array([dv_xJu, dv_yJu, dxJu, dyJu, dtheJu], float)


initJu = [0., 2.*pi, aJu, 0, 0]
solJu,outodeint = odeint(orbitJu, initJu, t, full_output=True)

vxxJu, vyyJu, xxJu, yyJu, thJu = solJu.T

#============================================
QSa= 9.174#UA Saturno
eSa= 0.0541
aSa= QSa/(1+eSa) #UA
bSa=sqrt((1-(eSa*eSa))*(aSa*aSa)) #AU
qSa= 2*aSa - QSa #UA
RSa= 0.26

thes= 45*pi/180.
omes= 0.

def orbitSa (initSa, t):
    dxSa = initSa[0]
    dySa = initSa[1]
    dthSa = initSa[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rSa=(aSa*(1-eSa**2))/(1+eSa*cos(thes))
    acSa = -GM / (rSa*rSa*rSa)

    dv_xSa = acSa*initSa[2]
    dv_ySa = acSa*initSa[3]

    dtheSa = (sqrt(((dv_xSa*dv_xSa)/((2/rSa)-(1/aSa))) + (dv_ySa*dv_ySa)/((2/rSa)-(1/bSa))))
    return array([dv_xSa, dv_ySa, dxSa, dySa, dtheSa], float)


initSa = [0., 2.*pi, aSa, 0, 0]
solSa,outodeint = odeint(orbitSa, initSa, t, full_output=True)

vxxSa, vyySa, xxSa, yySa, thSA = solSa.T

#_----------------------------------------------------
QUr= 19.23#UA urano
eUr= 0.0472
aUr= QUr/(1+eUr) #UA
bUr=sqrt((1-(eUr*eUr))*(aUr*aUr)) #AU
qUr= 2*aUr - QUr #UA
RUr= 0.26

thes= 45*pi/180.
omes= 0.

def orbitUr (initUr, t):
    dxUr = initUr[0]
    dyUr = initUr[1]
    dthUr = initUr[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rUr=(aUr*(1-eUr**2))/(1+eUr*cos(thes))
    acUr = -GM / (rUr*rUr*rUr)

    dv_xUr = acUr*initUr[2]
    dv_yUr = acUr*initUr[3]

    dtheUr = (sqrt(((dv_xUr*dv_xUr)/((2/rUr)-(1/aUr))) + (dv_yUr*dv_yUr)/((2/rUr)-(1/bUr))))
    return array([dv_xUr, dv_yUr, dxUr, dyUr, dtheUr], float)


initUr = [0., 2.*pi, aUr, 0, 0]
solUr,outodeint = odeint(orbitUr, initUr, t, full_output=True)

vxxUr, vyyUr, xxUr, yyUr, thUr = solUr.T

#======================0
QNe= 30.1#UA neptuno
eNe= 0.0472
aNe= QNe/(1+eNe) #UA
bNe=sqrt((1-(eNe*eNe))*(aNe*aNe)) #AU
qNe= 2*aNe - QNe #UA
RNe= 0.26

thes= 45*pi/180.
omes= 0.

def orbitNe (initNe, t):
    dxNe = initNe[0]
    dyNe = initNe[1]
    dthNe = initNe[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rNe=(aNe*(1-eNe**2))/(1+eNe*cos(thes))
    acNe = -GM / (rNe*rNe*rNe)

    dv_xNe = acNe*initNe[2]
    dv_yNe = acNe*initNe[3]

    dtheNe = (sqrt(((dv_xNe*dv_xNe)/((2/rNe)-(1/aNe))) + (dv_yNe*dv_yNe)/((2/rNe)-(1/bNe))))
    return array([dv_xNe, dv_yNe, dxNe, dyNe, dtheNe], float)


initNe = [0., 2.*pi, aNe, 0, 0]
solNe,outodeint = odeint(orbitNe, initNe, t, full_output=True)

vxxNe, vyyNe, xxNe, yyNe, thNe = solNe.T

#====================================================
Rsol=0.25

prtclJu =sphere(pos=vector(xp,yp,zp),radius=RJu,color=color.orange,make_trail=True,trail_type="points")

prtclSa =sphere(pos=vector(xp,yp,zp),radius=RSa,color=color.white,make_trail=True,trail_type="points")

prtclUr =sphere(pos=vector(xp,yp,zp),radius=RUr,color=color.cyan,make_trail=True,trail_type="points")

prtclNe =sphere(pos=vector(xp,yp,zp),radius=RNe,color=color.green,make_trail=True,trail_type="points")

prtclSOl =sphere(pos=vector(0,0,0),radius=Rsol,color=color.yellow,make_trail=True,trail_type="points")
#------------------------------------------

time_i = 0
t_run = 0
angles=array((1.303,2.489,0.773,1.7770))*pi/180 

MrotJu=matrix([[1,0,0],[0,cos(angles[0]),-sin(angles[0])],[0,sin(angles[0]),cos(angles[0])]])

MrotSa=matrix([[1,0,0],[0,cos(angles[1]),-sin(angles[1])],[0,sin(angles[1]),cos(angles[1])]])

MrotUr=matrix([[1,0,0],[0,cos(angles[2]),-sin(angles[2])],[0,sin(angles[2]),cos(angles[2])]])

MrotNe=matrix([[1,0,0],[0,cos(angles[3]),-sin(angles[3])],[0,sin(angles[3]),cos(angles[3])]])

while t_run<t_final:
    sleep(sleeptime)
    Ju = MrotJu*transpose(matrix([xxJu[time_i],yyJu[time_i], zp]))
    Sa = MrotSa*transpose(matrix([xxSa[time_i],yySa[time_i], zp]))
    Ur = MrotUr*transpose(matrix([xxUr[time_i],yyUr[time_i], zp]))
    Ne = MrotNe*transpose(matrix([xxNe[time_i],yyNe[time_i], zp]))
    prtclJu.pos = vector(Ju[0],Ju[1],Ju[2])
    prtclSa.pos = vector(Sa[0],Sa[1],Sa[2])
    prtclUr.pos = vector(Ur[0],Ur[1],Ur[2])
    prtclNe.pos = vector(Ne[0],Ne[1],Ne[2])
    t_run += t_delta
    time_i += 1
