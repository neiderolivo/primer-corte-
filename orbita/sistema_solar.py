from numpy import array, linspace, shape,matrix,transpose
from math import sin, cos, pi,sqrt,atan
from scipy.integrate import odeint
from vpython import sphere, scene, vector, color, arrow, text, sleep
scene.range = 70 # m

xp =1.
yp =0.
zp =0.

sleeptime = 0.0001
n_steps = 1000
t_start = 0.
t_final = 100.
t_delta = (t_final - t_start) / n_steps
t = linspace(t_start, t_final, n_steps)

#arrow_size = 3

#arrow_x = arrow(pos=vector(0,0,0), axis=vector(arrow_size,0,0), color=color.red)
#arrow_y = arrow(pos=vector(0,0,0), axis=vector(0,arrow_size,0), color=color.green)
#arrow_z = arrow(pos=vector(0,0,0), axis=vector(0,0,arrow_size))

QMe= 0.389#UA mercurio
eMe= 0.2056
aMe= QMe/(1+eMe) #UA
bMe=sqrt((1-(eMe*eMe))*(aMe*aMe)) #AU
qMe= 2*aMe - QMe #UA
RMe= 0.2#0.00279

thes= 45*pi/180.
omes= 0.

def orbitMe (initMe, t):
    dxMe = initMe[0]
    dyMe = initMe[1]
    dthMe = initMe[4]
	            
    GM = 4*pi*pi # UA^3 / yr^2
    rMe=(aMe*(1-eMe**2))/(1+eMe*cos(thes))
   # r = sqrt(init[2]*init[2] + init[3]*init[3])
    acMe = -GM / (rMe*rMe*rMe)

    dv_xMe = acMe*initMe[2]
    dv_yMe = acMe*initMe[3]
					    
    dtheMe = (sqrt(((dv_xMe*dv_xMe)/((2/rMe)-(1/aMe))) + (dv_yMe*dv_yMe)/((2/rMe)-(1/bMe)))) 
    #print(dx,dy)
    return array([dv_xMe, dv_yMe, dxMe, dyMe, dtheMe], float) 


initMe = [0., 2.*pi, aMe, 0, 0]

solMe,outodeint = odeint(orbitMe, initMe, t, full_output=True)

vxxMe, vyyMe, xxMe, yyMe, thMe = solMe.T

# ====================
QVe= 0.719#UA Venus
eVe= 0.0068
aVe= QVe/(1+eVe) #UA
bVe=sqrt((1-(eVe*eVe))*(aVe*aVe)) #AU
qVe= 2*aVe - QVe #UA
RVe= 0.2#0.007

thes= 45*pi/180.
omes= 0.

def orbitVe (initVe, t):
    dxVe = initVe[0]
    dyVe = initVe[1]
    dthVe = initVe[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rVe=(aVe*(1-eVe**2))/(1+eVe*cos(thes))
    acVe = -GM / (rVe*rVe*rVe)

    dv_xVe = acVe*initVe[2]
    dv_yVe = acVe*initVe[3]

    dtheVe = (sqrt(((dv_xVe*dv_xVe)/((2/rVe)-(1/aVe))) + (dv_yVe*dv_yVe)/((2/rVe)-(1/bVe))))
    return array([dv_xVe, dv_yVe, dxVe, dyVe, dtheVe], float)


initVe = [0., 2.*pi, aVe, 0, 0]
solVe,outodeint = odeint(orbitVe, initVe, t, full_output=True)

vxxVe, vyyVe, xxVe, yyVe, thVe = solVe.T
# ====================
QTe= 1#UA Tierra
eTe= 0.0167
aTe= QTe/(1+eTe) #UA
bTe=sqrt((1-(eTe*eTe))*(aTe*aTe)) #AU
qTe= 2*aTe - QTe #UA
RTe= 0.2#0.007

thes= 45*pi/180.
omes= 0.

def orbitTe (initTe, t):
    dxTe = initTe[0]
    dyTe = initTe[1]
    dthTe = initTe[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rTe=(aTe*(1-eTe**2))/(1+eTe*cos(thes))
    acTe = -GM / (rTe*rTe*rTe)

    dv_xTe = acTe*initTe[2]
    dv_yTe = acTe*initTe[3]

    dtheTe = (sqrt(((dv_xTe*dv_xTe)/((2/rTe)-(1/aTe))) + (dv_yTe*dv_yTe)/((2/rTe)-(1/bTe))))
    return array([dv_xTe, dv_yTe, dxTe, dyTe, dtheTe], float)


initTe = [0., 2.*pi, aTe, 0, 0]
solTe,outodeint = odeint(orbitTe, initTe, t, full_output=True)

vxxTe, vyyTe, xxTe, yyTe, thTe = solTe.T
#=========================================
QMa= 1.524#UA MArte
eMa= 0.0934
aMa= QMa/(1+eMa) #UA
bMa=sqrt((1-(eMa*eMa))*(aMa*aMa)) #AU
qMa= 2*aMa - QMa #UA
RMa= 0.2#0.0038

thes= 45*pi/180.
omes= 0.

def orbitMa (initMa, t):
    dxMa = initMa[0]
    dyMa = initMa[1]
    dthMa = initMa[4]

    GM = 4*pi*pi # UA^3 / yr^2
    rMa=(aMa*(1-eMa**2))/(1+eMa*cos(thes))
    acMa = -GM / (rMa*rMa*rMa)

    dv_xMa = acMa*initMa[2]
    dv_yMa = acMa*initMa[3]

    dtheMa = (sqrt(((dv_xMa*dv_xMa)/((2/rMa)-(1/aMa))) + (dv_yMa*dv_yMa)/((2/rMa)-(1/bMa))))
    return array([dv_xMa, dv_yMa, dxMa, dyMa, dtheMa], float)


initMa = [0., 2.*pi, aMa, 0, 0]
solMa,outodeint = odeint(orbitMa, initMa, t, full_output=True)

vxxMa, vyyMa, xxMa, yyMa, thMa = solMa.T

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
RSa= 0.2#0.0666

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
RUr= 0.2#0.0666

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
RNe= 0.2#0.0666

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

prtclMe = sphere(pos=vector(xp,yp,zp),radius=RMe,color=color.yellow,make_trail=True,trail_type="points")

prtclVe = sphere(pos=vector(xp,yp,zp),radius=RVe,color=color.green,make_trail=True,trail_type="points")

prtclTe =sphere(pos=vector(xp,yp,zp),radius=RTe,color=color.blue,make_trail=True,trail_type="points")

prtclMa =sphere(pos=vector(xp,yp,zp),radius=RMa,color=color.red,make_trail=True,trail_type="points")

prtclJu =sphere(pos=vector(xp,yp,zp),radius=RJu,color=color.orange,make_trail=True,trail_type="points")

prtclSa =sphere(pos=vector(xp,yp,zp),radius=RSa,color=color.white,make_trail=True,trail_type="points")

prtclUr =sphere(pos=vector(xp,yp,zp),radius=RUr,color=color.cyan,make_trail=True,trail_type="points")

prtclNe =sphere(pos=vector(xp,yp,zp),radius=RNe,color=color.magenta,make_trail=True,trail_type="points")

#------------------------------------------

time_i = 0
t_run = 0
angles=array((7,3.39,0,1.85,1.31,2.49,0.773,1.7770))*pi/180 

MrotMe=matrix([[1,0,0],[0,cos(angles[0]),-sin(angles[0])],[0,sin(angles[0]),cos(angles[0])]])

MrotVe=matrix([[1,0,0],[0,cos(angles[1]),-sin(angles[1])],[0,sin(angles[1]),cos(angles[1])]])

MrotTe=matrix([[1,0,0],[0,cos(angles[2]),-sin(angles[2])],[0,sin(angles[2]),cos(angles[2])]])

MrotMa=matrix([[1,0,0],[0,cos(angles[3]),-sin(angles[3])],[0,sin(angles[3]),cos(angles[3])]])

MrotJu=matrix([[1,0,0],[0,cos(angles[4]),-sin(angles[4])],[0,sin(angles[4]),cos(angles[4])]])

#MrotSa=matrix([[1,0,0],[0,cos(angles[5]),-sin(angles[5])],[0,sin(angles[5]),cos(angles[5])]])

#MrotUr=matrix([[1,0,0],[0,cos(angles[6]),-sin(angles[6])],[0,sin(angles[6]),cos(angles[6])]])

#MrotNe=matrix([[1,0,0],[0,cos(angles[7]),-sin(angles[7])],[0,sin(angles[7]),cos(angles[7])]])


while t_run<t_final:
    sleep(sleeptime)
    Me = MrotMe*transpose(matrix([xxMe[time_i],yyMe[time_i], zp]))
    Ve = MrotVe*transpose(matrix([xxVe[time_i],yyVe[time_i], zp]))
    Te = MrotTe*transpose(matrix([xxTe[time_i],yyTe[time_i], zp]))
    Ma = MrotMa*transpose(matrix([xxMa[time_i],yyMa[time_i], zp]))
    Ju = MrotJu*transpose(matrix([xxJu[time_i],yyJu[time_i], zp]))
    #Sa = MrotSa*transpose(matrix([xxJu[time_i],yyJu[time_i], zp]))
    #Ur = MrotUr*transpose(matrix([xxJu[time_i],yyJu[time_i], zp]))
    #Ne = MrotNe*transpose(matrix([xxJu[time_i],yyJu[time_i], zp]))


    prtclMe.pos = vector( Me[0], Me[1],Me[2])
    prtclVe.pos = vector( Ve[0], Ve[1],Ve[2])
    prtclTe.pos = vector( Te[0], Te[1],Te[2])
    prtclMa.pos = vector( Ma[0], Ma[1],Ma[2])
    prtclJu.pos = vector( Ju[0],Ju[1],Ju[2])
    #prtclSa.pos = vector( Sa[0],Sa[1],Sa[2]) 
    #prtclUr.pos = vector( Ur[0],Ur[1],Ur[2])
    #prtclNe.pos = vector( Ne[0],Ne[1],Ne[2])

    prtclSa.pos = vector( xxSa[time_i],yySa[time_i], zp )
    prtclUr.pos = vector( xxUr[time_i],yyUr[time_i], zp )
    prtclNe.pos = vector( xxNe[time_i],yyNe[time_i], zp )
    t_run += t_delta
    time_i += 1
