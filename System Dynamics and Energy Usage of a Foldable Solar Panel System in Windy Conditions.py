# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 13:57:49 2024

@author: kaan2
"""
from scipy.integrate import simps
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#from Pendulum_Animations import animate_pendulum

##INPUTS######################################################################################

##MECHANISM PARAMETERS
m = 8.5 #Mass of solar panel, kg
g = 9.81 #gravitational acceleration, m/s^2
L = 0.55 #Width of panel, m
w = 1.2 #Length of panel, m
I = ((4*m)*(2*L)**2) #Moment of inertia, m^4
Q1 = 5.8 #Initial angle, degrees theta
Q2 = 48.6 #Final angle first phase, degrees psi
Q3 = 33 #Final angle, degrees
a = 0.4 #Length to joint, m , fixed
b = 0.305 #Length to joint, m , moving
d = 0.125 #Damper distance from pivot point, m
n = b/a

#DAMPER DA500ADJ
c = 240.3 #Damping coefficient of damper, Ns/m

##LEAD SCREW PARAMETERS
l = 0.005 #lead length, m
v = 15*np.pi/180 #Thread angle, radians
D = 0.024 #Thread diameter,m
C = 0.14 #Coefficient of friction between thread and nut
Cl=np.pi*D #Circumfernce of the lead screw
R = l/Cl #Ratio of the lead to the circumfernce of the lead screw
h = np.arctan(R)
Ef = (np.tan(h))/ np.tan(h + np.arctan(C)) #Lead screw efficiency

#WIND PARAMETERS
V = 10 #max wind speed, m/s
p = 1.293 #density, kg/m^3

##MOTOR AND GEAR BOX PARAMETERS
#Part no. 0 390 203 389
Ts = 0.43 #Stall Torque, Nm
Gr = 30.58 #Gear ratio
Tg = Ts*Gr #Gear torque, Nm
Wnl = 2700*2*np.pi/60 #No load speed, rad/s
SF = C*np.pi*D*(1/np.cos(v))-l #Self locking in lead screw

VR = 1/R #velocity ratio

F= 2*Tg /(D*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v)))) #Force input from motor, N
Tf = (F*l)/(2*np.pi*Ef) #Actual torque on lead screw
Fn = 2*Tf /(D*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v)))) #Actual force from motor

########################EQUATION OF MOTION#######################################
crossed_Q2 = False

def f(t, z):
    global crossed_Q2 
   
    if z[0] <= Q2 * np.pi / 180 and not crossed_Q2:
        I = ((4*m)*(a)**2) #Moment of inertia, m^4
        x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(z[0]))
        dx = (z[1])*(d*a/x)*np.sin(z[0])
        Fd = 2*c*dx
       
        Ta = a*(Fn*((np.sin(z[0]) + (np.sin(2*z[0]))/(2*np.sqrt((n**2)-(np.sin(z[0]))**2))) - Fd*(d/x)*np.sin(z[0]))) #Accelerating torque
        Tr = 8*m*g*L*(np.cos(z[0])) #Resistive torque
        Tw = 4*w*p*(L*V*np.sin(z[0]))**2 #Wind load
        dz = [z[1]/VR, (Ta - Tr - Tw)/I]
       
        # Reset the flag if z[0] goes below Q2*np.pi/180 again
        crossed_Q2 = False
       
    if z[0] > Q2 * np.pi / 180 or crossed_Q2:
        I = ((4*m)*(2*L)**2) #Moment of inertia, m^4
        x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(z[0]))
        dx =(z[1])*(d*a/x)*np.sin(z[0])
        Fd = 2*c*dx
       
        Ta = a*(-Fn*((np.sin(z[0]) + (np.sin(2*z[0]))/(2*np.sqrt((n**2)-(np.sin(z[0]))**2))) - Fd*(d/x)*np.sin(z[0])))
        Tr = 8*m*g*L*(np.cos(z[0])) #Resistive torque
        Tw = 4*w*p*(L*V*np.sin(z[0]))**2 #Wind load
        dz = [-z[1]/VR, -(Ta - Tr - Tw)/I]
       
        # Set the flag to True once z[0] crosses Q2*np.pi/180
        crossed_Q2 = True
       
    return dz

# Solver paramaters
T = 400
rtol = 1e-7
z0 = [Q1*np.pi/180,0]

##EVENT FUNCTION
# Variable to track the last angle reached
last_angle = Q1 * np.pi / 180

# Variable to count the number of times the event has occurred
def event1(t, z):
    # Check if the peak Q2 has been passed and z[0] is equal to or less than Q3
    if crossed_Q2 and z[0] <= Q3 * np.pi / 180:
        return 0  # Return 0 to indicate that the event condition is met
    else:
        return 1  # Return 1 to continue integration

# Set terminal property to True
event1.terminal = True

##IVP solver
sol = solve_ivp(f,(0,T),z0, rtol=rtol, events=[event1])

##OUTPUTS############################################################################################
t = sol.t
z = sol.y
thet=z[0,:]
dthet=z[1,:]
dthet1=z[1,:]/Gr

#CALCULATING HOLDING TORQUE
x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(thet))
dx = (dthet)*(d*a/x)*np.sin(thet)
Fd = 2*c*dx
Tw = 4*w*p*(L*V*np.sin(thet))**2 #Wind load
Tr = 8*m*g*L*(np.cos(thet)) #Resistive torque
T = (1/Gr)*(2*np.pi*Ef/l) * (((D/2)*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v))))**2) * ((((Tw+Tr)/a + Fd*(d/x)*np.sin(thet)))/ ((np.sin(thet) + (np.sin(2*thet))/(2*n*np.sqrt((n**2)-(np.sin(thet))**2)))))  #Holding torque as angle changes

W = (-Wnl/Ts)*(T-Ts) #Rotation rate of motor, rad/s
P = T*W #Power output of motor, W
i = (T/Ts) #Current input to motor, A

##VERTICAL LOAD ON WHEELS
# Fv = Fn*(np.sin(thet)/2*n*np.sqrt((n**2)-(np.sin(thet))**2)) #Vertical load on wheel
Tw = 4*w*p*(L*V*np.sin(thet))**2 #Wind load
Tr = 8*m*g*L*(np.cos(thet)) #Resistive torque due to mass
Fv = (np.sin(thet)/n) * (np.sin(thet)/(n*np.sqrt((n**2)-(np.sin(thet))**2))) * ((Tw+Tr)/a - Fd*(d/x)*np.sin(thet))

Ft=(Fn*n)/(np.sqrt((n**2)-(np.sin(thet))**2)) * (np.sin(thet) + np.sin(2*thet)/(2*n*np.sqrt((n**2)-(np.sin(thet))**2))) - (Tw+Tr)/a - Fd*(d/x)*np.sin(thet)

##ANGLE PLOT OVER TIME
fig, ax1 = plt.subplots()
ax1.plot(t,thet*180/np.pi,'.-',color='red')
ax1.set_ylim(bottom=0, top=(50))
ax1.set_xlabel('Time,t')
ax1.set_ylabel('Angle,θ (degrees)')

##VELOCITY PLOT OVER TIME
fig2, ax2 = plt.subplots()
ax2.plot(t,dthet1,'.-',color='blue')
ax2.set_xlabel('Time,t')
ax2.set_ylabel('Velocity,w (rad/s)')

#POWER OUTPUT AND CURRENT OVER TIME
fig2, ax3 = plt.subplots()
ax3.plot(t,P,'.-',color='green',label='Power')
ax3.plot(t,i,'.-',color='Brown',label='Current')
ax3.set_xlabel('Time,t')
ax3.set_ylabel('Power,P (W) & Current,i (A)')
plt.legend()

#MOTOR TORQUE OVER TIME
fig2, ax4 = plt.subplots()
ax4.plot(T,W*30/np.pi,'.-',color='purple')
ax4.set_ylabel('Motor speed w (rpm)')
ax4.set_xlabel('Torque, T (Nm)')

VO=P/i

fig2, ax5 = plt.subplots()
ax5.plot(t,VO/0.33,'.-',color='turquoise')
ax5.set_xlabel('Time,t')
ax5.set_ylabel('Wheel laod,Fv (N)')


print('Max Velocity = ' + '{:.2f}'.format(dthet.max()) + 'rad/s'' Force input = '+ '{:.1f}'.format(Fn) + 'N')
print('Deployment time ' + '{:.2f}'.format(t.max()) + 's') #Time taken to deploy
print('Load on load sensor ' + '{:.2f}'.format(Fv[-1]/2) + 'N')
print('Max load on wheel ' + '{:.2f}'.format(Fv.max()/2.5) + 'N')
print('Max torque on motor ' + '{:.2f}'.format(T.max()) + 'Nm')
print(W.max()*60/(2*np.pi)/Gr)
#print(P.max())

##BUCKLING LOAD and CRITICAL SPEED of lead screw
fb=2 #Fixed Fixed
E=196*1000 #Youngs modulus of 300 series stainless steel GPa
Dr=24 #lead screw diameter, mm
I=np.pi*Dr**4/64 #Momnet of inertia, mm^4
Lc=643 #max length of lead screw, mm
BL=fb*(np.pi)**2*(E*I)/Lc**2 #Buckling load, N

fc=15.1 #Fixed Fixed
Nc=fc*(Dr/(Lc**2))*(10**7) #Critical speed, rpm

print('Buckling load at max length ' + '{:.1f}'.format(BL/1000) + 'kN') 
print('80% of critical speed ' + '{:.1f}'.format(Nc*0.8) + 'rpm')

# Find the index of the maximum value of dthet*Gr
max_index = np.argmax(W*60/(2*np.pi))
max_index2 = np.argmax(P)

# Print the corresponding value of T
#print("Power at max rpm:", (P[max_index]))
#print("Rpm at max Power:", (W[max_index2])*60/(2*np.pi))

Energy = simps(P,t)

print('Energy usage ' + '{:.1f}'.format(Energy) + 'J')
print('Peak current ' + '{:.2f}'.format(i.max()) + 'A')
print(W.max())
print(SF)
print(Ef)
print(Ft.max()/50)
print(VO.max()/0.85)



