from scipy.integrate import simps
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#from Pendulum_Animations import animate_pendulum

##INPUTS######################################################################################

plt.figure()

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
c = 240.3 #Damping coefficient of damper, Ns/m

##LEAD SCREW PARAMETERS
l = 0.005 #lead length, m
v = 15*np.pi/180 #Thread angle, radians
D = 0.024 #Thread diameter,m
C = 0.14 #Coefficient of friction between thread and nut
Cl=np.pi*D #Circumfernce of the lead screw
R=l/Cl #Ratio of the lead to the circumfernce of the lead screw
h=np.arctan(R)
Ef = (np.tan(h))/ np.tan(h + np.arctan(C)) #Lead screw efficiency


deployment_times = []
max_dthet_list = []
Energy = []
j=1

T = 400  # Defined outside the loop
rtol = 1e-6
z0 = np.array([Q1 * np.pi / 180, 0])  # Defined outside the loop

#WIND PARAMETERS
V = 10 #max wind speed, m/s
p = 1.293 #density, kg/m^3

motor_configurations = [
    {'Ts': 0.19, 'Wnl': 5000*2*np.pi/60},  # Motor configuration 1
    {'Ts': 0.43, 'Wnl': 2700*2*np.pi/60},  # Motor configuration 2
    {'Ts': 0.48, 'Wnl': 8000*2*np.pi/60},   # Motor configuration 3
    {'Ts': 0.69, 'Wnl': 3700*2*np.pi/60}    # Motor configuration 4
]

# Loop over each motor configuration and its corresponding gear ratio range
for motor_config in motor_configurations:
    Ts = motor_config['Ts']
    Wnl = motor_config['Wnl']
    deployment_times = []
    energy_values =[]
    if j==1:
        gear_ratios = np.arange(55, 63, 0.5)
    if j==2:
        gear_ratios = np.arange(24, 32, 0.5)
    if j==3:
        gear_ratios = np.arange(22, 30, 0.5)
    if j==4:
        gear_ratios = np.arange(15, 23, 0.5)
    for gear_ratio in gear_ratios:
    ##MOTOR AND GEAR BOX PARAMETERS
        Gr = gear_ratio #Gear ratio
        Tg = Ts*Gr #Gear torque, Nm 
        SF = C*np.pi*D*(1/np.cos(v))-l #Self locking in lead screw
        F= 2*Tg /(D*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v)))) #Force input from motor, N
        Tf = (F*l)/(2*np.pi*Ef) #Actual torque on lead screw
        Fn = 2*Tf /(D*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v)))) #Actual force from motor
    
    ######################################################################################################
    ##EQUATION OF MOTION
        crossed_Q2 = False
        
        def f(t, z):
            global crossed_Q2
            global Ta
       
            if z[0] <= Q2 * np.pi / 180 and not crossed_Q2:
                x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(z[0]))
                dx = (z[1])*(-d*a/x)*np.sin(z[0])
                Fd = 2*c*dx
           
                Ta = a*(Fn*((np.sin(z[0]) + (np.sin(2*z[0]))/(2*n*np.sqrt((n**2)-(np.sin(z[0]))**2))) + Fd*(d/x)*np.sin(z[0]))) #Accelerating torque
                Tr = 8*m*g*L*(np.cos(z[0])) #Resistive torque
                Tw = 4*w*p*(L*V*np.sin(z[0]))**2 #Wind load
                dz = [z[1]*R, (Ta - Tr - Tw)/I]
           
            # Reset the flag if z[0] goes below Q2*np.pi/180 again
                crossed_Q2 = False
           
            if z[0] > Q2 * np.pi / 180 or crossed_Q2:
                x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(z[0]))
                dx =(z[1])*(-d*a/x)*np.sin(z[0])
                Fd = 2*c*dx
           
                Ta = a*(-Fn*((np.sin(z[0]) + (np.sin(2*z[0]))/(2*n*np.sqrt((n**2)-(np.sin(z[0]))**2))) + Fd*(d/x)*np.sin(z[0])))
                Tr = 8*m*g*L*(np.cos(z[0])) #Resistive torque
                Tw = 4*w*p*(L*w*np.sin(z[0]))**2 #Wind load
                dz = [-z[1]*R, -(Ta - Tr - Tw)/I]
           
            # Set the flag to True once z[0] crosses Q2*np.pi/180
                crossed_Q2 = True
           
            return dz
            
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
        deployment_times.append(t.max())
        max_dthet_list.append(dthet.max())
        #print('{:.2f}'.format(t.max()) + 's') #Time taken to deploy
    #CALCULATING HOLDING TORQUE
        x = np.sqrt(d**2 + a**2 - 2*d*a*np.cos(thet))
        dx = (dthet1)*(-d*a/x)*np.sin(thet)
        Fd = 2*c*dx
        Tw = 4*w*p*(L*V*np.sin(thet))**2 #Wind load
        Tr = 8*m*g*L*(np.cos(thet)) #Resistive torque
        Th = (1/Gr)*(2*np.pi*Ef/l) * (((D/2)*((l + C*np.pi*D*(1/np.cos(v)))) / (np.pi*D - C*l*(1/np.cos(v))))**2) * ((((Tw+Tr)/a - Fd*(d/x)*np.sin(thet)))/ ((np.sin(thet) + (np.sin(2*thet))/(2*n*np.sqrt((n**2)-(np.sin(thet))**2)))))  #Holding torque as angle changes
        
        W = (-Wnl/Ts)*(Th-Ts) #Rotation rate of motor, rad/s
        P = Th*W #Power output of motor, W
        A = simps(P,t)
        energy_values.append(A)
        print(A)
        
    Energy.append(energy_values)
    j=j+1

    # Plot energy values for each motor configuration
    plt.plot(gear_ratios,deployment_times, '.-', label=f'Motor {j-1}')

# Add labels and legend to the plot
plt.xlabel('Gear Ratio')
plt.ylabel('Deployment Time (s)')

plt.legend()
plt.grid(True)
plt.show()
