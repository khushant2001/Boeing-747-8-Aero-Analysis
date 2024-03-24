
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sp
import math

plt.close('all')

# Given parameters. 
chord = 14.84 #m
wingspan = 68.4 #m
max_weight = 447700*9.81 #N
cruise_mach = .84
wing_area = 554 # m^2
density = 1.225 # kg/m^3

# Derived parameters
crusie_alt_speed_of_sound = math.sqrt(1.4*287*218.8)
crusie_speed = cruise_mach*crusie_alt_speed_of_sound

# Getting points
bacxxx = pd.read_excel("bacxxx cordinates.xlsx", header = 0)
x_cord = chord*bacxxx.iloc[:,0]
y_cord = chord*bacxxx.iloc[:,1]


mean_camber_line = []
for i in range(int(len(y_cord)/2)):
    mean_camber_line.append((y_cord[i] + y_cord[i+len(y_cord)/2])/2)

plt.figure()
plt.scatter(x_cord,y_cord, label = 'Airfoil geometry')
plt.scatter(x_cord[0:int(len(x_cord)/2)], mean_camber_line, label = 'Mean Camber Line')
plt.title("BACXXX Airfoil")
plt.legend()
plt.xlabel("x/c")
plt.ylabel("z/c")
plt.grid()

# Parameters of the airfoil
m = 1.4/100
p =  15/100
t = 11.3

def thin_arifoil(m,p, alpha):
    # Calculating coefficients!
    phi = sp.symbols('phi')
    x = sp.symbols('x')
    mean_camber_func1 = (m/p**2)*(2*p*x - x**2)
    mean_camber_func2 = (m/(1-p)**2)*((1-2*p) + 2*p*x - x**2)

    dz_dx_1_1 = sp.diff(mean_camber_func1, x)
    dz_dx_2_1 = sp.diff(mean_camber_func2, x)

    dz_dx_1 = dz_dx_1_1.subs(x,.5*(1-sp.cos(phi)))
    dz_dx_2 = dz_dx_2_1.subs(x,.5*(1-sp.cos(phi)))

    break_point = np.arccos(1 - p*.5)
    a0_minus_alpha_1 = (1/np.pi)*sp.integrate(dz_dx_1, (phi, 0, break_point))
    a0_minus_alpha_2 = (1/np.pi)*sp.integrate(dz_dx_2, (phi, break_point, np.pi))
    a0_minus_alpha = a0_minus_alpha_1 + a0_minus_alpha_2

    a1_1 = (2/np.pi)*sp.integrate(dz_dx_1*sp.cos(phi), (phi, 0, break_point))
    a1_2 = (2/np.pi)*sp.integrate(dz_dx_2*sp.cos(phi), (phi, break_point, np.pi))
    a1 = a1_1 + a1_2

    a2_1 = (2/np.pi)*sp.integrate(dz_dx_1*sp.cos(2*phi), (phi, 0, break_point))
    a2_2 = (2/np.pi)*sp.integrate(dz_dx_2*sp.cos(2*phi), (phi, break_point, np.pi))
    a2 = a2_1 + a2_2

    cl = 2*np.pi*(alpha - a0_minus_alpha +a1/2)
    cm_c_4 = -cl*.25 + np.pi*.25*(a1-a2)

    return cl, cm_c_4

alpha = np.linspace(math.radians(-10),math.radians(10),100)
alpha_deg = np.linspace(-10,10,100)
cl,cm_c_4 = thin_arifoil(m,p,alpha)

print((cl[-1]-cl[0])/(alpha[-1]-alpha[0]))
plt.figure()
plt.plot(alpha_deg,cl)
idx = np.argwhere(np.diff(np.sign(cl - 0))).flatten()
plt.plot(alpha_deg[idx], cl[idx], 'ro')
plt.grid()
plt.xlabel("Angle of attack [deg].")
plt.title(" Cl vs AOA. Aplha_0 = " + str(alpha_deg[idx]) + "deg")
plt.ylabel("C_L")

plt.figure()
plt.plot(alpha_deg,cm_c_4)
plt.grid()
plt.xlabel("Angle of attack [deg]")
plt.ylabel("C_M_Leading Edge")

lift_force = .5*cl*density*(crusie_speed**2)*wing_area

plt.figure()
plt.plot(alpha_deg,lift_force, label = 'Lift force', color = 'green')
plt.axhline(y = max_weight, xmin = -10,xmax=10, label = 'Max weight', linestyle = '--')
plt.grid()
plt.xlabel("Angle of attack [deg]")
plt.legend()
idx = np.argwhere(np.diff(np.sign(lift_force - max_weight))).flatten()
plt.plot(alpha_deg[idx], lift_force[idx], 'ro')
plt.title("Lift force vs angle of degrees. AOA at cruise = " + str(alpha_deg[idx]) + "deg")
plt.show()