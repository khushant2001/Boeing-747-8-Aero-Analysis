import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import math

plt.close('all')

# Given parameters. 
chord_root = 14.7 #m
chord_tip = 3.63 #m
# Chord parameters derived from https://aviation.stackexchange.com/questions/82188/what-is-the-exact-wing-chord-length-and-thickness-for-boeing-747#:~:text=As%20you%20can%20see%2C%20the,and%20wider%20with%20raked%20wingtips.

wingspan = 68.58 #m
AR = 6.97
sweep = 37 #deg
twist = 3 #deg
taper_ratio = chord_tip/chord_root
max_weight = 447700*9.81 #N
cruise_mach = .84
wing_area = 554 # m^2
density = 1.225 # kg/m^3

# Derived parameters
crusie_alt_speed_of_sound = math.sqrt(1.4*287*218.8)
crusie_speed = cruise_mach*crusie_alt_speed_of_sound

# Things derived from thin airfoil theory!
m0 = 6.28 # 1/rad
alpha_l_0 = -1.31 #deg

# Parameter to vary
n_coefficients = 10

# Creating temporary array to store the odd indexes since the wing is symmetric.
odd = np.arange(1,n_coefficients*2,step=2)
theta = np.linspace(0.01,np.pi/2,n_coefficients)

# Function definitions!
def theta_to_y(theta):
    return .5*wingspan*np.cos(theta)

def aoa_vs_span(alpha,y):
    return math.radians(alpha - twist*y/(wingspan/2))

def chord_vs_span(y):
    return chord_root - (chord_root-chord_tip)*y/(wingspan/2)

def find_row(n,t):
    a_row = []
    for i in range(n):
        a_row.append((4*wingspan/(m0*chord_vs_span(theta_to_y(t))))*np.sin(odd[i]*t) + odd[i]*np.sin(odd[i]*t)/np.sin(t))
    return a_row

def find_coefficients(alpha):
    b = np.zeros((n_coefficients,1))
    a = np.zeros((n_coefficients,n_coefficients))
    for n in range(n_coefficients):
        b[n] = aoa_vs_span(alpha,theta_to_y(theta[n])) - math.radians(alpha_l_0)
        a[n:,] = find_row(n_coefficients,theta[n])
    coefficients = np.linalg.solve(a,b)
    return coefficients

def find_constants(alpha):
    coefficients = find_coefficients(alpha)
    cl = coefficients[0]*np.pi*AR
    delta = 0
    induced_aoa = 0
    for i in range(n_coefficients):
        delta = delta + odd[i]*(coefficients[i]/coefficients[0])**2
        induced_aoa = induced_aoa + odd[i]*coefficients[i]*np.sin(odd[i]*theta[i])/np.sin(theta[i])
    delta = delta-odd[0]
    e = 1/(1+delta)
    cd = cl**2/(np.pi*AR*e)
    return cl,delta,e,cd,cl/cd,induced_aoa

# Using the functions to get data.
alpha_deg = np.linspace(-20,20,100)
cl = []
cd = []
cl_cd = []
for i in range(len(alpha_deg)):
    constants = find_constants(alpha_deg[i])
    #print(constants)
    cl.append(constants[0].item())
    cd.append(constants[3].item())
    cl_cd.append(constants[4].item())

# Plotting

plt.figure()
plt.plot(alpha_deg,cl)
plt.grid()
plt.xlabel("Alpha")
plt.ylabel("C_L")
plt.title("C_L vs Alpha")

plt.figure()
plt.plot(alpha_deg,cd)
plt.grid()
plt.xlabel("Alpha")
plt.ylabel("C_D")
plt.title("C_D vs Alpha")

plt.figure()
plt.plot(cd,cl)
plt.xlabel("C_D")
plt.ylabel("C_L")
plt.title("C_L vs C_D")
plt.grid()

plt.figure()
plt.plot(alpha_deg,cl_cd)
plt.grid()
plt.title("C_L/C_D vs Alpha")
plt.xlabel("Alpha")
plt.ylabel("C_L/C_D")
plt.show()