import numpy as np
import matplotlib.pyplot as plt
import math

plt.close('all')

# Characteristics of the boeing 747 wing. 

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
mu = 1.8 * 10**(-5)

# Derived parameters
crusie_alt_speed_of_sound = math.sqrt(1.4*287*218.8)
crusie_speed = cruise_mach*crusie_alt_speed_of_sound

# Critical reynolds number for skin drag
Re_crit = 3 * 10**5


def find_chord(section_points):
    chord_per_point =[]
    for i in range(len(section_points)):
        chord_per_point.append(chord_root - (chord_root-chord_tip)*section_points[i]/(wingspan/2))
    return chord_per_point

def find_total_drag(n):
    section_points = []
    for i in range(1,n+1):
        section_points.append(wingspan/(2*i))
    chord_per_point = find_chord(section_points)
    l_crit = Re_crit*mu/(density*crusie_speed)
    #print(section_points)
    #print(chord_per_point)
    total_drag = 0
    for i in range(len(chord_per_point)):
        chord = chord_per_point[i]
        re = density*crusie_speed*chord/mu
        turbulent_drag_whole = .036*density*(crusie_speed**2)*chord*section_points[i]/(re**(1/5))
        turbulent_laminar = .036*density*(crusie_speed**2)*l_crit*section_points[i]/(Re_crit**(1/5))
        laminar_drag = .664*density*(crusie_speed**2)*l_crit*Re_crit**(-1/2)
        total_drag = total_drag + 2*(turbulent_drag_whole-turbulent_laminar+laminar_drag)
    cd = total_drag/(.5*density*(crusie_speed**2)*(wingspan/2)*chord_root)
    return total_drag, cd

drag,cd = find_total_drag(2000)

print("Total viscous drag = ", drag)
print("Drag coefficient derived from viscous drag = ", cd)