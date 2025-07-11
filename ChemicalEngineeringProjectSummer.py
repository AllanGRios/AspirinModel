import math
import matplotlib.pyplot as plt
import numpy as np
#Engineering Project Outline - Chemical Reactor Simulator
pi = 3.14159265

def Maxwell_Boltzman_distribution(mass,temperature,velocity): #kg, K, m/s, J
    # Boltzman Constant
    kb = 1.380649e-23
    
    #Most probable velocity
    v_p = ((2*kb*temperature)/(mass))**(0.5)
    
    #Velocity Distribution
    f_v =4*pi*((mass/(2*pi*kb*temperature))**3/2)*(velocity**2)*math.exp(-(mass*velocity**2)/(2*kb*temperature))
    distribution = []
    velocities = np.linspace(1,5000)
    for v in velocities:
        f = 4 * pi * ((mass / (2 * pi * kb * temperature)) ** 3 / 2) * (v ** 2) * math.exp(
            -(mass * v ** 2) / (2 * kb * temperature))
        distribution.append(f)
        
    #Distribution at most probable velocity
    f_vp = 4*pi*((mass/(2*pi*kb*temperature))**3/2)*(v_p**2)*math.exp(-(mass*v_p**2)/(2*kb*temperature))

    #Graphing
    plt.plot(velocities, distribution)
    plt.vlines(v_p,ymax=f_vp, ymin=0, colors="red")
    plt.title("Maxwell Boltzman Distribution")
    plt.xlabel("Velocity m/s")
    plt.ylabel("probability")
    plt.show()
    return f_v

def salycilic_acid(mass=1,temperature=25):
    Molecular_weight = 0.13812 #kg/mol
    mol = (0.13812/mass) #mol
    return mol, Melting_point, Molecular_weight

def acetic_anhydride(mass=1,temperature=25):
    Molecular_weight = 0.10209 #kg/mol
    mol = (0.10209/mass) #mol
    return mol, Melting_point, Molecular_weight

def CSTR(mass_inflow=1, accumulation = 0, mass_outflow=1, generation=0):
    # dm/dt = Σφ_in - Σφ_out + rV
    return



