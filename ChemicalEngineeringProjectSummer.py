import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.integrate import quad
import tkinter as tk
from tkinter import ttk
#Engineering Project Outline - Chemical Reactor Simulator
pi = 3.14159265
kb = 1.380649e-23
Av = 6.02214e+23

def START(Mol, Temp, Time, V_threshold):
    Mol = float(Mol)
    Temp = float(Temp)
    Time = float(Time)
    V_Threshold = float(V_threshold)
    return Mol, Temp, Time, V_threshold

#Classical model based on translational energy
def Maxwell_Boltzman_distribution(mass,temperature, v_max=1000): #kg, K, m/s, J
    #Most probable velocity
    v_p = ((2*kb*temperature)/(mass))**(0.5)
    
    #Velocity Distribution
    #f_v =4*pi*((mass/(2*pi*kb*temperature))**3/2)*(velocity**2)*math.exp(-(mass*velocity**2)/(2*kb*temperature))
    distribution = []
    velocities = np.linspace(0, v_max, v_max)
    for v in velocities:
        f = 4 * pi * ((mass / (2 * pi * kb * temperature)) ** 3 / 2) * (v ** 2) * math.exp(
            -(mass * v ** 2) / (2 * kb * temperature))
        distribution.append(f)

    #Distribution at most probable velocity
    f_vp = 4*pi*((mass/(2*pi*kb*temperature))**3/2)*(v_p**2)*math.exp(-(mass*v_p**2)/(2*kb*temperature))
    return distribution,velocities ,v_p, f_vp

def Energy_Distribution(mol=1, mass=2.29e-25, T=295): #salicylic acid kg/mol
    V_min = 500 # Activation energy used as kinetic to find velocity Ek = 1/2 mv^2
    C = mass/(2*kb*T)
    N_total = mol*Av
    def f(v):
        dist = (4*pi*v**2)*((C/pi)**(1.5))*math.exp(-C*v**2)
        return dist

    integral, _ = quad(f, V_min, 6000)
    N = N_total * integral
    Percentage_yield = ( N/N_total ) * 100
    return Percentage_yield, N

def salycilic_acid(volumetric_flow):
    density = 1.44 #kg/L
    mass = volumetric_flow * density #kg/hr
    Molecular_weight = 0.13812 #kg/mol
    cost = 80.55 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def acetic_anhydride(volumetric_flow):
    density = 1.08 #kg/L
    mass = volumetric_flow * density # kg/hr
    Molecular_weight = 0.10209 #kg/mol
    cost = 43.78 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def DMSO(volumetric_flow): # L/hr
    density = 1.1 #kg/L
    mass = volumetric_flow * density #kg/hr
    Molecular_weight = 0.07813 #kg/mol
    cost = 22.83 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def CSTR(T=295, V=50): #kg/hr, K
    # GOAL: Find molar flow of Aspirin, conversion and Time extinction coefficient
    # Assuming solvent is saturated with reactants
    vol_DMSO = V/1.225
    vol_AA = 0.125 * vol_DMSO
    vol_SA = 0.1 * vol_DMSO
    R = 8.31
    Ea = 40000 # j/mol Assumption based on common esterification values
    A = 1e+11 # Assumption based on similar Liquid Phase Esterificaiton Reactions
    #Step 1: Finding Concentration of inlet streams SA & AA
    DMSO_C,DMSO_M, DMSO_F = DMSO(vol_DMSO)
    AA_C, AA_M, AA_F = acetic_anhydride(vol_AA)
    SA_C, SA_M, SA_F = salycilic_acid(vol_SA)
    C_SA0 = SA_F/vol_DMSO
    C_AA0 = AA_F/vol_DMSO
    #Step 2: Find Rate of consumption
    K_r = A*(math.exp(-Ea/(R*T)))
    r_consumption_SA = K_r * C_SA0 * C_AA0


    return

def open_new_window():
    new_win = tk.Toplevel(root)
    new_win.title("Secondary Graph Window")
    new_fig, new_ax = plt.subplots()
    new_canvas = FigureCanvasTkAgg(new_fig, master=new_win)
    new_canvas.get_tk_widget().pack()

    # Plot example graph
    v = np.linspace(1, 5000, 500)
    new_ax.plot(v, np.sin(v / 500))
    new_ax.set_title("Sine Wave Example")
    new_canvas.draw()

def Window():

    # main window
    root = tk.Tk()
    root.title("System Visual")
    root.geometry("1500x1000")
    frame = tk.Frame(root)
    frame.pack()

    # 4 sub plots
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.get_tk_widget().pack()

    # slider
    def update_temp(temp, mol=0.03):
        mol = float(mol)
        temp = float(temp)
        # Maxwell-Boltzman Graph
        m = 2.29e-25
        PY, N = Energy_Distribution(mol, T=temp)
        distribution, velocities, v_p, f_vp = Maxwell_Boltzman_distribution(m, temp)
        axs[0][0].clear()
        axs[0][0].plot(velocities, distribution)
        axs[0][0].set_title(f"T = {temp:.0f} K \n Most Probable Velocity = {v_p:.1f} ms⁻¹ \n N% Above Φ = {PY:.1f} %")
        canvas.draw()

    def update_cost(mass_inflow):
        temp = slider_1.get()

        return


    slider_1 = tk.Scale(root, from_=200, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp) #Reactor Temperature Slider
    slider_1.set(273)
    slider_1.pack(pady=20)

    slider_2 = tk.Scale(root, from_=0, to= 100, orient=tk.HORIZONTAL, label="Mass Inflow Salicylic Acid (kg/hr)", command=update_cost) #Mass Inflow Slider
    slider_2.set(5)
    slider_2.pack(pady=50)

    btn = tk.Button(root, text="Open New Window", command=open_new_window)
    btn.pack(pady=20)
    root.mainloop()

x=100
CSTR(V=x)