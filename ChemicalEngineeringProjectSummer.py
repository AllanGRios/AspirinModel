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
    cost = 4.73 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def acetic_anhydride(volumetric_flow):
    density = 1.08 #kg/L
    mass = volumetric_flow * density # kg/hr
    Molecular_weight = 0.10209 #kg/mol
    cost = 0.76 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def DMSO(volumetric_flow): # L/hr
    density = 1.1 #kg/L
    mass = volumetric_flow * density #kg/hr
    Molecular_weight = 0.07813 #kg/mol
    cost = 4.67 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return cost, mass, molar_flow

def CSTR(v_SA,  V_CSTR, T=295): #kg/hr, K
    # GOAL: Find molar flow of Aspirin at the outlet, find the conversion, plot a graph of Aspirin Production vs Time at that specific temperature and that specific volume
    # Assuming: solvent is always saturated with reactants and V_CSTR Remains Constant
    # in terms of limiting reagent, salicylic acid
    V_CSTR, T, v_SA = (float(V_CSTR),float(T), float(v_SA))
    A, R, Ea = [1e+8, 8.31, 40000]
    # Relationship between volumetric flows such that solvent is always saturated and Volume always constant
    v_DMSO = -2.25 * v_SA + 1000
    v_AA = 1.25 * v_SA
    # Reactor Concentrations
    C_AA = ((1.08 * v_AA)/(0.10209 * (v_DMSO + v_AA + v_SA)))
    C_SA = ((1.44 * v_SA)/(0.13812 * (v_DMSO + v_AA + v_SA)))
    # Consumption
    r_SA = (A * math.exp((-Ea)/(R * T))) * C_AA * C_SA
    # Conversion & Aspirin Flow
    X_SA = 100 * ((-r_SA) * V_CSTR * 0.13812)/(1.44 * v_SA)
    n_AS = (r_SA * V_CSTR)

    def net_cash_flow():
        # Q = mc∆T   Kj/hr = Kwh
        mass = (1.1 * v_DMSO) + (1.08 * v_AA) + (1.44 * v_SA)
        Power = mass * 412 * T
        # Operational Costs
        Total_Costs = -((0.76 * 1.08 * v_AA) + (4.73 * 1.44 * v_SA) + (4.67 * 1.1 * v_DMSO) )
        Total_Profit = (90.90 * n_AS) + Total_Costs
        return Total_Profit

    return n_AS, net_cash_flow(), v_AA, v_DMSO # NEED DISTRIBUTION TO FORM GRAPH

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
        volume = slider_2.get()
        update_volume(volume)
        # Maxwell-Boltzman Graph
        m = 2.29e-25
        PY, N = Energy_Distribution(mol, T=temp)
        distribution, velocities, v_p, f_vp = Maxwell_Boltzman_distribution(m, temp)
        axs[0][0].clear()
        axs[0][0].plot(velocities, distribution)
        axs[0][0].set_title(f"T = {temp:.0f} K \n Most Probable Velocity = {v_p:.1f} ms⁻¹ \n N% Above Φ = {PY:.1f} %")
        canvas.draw()

    def update_volume(reactor_volume):
        temp = slider_1.get()
        AS, profit, v_AA, v_DMSO = CSTR(reactor_volume, 1000, temp)
        Time = np.linspace(0,1000)
        Outlet = []
        for i in Time:
            Outlet.append(AS * i)
        # Aspirin mol vs Time
        axs[0][1].clear()
        axs[0][1].set_ylim(0,10000000)
        axs[0][1].plot(Time, Outlet)
        axs[0][1].set_title(f"Production of aspirin in mol/hr")
        axs[0][1].set_ylabel(f"mol")
        canvas.draw()

    slider_1 = tk.Scale(root, from_=200, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp) #Reactor Temperature Slider
    slider_1.set(273)
    slider_1.pack(pady=5)


    slider_2 = tk.Scale(root, from_=1, to= 100, orient=tk.HORIZONTAL, label="Inlet Stream SA (L/hr)", command=update_volume) #Mass Inflow Slider
    slider_2.set(50)
    slider_2.pack(pady=5)

    btn = tk.Button(root, text="Open New Window", command=open_new_window)
    btn.pack(pady=5)
    root.mainloop()

Window()
