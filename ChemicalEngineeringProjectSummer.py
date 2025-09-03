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

def salycilic_acid(V_SA):
    density = 1.44 #kg/L
    mass = V_SA * density #kg/hr
    Molecular_weight = 0.13812 #kg/mol
    cost = 4.73 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return molar_flow, cost, mass

def acetic_anhydride(V_SA, T):
    density = 1.08 #kg/L
    volumetric_flow = 1.25 * V_SA
    mass = volumetric_flow * density # kg/hr
    Molecular_weight = 0.10209 #kg/mol
    cost = 0.76 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr

    if T >= 400:
        molar_flow = molar_flow / (1 + (math.exp(T - 412)))
    return molar_flow, cost, mass

def DMSO(V_SA): # L/hr
    density = 1.1 #kg/L
    volumetric_flow = (-2.25 * V_SA) + 5000
    mass = volumetric_flow * density #kg/hr
    Molecular_weight = 0.07813 #kg/mol
    cost = 4.67 * mass #€/hr
    molar_flow = mass/Molecular_weight #mol/hr
    return molar_flow, cost, mass

def aspirin(Molar_flow):
    density = 1.4 #kg/L
    molecular_weight = 0.180158 #kg/mol
    mass = Molar_flow * molecular_weight
    revenue = 90.90 * mass
    return revenue, mass

def CSTR(v_SA, V_CSTR, T=295):#kg/hr, K
    # Assuming: solvent is always saturated with reactants and V_CSTR Remains Constant
    # in terms of limiting reagent, salicylic acid
    """Constants"""
    V_CSTR, T, v_SA = (float(V_CSTR),float(T), float(v_SA))
    A, R, Ea = [92245, 8.31, 40000]
    cp = 0.7025  # kJ/kg-K, specific heat capacity of mixture
    power_coeff = 0.0032  # €/J

    def reactor_balance(temp):
        # Feed streams and relationships between them
        n_AA, cost_AA, mass_AA = acetic_anhydride(v_SA, temp)
        n_SA, cost_SA, mass_SA = salycilic_acid(v_SA)
        n_DMSO, cost_DMSO, mass_DMSO = DMSO(v_SA)

        # Reactor concentrations [mol/L]
        C_AA = n_AA / V_CSTR
        C_SA = n_SA / V_CSTR

        # Reaction rate [mol/(L·hr)]
        r_AS = A * math.exp(-Ea / (R * temp)) * C_AA * C_SA

        # Aspirin production [mol/hr]
        n_AS = r_AS * V_CSTR
        revenue, mass_AS = aspirin(n_AS)

        # Power cost
        mass_flow = mass_AA + mass_DMSO + mass_SA + mass_AS
        pc = mass_flow * cp * temp * power_coeff

        # Net cash flow
        materials = revenue - (cost_AA + cost_DMSO + cost_SA)
        cash = materials - pc

        return n_AS, n_AA, n_DMSO, cash

    # Single-point results at input T
    n_AS, n_AA, n_DMSO, _ = reactor_balance(T)

    # Sweep over range of temperatures for profit curve
    Cash = {temp: reactor_balance(temp)[3] for temp in range(295, 500)}

    return n_AS, n_AA, n_DMSO, Cash

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
        axs[0][0].annotate(f"T = {temp:.0f} K\nvₚ = {v_p:.1f} m/s\nN% > Φ = {PY:.1f} %", xy=(0.7*max(velocities), max(distribution)*0.8), fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray"))
        axs[0][0].set_title("Maxwell-Boltzman Distribution",fontsize=12, weight="bold")
        axs[0][0].grid(True, alpha=0.3)
        axs[0][0].set_ylim(0,7e-11)
        canvas.draw()

    def update_volume(reactor_volume):
        temp = slider_1.get()
        AS, n_AA, n_DMSO, Cash = CSTR(reactor_volume, 5000, temp)
        Time = np.linspace(0,1000)
        Outlet = []
        for i in Time:
            Outlet.append(AS * i)
        # Aspirin mol vs Time
        axs[0][1].clear()
        axs[0][1].plot(Time, Outlet, color="tab:blue", linewidth=2)
        axs[0][1].set_ylim(0, 1e6)
        axs[0][1].set_title(f"Aspirin Production [{round(AS)} mol/hr]", fontsize=12, weight="bold")
        axs[0][1].set_ylabel("Aspirin (mol)", fontsize=10)
        axs[0][1].grid(True, alpha=0.3)

        # Profit vs temperature
        axs[1][1].clear()
        axs[1][1].plot(list(Cash.keys()), list(Cash.values()), color="tab:green", linewidth=2)
        axs[1][1].axvline(temp, color="red", linestyle="--", linewidth=1.5, label=f"T = {temp:.0f} K")
        axs[1][1].set_title("Profit vs Temperature", fontsize=14, weight="bold")
        axs[1][1].set_xlabel("Temperature (K)", fontsize=12)
        axs[1][1].set_ylabel("Profit (EUR/hr)", fontsize=12)
        axs[1][1].legend()
        axs[1][1].grid(True, alpha=0.3)
        canvas.draw()

    slider_1 = tk.Scale(root, from_=295, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp) #Reactor Temperature Slider
    slider_1.set(298)
    slider_1.pack(pady=5)


    slider_2 = tk.Scale(root, from_=1, to= 500, orient=tk.HORIZONTAL, label="Inlet Stream SA (L/hr)", command=update_volume) #Mass Inflow Slider
    slider_2.set(50)
    slider_2.pack(pady=5)

    btn = tk.Button(root, text="Open New Window", command=open_new_window)
    btn.pack(pady=5)
    root.mainloop()

Window()