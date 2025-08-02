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

def CSTR(v_SA,  V_CSTR, T=295): #kg/hr, K
    # GOAL: Find molar flow of Aspirin at the outlet, find the conversion, plot a graph of Aspirin Production vs Time at that specific temperature and that specific volume
    # Assuming: solvent is always saturated with reactants and V_CSTR Remains Constant
    # in terms of limiting reagent, salicylic acid
    V_CSTR, T, v_SA = (float(V_CSTR),float(T), float(v_SA))
    A, R, Ea = [1e-2, 8.31, 40000]
    # Relationship between volumetric flows such that solvent is always saturated and Volume always constant
    """v_DMSO = (((0.775 * V_CSTR) - (V_CSTR))/(0.1 * V_CSTR)) * v_SA + V_CSTR
    v_AA = (((0.125 * V_CSTR) - (V_CSTR))/(0.1 * V_CSTR)) * v_SA"""
    v_DMSO = -2.25 * v_SA + 1000
    v_AA = 1.25 * v_SA
    print(f"v_AA {v_AA}, v_DMSO {v_DMSO}")
    print(f"total volume {v_SA + v_AA + v_DMSO} / {V_CSTR} L")
    # Inlet Concentrations
    C_AA0 = ((1.08 * v_AA)/(v_DMSO * 0.10209))
    C_SA0 = ((1.44 * v_SA)/(v_DMSO* 0.13812))
    print(f"concentrations{C_SA0, C_AA0}")
    # Consumption
    r_SA = (A * math.exp((-Ea)/(R * T))) * C_AA0 * C_SA0
    print(f"rate of consumption {r_SA} [units]")
    # Conversion & Aspirin Flow
    X_SA = ((-r_SA) * V_CSTR * 0.13812)/(1.44 * v_SA)
    n_AS = (r_SA * V_CSTR)
    print(f"aspirin {n_AS} and conv {X_SA}")

    def net_cash_flow():
        # Q = mc∆T   Kj/hr = Kwh
        mass = (1.1 * v_DMSO) + (1.08 * v_AA) + (1.44 * v_SA)
        Power = mass * 412 * T
        # Operational Costs
        Total_Costs = -((43.78 * 1.08 * v_AA) + (80.55 * 1.44 * v_SA) + (22.83 * 1.1 * v_DMSO) + (0.32 * Power))
        Total_Profit = (90.90 * n_AS) + Total_Costs
        print(Total_Profit)
        return Total_Profit

    return n_AS, X_SA, net_cash_flow() # NEED DISTRIBUTION TO FORM GRAPH

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

    def update_volume(reactor_volume):
        temp = slider_1.get()
        AS, Volume, profit = CSTR(reactor_volume, 1000, temp)
        Time = np.linspace(0,1000)
        axs[1][0].clear()
        axs[0][1].plot(Time, Outlet)
        axs[0][1].set_title(f"Production of aspirin in mol/hr")
        canvas.draw()

    slider_1 = tk.Scale(root, from_=200, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp) #Reactor Temperature Slider
    slider_1.set(273)
    slider_1.pack(pady=10)

    slider_2 = tk.Scale(root, from_=1, to= 100, orient=tk.HORIZONTAL, label="Inlet Stream SA (L/hr)", command=update_volume) #Mass Inflow Slider
    slider_2.set(50)
    slider_2.pack(pady=10)

    btn = tk.Button(root, text="Open New Window", command=open_new_window)
    btn.pack(pady=20)
    root.mainloop()

Window()
"""
ABOUT.me
# AspirinModel

***

## Project Description
> AspirinModel is a python software developed by an undergraduate Chemical Engineering student in order to simulate different aspects of a Continuous Stirr Tank Reactor (CSTR), patching together my understanding of different fundamental concepts. The project aims at creating a model that approximates the real world physics of running a CSTR reactor. The project will be looking at the industrial process at synthesizing Aspirin starting from Acetic Anhydride and Salicylic Acid under an acid catalyst.

## Purpose
> The purpose of the project is the development of a comprehensive software model of a CSTR that can predict the production of Aspirin and adapt to changes in reaction conditions, giving results in the form of graphs, helping out reactor designers in choosing optimal conditions for the industrial level reaction.

## Scope
> The model must be able to visualize data of physical quantities on graphs including Maxwell-Boltzman Distribution, Operating Cost, Molar Flow of Aspirin at the outlet and the Time extinction coefficient
> The model should be able to take raw HNMR or FTIR data and convert it into a standard spectrum able to be interpreted by analysts.
> The model could include additional visualization in the form of an animation via the Pygame module.
> The model won't consider non-steady state, non-idealized reactors or recycle streams, additionally the model is only focused on material flows through the CSTR stage of the industrial production so won't consider any flows outside the CSTR stage.

## Limitations
> A few limitations can be pointed out such as:
> > - Lack of real world kinetic data therefore inaccurate constants must be assumed, or alternatively constants from similar esterification reactions should be taken.
> > - Considering only idealized conditions makes the model inaccurate to the majority of real world scenarios.
> > - Heat transfers and Mass transfers are not considered.
> > - Maxwell-Boltzman distribution is only considered in one plane of motion.

### About the Creator
> This project was a self assigned endevour aimed at developing understanding and skills in python, Data Visualization, Physical Chemistry, Organic Chemistry, Engineering Concepts, Project Management, Research Skills, Problem Solving aswell as Documentation and Reporting.
> The creator of the project is an undergraduate student (A. G. Rios Frattaroli) at the University of Groningen in The Netherlands studying Bsc Chemical Engineering.
> LinkedIn :

"""
