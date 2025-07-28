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
    return

#Classical model based on translational energy
def Maxwell_Boltzman_distribution(mass,temperature, v_max=6000): #kg, K, m/s, J
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

def Energy_Distribution(mol=1, mass=2.21e-25, T=295): #salicylic acid kg/mol
    V_min = 1000 # Activation energy used as kinetic to find velocity Ek = 1/2 mv^2
    C = mass/(2*kb*T)
    N_total = mol*Av
    def f(v):
        dist = (4*pi*v**2)*((C/pi)**(1.5))*math.exp(-C*v**2)
        return dist

    integral, _ = quad(f, V_min, 6000)
    N = N_total * integral
    Percentage_yield = ( N/N_total ) * 100
    print(V_min,N_total, integral,N,Percentage_yield)
    return Percentage_yield, N

def Cost_Graph(mol=1, Temp=295):
    return

def salycilic_acid(mass=1,temperature=25):
    density = 1.44 #kg/L
    Molecular_weight = 0.13812 #kg/mol
    mol = (0.13812/mass) #mol
    cost = 80.55 * mass #€
    return mol, Molecular_weight, cost

def acetic_anhydride(mass=1,temperature=25):
    density = 1.08 #kg/L
    Molecular_weight = 0.10209 #kg/mol
    mol = (0.10209/mass) #mol
    cost = 138.38 * mass #€
    return mol, Molecular_weight, cost

def CSTR(mol_inflow=1, mol_outflow=1, T=25, accumulation_not_outflow=True, steady_state=True, accumulation=0): #mol/s, mol/s
    # dn/dt = Σφ_in - Σφ_out + rV
    #dn/dt = Σφ_in - Σφ_out + k[Ca]V --> dn/dt = Σφ_in - Σφ_out + A*n*exp(-Ea/RT)
    R = 8.31
    net_flow = mol_inflow - mol_outflow
    Ea = 75312 # j/mol Assuming same mechanism for all changes in condition
    A = 1e+11 # Assumption based on similar Liquid Phase Esterificaiton Reactions
    n, MW, Cost = salicylic_acid()
    if accumulation_not_outflow == True and steady_state == False:
        if accumulation != 0:
            accumulation = net_flow + A*n*math.exp(-Ea/(R*T))
            return accumulation #mol/s
    elif accumulation_not_outflow == False and steady_state == True:
        mol_outflow = mol_inflow + A*n*math.exp(-Ea/(R*T))
        return mol_outflow
    else:
        print("There was an error")

def Batch_reactor(mol_inflow=1, accumulation=0, mol_outflow=1, generation=0):
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
        axs[0][0].set_title(f"T = {temp:.0f} K \n Most Probable Velocity = {v_p:.1f} ms⁻¹ \n N% Above Φ = {N:.10f} %")
        canvas.draw()

    slider = tk.Scale(root, from_=200, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp)
    slider.set(273)
    slider.pack(pady=20)

    btn = tk.Button(root, text="Open New Window", command=open_new_window)
    btn.pack(pady=20)
    root.mainloop()
