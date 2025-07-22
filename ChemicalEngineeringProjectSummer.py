import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import tkinter as tk
from tkinter import ttk
#Engineering Project Outline - Chemical Reactor Simulator
pi = 3.14159265
kb = 1.380649e-23

def Maxwell_Boltzman_distribution(mass,temperature, v_max=5000): #kg, K, m/s, J
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

def CSTR(mol_inflow=1, mol_outflow=1, T=25, V=1): #mol/s, mol/s
    # dn/dt = Σφ_in - Σφ_out + rV
    R = 8.31
    net_flow = mol_inflow - mol_outflow
    Ea = 75312 # j/mol Assuming same mechanism for all changes in condition
    A = 1e+11
    n, MW, Cost = salicylic_acid()
    accumulation = net_flow + A*n*math.exp(-Ea/(R*T))

    return accumulation

def Batch_reactor(mass_inflow=1, accumulation=0, mass_outflow=1, generation=0):
    return

#main window
root = tk.Tk()
root.title("System Visual")
root.geometry("1500x1000")
frame = tk.Frame(root)
frame.pack()

#4 sub plots
fig, axs = plt.subplots(2,1,figsize=(10, 7))
canvas = FigureCanvasTkAgg(fig,master=frame)
canvas.get_tk_widget().pack()

#slider
def update_temp(temp, vel=0):
    vel = float(vel)
    temp = float(temp)
    m = 1.67e-27
    distribution, velocities, v_p, f_vp = Maxwell_Boltzman_distribution(m,temp)
    axs[0].clear()
    axs[0].plot(velocities, distribution)
    axs[0].set_title(f"T = {temp:.0f} K \n Most Probable Velocity = {v_p:.1f}")
    canvas.draw()


slider = tk.Scale(root, from_=200, to=500, orient=tk.HORIZONTAL, label="Temperature (K)", command=update_temp)
slider.set(273)
slider.pack(pady=20)

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

btn = tk.Button(root, text="Open New Window", command=open_new_window)
btn.pack(pady=20)
root.mainloop()