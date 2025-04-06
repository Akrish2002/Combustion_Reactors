import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#import PFR.constants as c
#import PFR.util_fns as fns

import constants as c
import util_fns as fns

def PFR(gas, T_i=c.T_i, P_i=c.P_i, phi=c.phi, u_i=c.u_i, area=c.area, length=c.length, n=2000):
    
    gas = fns.set_mixture(gas, T_i, P_i, phi) 

    #Assumption --> Mass flow rate stays constant through the PFR 
    #Will use it to recalculate velocity with changing density
    mass_flow_rate = u_i * gas.density * area 

    #Reactor
    #Pseudo-PFR
    reactor = ct.IdealGasConstPressureReactor(gas)
     
    #Reactor network to advance initialized reactor through time
    sim = ct.ReactorNet([reactor])

    #Total residence time
    t = length/u_i
    dt = t/n
    
    #Array of time steps we are simulating
    t_sim = (np.arange(n) + 1) * dt #Skipping the first timestep, t = 0
    u_sim = np.zeros_like(t_sim)
    z_sim = np.zeros_like(t_sim)
    states = ct.SolutionArray(gas) #Container to store gas's thermodynamic property at each time step
    
    for i, t in enumerate(t_sim):
        sim.advance(t)
        
        u_sim[i] = mass_flow_rate/(reactor.thermo.density * area)
        z_sim[i] = z_sim[i - 1] + u_sim[i] * dt #Computing the new axial location in the reactor
        states.append(reactor.thermo.state)
        
    return z_sim, t_sim, states 

def computesoln(gas, T_list, phi_list, P_list, cycle):
    
    #Q1.a
    plt.figure()
    for P_i in P_list:
        for phi_i in phi_list:
            print("--Running case: P =", P_i, "phi =", phi_i, "\n") 
            for T_i in T_list:
                _, t_sim, states = PFR(gas, T_i, P_i, phi_i)
                T_profile = states.T
                plt.plot(t_sim, T_profile, label = f"T₀ = {T_i} K, phi₀ = {phi_i}, P₀ = {P_i} atm")
         
            plt.xlabel("Time [s]")
            plt.ylabel("Temperature [K]")
            plt.grid(True)
            
            filename = f"T_vs_time_P{int(P_i/ct.one_atm)}atm_phi{phi_i}.png"
            plt.savefig(filename)
            plt.close()
     
            if(cycle == 1): return   

T_list = [1200, 1400, 1600, 1800]
phi_list = [0.5, 1.0, 1.5]
P_list = [1, 10]
cycle = 0

gas = fns.load_mechanism("mech-FFCM1.yaml")
computesoln(gas, T_list, phi_list, P_list, cycle)    
