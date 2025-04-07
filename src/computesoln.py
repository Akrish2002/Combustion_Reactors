import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#import PFR.constants as c
#import PFR.util_fns as fns

import constants as c
import util_fns as fns
import reactors as r

def computesoln(gas, T, P, phi, Q, cycle=0):
    
    #Q1.a.1
    if(Q == "1.a.1"):
        plt.figure()
        for P_i in P_list:
            for phi_i in phi_list:
                print("--Running case: P =", P_i, "phi =", phi_i) 
                for T_i in T_list:
                    _, t_sim, states = r.PFR(gas, T_i, P_i, phi_i)
                    T_profile = states.T
                    plt.plot(t_sim, T_profile, label = f"T₀ = {T_i} K, phi₀ = {phi_i}, P₀ = {P_i} atm")
             
                plt.xlabel("Time [s]")
                plt.ylabel("Temperature [K]")
                plt.grid(True)
                
                filename = f"T_vs_time_P{int(P_i/ct.one_atm)}atm_phi{phi_i}.png"
                plt.savefig(filename)
                plt.close()

                if(cycle == 0): return   
    
    #Q1.a.2
    elif(Q == "1.a.2"):
        plt.figure()
        z_sim, _, states = r.PFR(gas, T, P, phi)
        
        #TvsZ
        T_profile = states.T
        mask = z_sim < 0.025
        z_sim_T = z_sim[mask]

        plt.plot(z_sim_T, T_profile[mask])
        plt.xlabel("Location [m]")
        plt.ylabel("Temperature [K]")
        plt.grid(True)
              
        filename = "T_vs_Z.png"
        plt.savefig(filename)
        print("\n--Temp vs Location profile has been plotted!")
        plt.close()

        #SpeciesvsZ
        mask = z_sim < 0.0075
        z_sim_masked = z_sim[mask]
        X_CH4 = states('CH4').X[mask]
        X_O2  = states('O2').X[mask]
        X_OH  = states('OH').X[mask]
        X_H2  = states('H2').X[mask]
        X_CO  = states('CO').X[mask]
        X_H2O = states('H2O').X[mask]
        X_CO2 = states('CO2').X[mask]
        
        plt.figure(figsize=(10, 8))
        plt.plot(z_sim_masked, X_CH4, label='CH₄')
        plt.plot(z_sim_masked, X_O2,  label='O₂')
        plt.plot(z_sim_masked, X_OH,  label='OH')
        plt.plot(z_sim_masked, X_H2,  label='H₂')
        plt.plot(z_sim_masked, X_CO,  label='CO')
        plt.plot(z_sim_masked, X_H2O, label='H₂O')
        plt.plot(z_sim_masked, X_CO2, label='CO₂')
        plt.xlabel("Axial Location [m]")
        plt.ylabel("Mole Fraction")
        
        filename = "Species_vs_Location.png"
        plt.savefig(filename)
        print("--Species Conc vs Location profile has been plotted!\n")
        plt.close()

    #Q1.b
    elif(Q == "1.b"):
        oxidizer_list = [{'O2': 1.0, 'N2': 3.76}, {'O2': 1.0}, {'O2': 0.207, 'N2': 0.79, 'OH': 0.003}, {'O2': 0.207, 'N2': 0.79, 'H2O': 0.003}]
  
        plt.figure()

#       _, t_sim, states = r.PFR(gas, T, P, phi, ox=oxidizer_list[0])
#       T_profile = states.T
#
#       plt.plot(t_sim, T_profile)
#       print("--Q1.b gas mixture has been plotted for:", oxidizer_list[0])

        for ox in oxidizer_list:
  
            #Running PFR
            _, t_sim, states = r.PFR(gas, T, P, phi, ox=ox)
            T_profile = states.T
            
            plt.plot(t_sim, T_profile)
            print("--Q1.b gas mixture has been plotted for:", ox)
        
        plt.xscale("log")
        plt.savefig("TvsTime_ox.png")

    elif(Q == "2"):
        r.WSR(gas, T, P, phi)

gas = fns.load_mechanism("mech-FFCM1.yaml")

#Q1.a.1
T_list = [1200, 1400, 1600, 1800]
phi_list = [0.5, 1.0, 1.5]
P_list = [ct.one_atm, 10*ct.one_atm]
cycle = 1
#computesoln(gas, T_list, P_list, phi_list, "1.a.1", cycle)    

#Q1.a.2
#computesoln(gas, 1400, ct.one_atm, 1.5, "1.a.2")

#Q1.b
#computesoln(gas, 1400, ct.one_atm, 0.5, "1.b")

#Q2
gas = ct.Solution('gri30.yaml')
computesoln(gas, 900, ct.one_atm, 0.85, "2")

