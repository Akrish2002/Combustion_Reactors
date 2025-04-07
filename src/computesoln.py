import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

import constants as c
import util_fns as fns
import reactors as r

def computesoln(gas, T, P, phi, Q, cycle=0):
    
#   #Q1.a.1
#   if(Q == "1.a.1"):
#       plt.figure()
#       for P_i in P_list:
#           for phi_i in phi_list:
#               print("--Running case: P =", P_i, "phi =", phi_i) 
#               for T_i in T_list:
#                   _, t_sim, states = r.PFR(gas, T_i, P_i, phi_i)
#                   T_profile = states.T
#                   plt.plot(t_sim, T_profile, label = f"T₀ = {T_i} K, phi₀ = {phi_i}, P₀ = {P_i} atm")
#            
#               plt.xlabel("Time [s]")
#               plt.ylabel("Temperature [K]")
#               plt.grid(True)
#               
#               filename = f"T_vs_time_P{int(P_i/ct.one_atm)}atm_phi{phi_i}.png"
#               plt.savefig(filename)
#               plt.title
#               plt.close()
#
#               if(cycle == 0): return   
        
    if(Q == "1.a.1"):
        for P_i in P_list:
            for phi_i in phi_list:
                print("--Running case: P =", P_i, "phi =", phi_i)
                
                # Lists to store ignition delays and inverse initial temperatures
                tau_list_condition = []
                invT_list_condition = []
                
                # Create a new figure with two subplots side by side
                fig, axs = plt.subplots(1, 2, figsize=(14, 6))
                
                for T_i in T_list:
                    # Run the PFR simulation
                    _, t_sim, states = r.PFR(gas, T_i, P_i, phi_i)
                    T_profile = states.T
                    axs[0].plot(t_sim, T_profile, label=f"T₀ = {T_i} K")
                    
                    # Compute the temperature derivative with respect to time
                    dT_dt = np.gradient(T_profile, t_sim)
                    # Define ignition delay as the time when dT/dt is maximum
                    tau_ign = t_sim[np.argmax(dT_dt)]
                    tau_list_condition.append(tau_ign)
                    invT_list_condition.append(1000 / T_i)
                
                axs[0].set_xscale("log")
                axs[0].set_xlabel("Time [s] (log scale)")
                axs[0].set_ylabel("Temperature [K]")
                axs[0].set_title(f"Temperature vs Time (P = {int(P_i/ct.one_atm)} atm, ϕ = {phi_i})")
                axs[0].grid(True)
                axs[0].legend(title="Initial T₀")
                
                axs[1].plot(invT_list_condition, tau_list_condition, 'o-', color='tab:blue')
                axs[1].set_xlabel("1000 / T₀ [1/K]")
                axs[1].set_ylabel("Ignition Delay τ [s]")
                axs[1].set_title("Ignition Delay vs Inverse Initial Temperature")
                axs[1].set_yscale("log")
                axs[1].grid(True)
                
                plt.suptitle(f"Case: P = {int(P_i/ct.one_atm)} atm, ϕ = {phi_i}", fontsize=14)
                plt.tight_layout(rect=[0, 0, 1, 0.93])
                
                filename = f"T_vs_time_and_tau_P{int(P_i/ct.one_atm)}atm_phi{phi_i}.png"
                plt.savefig(filename)
                plt.close()
                
                if(cycle == 0):
                    return    

    elif(Q == "1.a.2"):
        # Temperature vs Axial Location
        plt.figure()
        z_sim, _, states = r.PFR(gas, T, P, phi)
        
        # Extract temperature profile and limit to the region of interest 
        T_profile = states.T
        mask = z_sim < 0.025
        z_sim_T = z_sim[mask]

        plt.plot(z_sim_T, T_profile[mask], label="Temperature")
        plt.xlabel("Axial Location [m]")
        plt.ylabel("Temperature [K]")
        plt.title("Temperature vs Axial Location")
        plt.legend()
        plt.grid(True)
        
        filename = "T_vs_Z.png"
        plt.savefig(filename)
        print("\n--Temp vs Location profile has been plotted!")
        plt.close()

        # Species vs Axial Location
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
        plt.title("Species Mole Fraction vs Axial Location")
        plt.legend(title="Species")
        plt.grid(True)
        
        filename = "Species_vs_Location.png"
        plt.savefig(filename)
        print("--Species Conc vs Location profile has been plotted!\n")
        plt.close()

    elif(Q == "1.b"):
        oxidizer_list = [
            {'O2': 1.0, 'N2': 3.76},
            {'O2': 1.0},
            {'O2': 0.207, 'N2': 0.79, 'OH': 0.003},
            {'O2': 0.207, 'N2': 0.79, 'H2O': 0.003}
        ]

        plt.figure(figsize=(10, 6))

        for ox in oxidizer_list:
            # Run the PFR simulation with the given oxidizer. 
            _, t_sim, states = r.PFR(gas, T, P, phi, ox=ox)
            T_profile = states.T

            plt.plot(t_sim, T_profile, label=f"Oxidizer: {ox}")
            print("--Q1.b gas mixture has been plotted for:", ox)

        plt.xscale("log")
        plt.xlabel("Time [s] (log scale)")
        plt.ylabel("Temperature [K]")
        plt.title("Temperature vs Time for Different Oxidizer Compositions")
        plt.legend(title="Oxidizer Mixture")
        plt.grid(True)
        plt.tight_layout()

        plt.savefig("TvsTime_ox.png")
        plt.show()

    elif(Q == "2"):

        tau_list = np.logspace(np.log10(1e8), np.log10(1e-4), num=45)  # 1 s to 0.1 ms

        # Storage for results
        T_out = []
        X_CO = []
        X_CO2 = []
        X_NO = []
     
        for tau in tau_list:
            gas = ct.Solution('gri30.yaml')
            gas = fns.set_mixture(gas, T, P, phi, "CH4", {"O2": 1.0, "N2": 3.76})
            
            try:
                T_ss, x_co, x_co2, x_no = r.WSR(gas, tau)
                print(f"τ = {tau:.5f} s → T = {T_ss:.1f} K, CO = {x_co:.2e}, CO₂ = {x_co2:.2e}, NO = {x_no:.2e}")
            except Exception as e:
                print(f"τ = {tau:.5f} s → Simulation failed: {e}")
                T_ss, x_co, x_co2, x_no = np.nan, np.nan, np.nan, np.nan
     
            T_out.append(T_ss)
            X_CO.append(x_co)
            X_CO2.append(x_co2)
            X_NO.append(x_no)
     
        # Convert to numpy arrays
        T_out = np.array(T_out)
        X_CO = np.array(X_CO)
        X_CO2 = np.array(X_CO2)
        X_NO = np.array(X_NO)
     
        # Plotting
        plt.figure(figsize=(10, 8))
     
        plt.subplot(2, 2, 1)
        plt.semilogx(tau_list, T_out, 'o-')
        plt.xlabel("Residence Time τ [s]")
        plt.ylabel("Temperature [K]")
        plt.title("T vs Residence Time")
        plt.grid(True)
     
        plt.subplot(2, 2, 2)
        plt.semilogx(tau_list, X_CO, 'o-', label="CO")
        plt.semilogx(tau_list, X_CO2, 'o-', label="CO₂")
        plt.xlabel("Residence Time τ [s]")
        plt.ylabel("Mole Fraction")
        plt.title("CO & CO₂ vs Residence Time")
        plt.legend()
        plt.grid(True)
     
        plt.subplot(2, 2, 3)
        plt.semilogx(tau_list, X_NO, 'o-')
        plt.xlabel("Residence Time τ [s]")
        plt.ylabel("NO Mole Fraction")
        plt.title("NO vs Residence Time")
        plt.grid(True)
     
        plt.tight_layout()
        plt.savefig("WSR_tau_sweep.png")
        plt.show()
        
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
#gas = ct.Solution('gri30.yaml')
computesoln(gas, 300, ct.one_atm, 0.85, "2")

