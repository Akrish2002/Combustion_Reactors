import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

import constants as c
import util_fns as fns

def PFR(gas, T_i=c.T_i, P_i=c.P_i, phi=c.phi, fuel='CH4', ox={'O2': 1.0, 'N2': 3.76}, u_i=c.u_i, area=c.area, length=c.length, n=2000):
    
    gas = fns.set_mixture(gas, T_i, P_i, phi, ox=ox) 

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

def WSR(gas, tau=c.tau, T_i=c.T_i, P_i=c.P_i, phi=c.phi, fuel='CH4', ox={'O2': 1.0, 'N2': 3.76}):

    '''
    Following arrangement is considered for the well sttired reactor :
    Inlet reservoir [Mixture tank] --> Mass flow controller --> Stirred Reactor --> Pressure Valve --> Exhaust reservoir
    
    '''

    gas = fns.set_mixture(gas, T_i, P_i, phi)

     
    inlet = ct.Reservoir(gas)
    #gas.equilibrate('HP') 
    stirred_reactor = ct.IdealGasReactor(gas, energy='on', volume=c.V)
    exhaust = ct.Reservoir(gas)
    

    #Massflow controller
    mfr = fns.massflowrate(stirred_reactor.density, tau=tau)
    mass_flow_controller = ct.MassFlowController(upstream=inlet, 
                                                 downstream=stirred_reactor,
                                                 mdot=mfr)

    pressure_regulator = ct.PressureController(upstream=stirred_reactor,
                                               downstream=exhaust,
                                               primary=mass_flow_controller,
                                               K=1e-2) 
    sim = ct.ReactorNet([stirred_reactor])

    time = 0.0
    t_end = 10 * tau
    T1 = stirred_reactor.T
    
    while time < t_end:
        time = sim.step()
        T2 = stirred_reactor.T
        
        #Checking for steady state
        if time > tau and abs(T2 - T1)/T1 < 1e-4:
            break
        
        T1 = T2
        

    T_ss = stirred_reactor.T
    X_CO = stirred_reactor.thermo['CO'].X[0]
    X_CO2 = stirred_reactor.thermo['CO2'].X[0]
    X_NO = stirred_reactor.thermo['NO'].X[0] if 'NO' in gas.species_names else 0.0
    
    return T_ss, X_CO, X_CO2, X_NO
