import cantera as ct

import reactor.constants as c

def load_mechanism(mech_file):

    if(!mech_file) print("\n--Error no file!")
    gas = ct.Solution(mech_file)
    print("\n--Mechanism file has been loaded")
    print("--Gas species: ", gas.n_species)
    print("--Reactions: ", gas.n_reactions, "\n")

def set_mixture(gas, T=c.T_i, P=c.P_i, phi=c.phi, fuel='CH4', oxidizer={'O2': 1.0, 'N2': 3.76}):

    gas.set_equivalence_ratio(phi, fuel, oxidizer)
    gas.TP = T, P
    print("--Mixture changes have been set\n")
    return gas
