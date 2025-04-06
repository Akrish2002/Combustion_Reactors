import cantera as ct

def load_mechanism(mech_file):
    gas = ct.Solution(mech_file)
    print("\n--Mechanism file has been loaded")
    print("--Gas species: ", gas.n_species)
    print("--Reactions: ", gas.n_reactions, "\n")
    
    return gas

def set_mixture(gas, T, P, phi, fuel='CH4', oxidizer={'O2': 1.0, 'N2': 3.76}):
    gas.set_equivalence_ratio(phi, fuel, oxidizer)
    gas.TP = T, P
    print("--Mixture changes have been set\n")
    return gas

gas = load_mechanism("mech-FFCM1.yaml")
gas = set_mixture(gas, 1200, 1, 1.5)

    
