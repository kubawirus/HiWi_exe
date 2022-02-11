import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def gas_funct(mechanism, initial_state):

    # choosing mechanism of reaction
    # gas_mech = 'mech_13.yaml'

    if mechanism == 2:
        gas_mech = 'mech_13.yaml'
    elif mechanism == 1:
        gas_mech = 'gri30.yaml'
    print('Reaction mechanism selected to: ', gas_mech, "\n")

    gas = ct.Solution(gas_mech)
    # gas
    # gas()
    # Here is the initial conditions are defined
    # initial_state = 873.0, ct.one_atm, 'C3H8:10, H2:1'

    # Run a simulation with the full mechanism
    gas.TPX = initial_state
    n_spec = gas.n_species
    # gas()
    # print(gas.n_reactions)
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])

    tt = []
    TT = []
    t = 0.0
    # Rmax is the maximum relative reaction rate at any timestep. An Array of size n_reactions full of zeros is created
    Rmax = np.zeros(gas.n_reactions)
    # It does not print any reactions!!! Maybe due to 2 phases inside this mechanism. For gri the reactions are printed!
    # print(Rmax)
    while t < 0.5:
        t = sim.step()
        tt.append(1000 * t)
        TT.append(r.T)
        rnet = abs(gas.net_rates_of_progress)
        rnet /= max(rnet)
        Rmax = np.maximum(Rmax, rnet)

    plt.plot(tt, TT, label='K=All_species, R=All_reac.', color='k', lw=2)

    # Get the reaction objects, and sort them so the most active reactions are first
    R = sorted(zip(Rmax, gas.reactions()), key=lambda x: -x[0])

    # Smallest reactions rate accepted
    kmin = 10.0 ** (-15)
    # Create an empty list for reactions
    reactions = []

    # Iterating through every reaction (these are sorted that's why break after else)
    # and writing to new list all these reactions that have k > 10^-20
    for i, j in enumerate(R):
        if j[0] > kmin:
            reactions.append(j[1])
        else:
            break

    # print(reactions)

    # find the species involved in these reactions. At a minimum, include all
    # species in the reactant mixture`
    species_names = {'C3H8', 'H2'}
    for reaction in reactions:
        species_names.update(reaction.reactants)
        species_names.update(reaction.products)

    # Get the species objects
    species = [gas.species(name) for name in species_names]
    # print(species)

    # create the new reduced mechanism
    gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                       species=species, reactions=reactions, transport='mixture-averaged')
    # gas2 = ct.Solution(thermo='IdealGas', species=species, transport='mixture-averaged')
    gas_surf = ct.Solution(thermo='IdealGas', species=species, transport='mixture-averaged')
    # show gas2
    # gas2()

    # write a new txt file to review the existing species in gas2
    f = open("species_reduced.txt", "w")
    for red_spec in gas2.species_names:
        f.write(red_spec + "\n")
    f.close()

    # Re-run the ignition problem with the reduced mechanism (from there is not compulsory)
    gas2.TPX = initial_state
    n_spec2 = gas2.n_species
    r = ct.IdealGasConstPressureReactor(gas2)
    sim = ct.ReactorNet([r])

    t = 0.0

    tt = []
    TT = []
    while t < 0.5:
        t = sim.step()
        tt.append(1000 * t)
        TT.append(r.T)

    # Checking on the plot if the mechanism looks same with reduced number of reactions.
    plt.plot(tt, TT, 'r--', label='K= Reduced spec., R= Reduced react.', lw=1)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc='upper left')
    plt.title('Reduced mechanism ignition delay times\n'
              'K: number of species; R: number of reactions')
    plt.tight_layout()

    # plt.show()
    print(gas_mech, " mechanism was reduced from ", n_spec, " to ", n_spec2, " species\n" )
    return gas2, gas_mech, gas_surf