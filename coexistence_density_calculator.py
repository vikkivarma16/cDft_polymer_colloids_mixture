# density functional minimizer/bulk rho and mue generator...

# this part of the code generates the bulk values of the density and the chemical potential in uniform way for each r space point supplied by the r-space files for all kind of the interactions specified in the interaction data json data profile... 

import numpy as np
import sympy as sp
from sympy import log, diff
from scipy import integrate
from scipy.special import j0

def  data_reader():
    
    import numpy as np
    import json
    from scipy.integrate import quad
    import calculator_pair_potential_custom
    
    json_file_particles_interactions = "input_data_particles_interactions_parameters.json"
    json_file_simulation_thermodynamics = "input_data_simulation_thermodynamic_parameters.json"
    r_space_file = "supplied_data_r_space.txt"
    output_file = "supplied_data_bulk_mue_rho_r_space.txt"



    # Load thermodynamic properties
    with open(json_file_simulation_thermodynamics, "r") as file:
        data_thermodynamic = json.load(file)
    temperature =  data_thermodynamic["simulation_thermodynamic_parameters"]["temperature"]


    # Load interaction properties
    with open(json_file_particles_interactions, 'r') as file:
        data_interactions = json.load(file)
    interactions = data_interactions["particles_interactions_parameters"]["interactions"]
    data_species = data_interactions["particles_interactions_parameters"]
    species = {k: v["rho_frac"]  for k, v in data_species["species"].items()}  # Calculate rho for each species


    # for mean field treatment
    epsilonij = [[0.0]*len(species) for _ in range(len(species))]
    sigmaij = [[1.0]*len(species) for _ in range(len(species))]
    interaction_type_ij = [["gs"]*len(species) for _ in range(len(species))]

   
    # for hard core treatment
    sigmai_p = [0.0]*len(species)
    flag = [0]*len(species)

    
    i = 0 
    j = 0
    grand_rosenfeld_flag = 0
    
    for species_type, rho_value in species.items():
        
        # Initialize total chemical potential for the current species
        j = i
        for other_species, rho_other in species.items():
            # Determine the interaction data key
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["primary"]:
                interaction_data = interactions["primary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["primary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Unpack interaction parameters
            epsilon = interaction_data["epsilon"]
            sigma_ij = interaction_data["sigma"]
            interaction_type = interaction_data["type"]
            cutoff = interaction_data["cutoff"]

            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                    flag[i] = 1
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
            elif interaction_type == "ghc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
                    
            else :
                interaction_type_ij[i][j] = interaction_type
                epsilonij[i][j] = epsilon
                sigmaij[i][j] = sigma_ij
            
            j = j+1

        j = i 
        for other_species, rho_other in species.items():
            # Determine the interaction data key
            
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["secondary"]:
                interaction_data = interactions["secondary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["secondary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Unpack interaction parameters
            epsilon = interaction_data["epsilon"]
            sigma_ij = interaction_data["sigma"]
            interaction_type = interaction_data["type"]
            cutoff = interaction_data["cutoff"]

            
            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                    flag[i] = 1
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
            elif interaction_type == "ghc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
                    
            else :
                interaction_type_ij[i][j] = interaction_type
                epsilonij[i][j] = epsilon
                sigmaij[i][j] = sigma_ij
                
            j = j+1
                
        
        j = i
        for other_species, rho_other in species.items():
            # Determine the interaction data key
            
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["tertiary"]:
                interaction_data = interactions["tertiary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["tertiary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                    flag[i] = 1
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
            elif interaction_type == "ghc":
                grand_rosenfeld_flag = 1
                if species_type == other_species:
                # Approximate hard-core chemical potential contribution
                    sigmai_p[i] = sigma_ij
                else:
                    print("\n ...wrong parameters have been defined... \n")
                    exit(0)
                    
            else :
                interaction_type_ij[i][j] = interaction_type
                epsilonij[i][j] = epsilon
                sigmaij[i][j] = sigma_ij
            j = j+1
        
        #print("mue before ideal", mue_total/temperature)
        # Store the total chemical potential for the species
        i = i+1
        
    return species, epsilonij, sigmaij, interaction_type_ij, sigmai_p, grand_rosenfeld_flag, flag


def canonical_coexistence(species, epsilonij, sigmaij, interaction_type_ij, sigmai_p, grand_rosenfeld_flag, flag):
    import numpy as np
    import json
    from scipy.integrate import quad
    import calculator_pair_potential_custom
    
    
    
    import numpy as np
    import sympy as sp
    from sympy import log, diff
    from scipy import integrate
    from scipy.special import j0


    def hankel_transform(interaction_type, k, r_min, r_max, epsilon, sigma):
        integral, _ = integrate.quad(lambda r: 4 * np.pi * r**2 * interaction_potential(r, epsilon, sigma, interaction_type) * j0(k * r), r_min, r_max, limit=10000)
        return integral


    def free_energy_mean_field( epsilonij, sigmaij, interaction_type_ij, densities):
        k_values = np.linspace(0.00001, 5, 100)
        
        vij = []
        for i in range(len(epsilonij)):
            temp = []
            for j in range(len(epsilonij)):
                epsilon = epsilonij[i][j]  
                sigma = sigmaij[i][j]  
                r_min = 0
                r_max = 5*sigma
                interaction_type = interaction_type_ij[i][j]  

                Vk = [hankel_transform(interaction_type, k, r_min, r_max, epsilon, sigma) for k in k_values]
                temp.append(Vk[0])
                #print(Vk[0])
                
            vij.append(temp)

        

        
        
        fideal = densities[0] * sp.log(densities[0]) - densities[0]
        
        fsub = densities[0]
        for i in range(1, len(epsilonij)):
            fideal += densities[i] * sp.log(densities[i]) - densities[i]
            fsub += densities[i]
            
      

        ftotal = fideal
        for i in range(len(epsilonij)):
            for j in range(len(epsilonij)):
                ftotal += 0.5 * vij[i][j] * densities[i] * densities[j]
                
       
        mue_mf = [diff(ftotal, densities[i]) for i in range(len(epsilonij))]
        partial_pressure_mf =-ftotal + sum(densities[i]*mue_mf[i] for i in range(len(epsilonij)))
        
        #mue_mf = [sp.lambdify(densities, mue[i], 'numpy') for i in range(len(epsilonij))]
        #ftotal = sp.lambdify(densities, ftotal, 'numpy')
        #partial_pressure_mf = sp.lambdify(densities, partial_pressure_mf, 'numpy')

        return mue_mf, partial_pressure_mf, ftotal


    def hard_core_approach(sigmai,  flag, densities):
        measures = [[1, sigma/2, np.pi*sigma**2, sigma**3 * np.pi/6] for sigma in sigmai]
        
        

        etas = [sp.symbols(f"eta_{i}") for i in range(len(sigmai))]
        variables = [[densities[i] * measures[i][j] for j in range(4)] for i in range(len(sigmai))]

        fac1 = 1 - sum(etas)
        fac2 = sum(etas[i] for i in range(len(sigmai)) if flag[i] == 1)
        
        
        phi0 = fac1 * sp.log(1 - fac2) + fac2
        
     
        
        diff_1 = [sp.diff(phi0, etas[i]) for i in range(len(sigmai))]
        diff_2 = [[sp.diff(phi0, etas[i], etas[j]) for j in range(len(sigmai))] for i in range(len(sigmai))]
        diff_3 = [[[sp.diff(phi0, etas[i], etas[j], etas[k]) for k in range(len(sigmai))] for j in range(len(sigmai))] for i in range(len(sigmai))]

      
        
        

        phi1 = sum(variables[i][0] * diff_1[i] for i in range(len(sigmai)))
        phi2 = sum(variables[i][1] * variables[j][2] * diff_2[i][j] for i in range(len(sigmai)) for j in range(len(sigmai)))
        phi3 = sum((1/(24*np.pi)) * variables[i][2]*variables[j][2]*variables[k][2] * diff_3[i][j][k] for i in range(len(sigmai)) for j in range(len(sigmai)) for k in range(len(sigmai)))
        
        
        
         
        
        fhc = phi1 + phi2 + phi3

        
        
        
        for i in range(len(sigmai)):
            fhc = fhc.subs(etas[i], variables[i][3])
            phi1 = phi1.subs(etas[i], variables[i][3])
            phi2 = phi2.subs(etas[i], variables[i][3])
            phi3 = phi3.subs(etas[i], variables[i][3])
        
        vphi1 = sp.lambdify(densities, phi1, 'numpy')
        vphi2 = sp.lambdify(densities, phi2, 'numpy')
        vphi3 = sp.lambdify(densities, phi3, 'numpy')
        
        
        
        
        mue_hc = [diff(fhc, densities[i]) for i in range(len(sigmai))]
        partial_pressure_hc = -fhc + sum(densities[i]*mue_hc[i] for i in range(len(epsilonij)))
        
        
        #mue_hc = [sp.lambdify(densities, mue[i], 'numpy') for i in range(len(sigmai))]
        #fhc = sp.lambdify(densities, fhc, 'numpy')
        #partial_pressure_hc = sp.lambdify(densities, partial_pressure_hc, 'numpy')
        
        return mue_hc, partial_pressure_hc, fhc


    def interaction_potential(r, epsilon, sigma, interaction_type):
        
        if interaction_type == "wca":
            if r < 2**(1/6) * sigma:
                return -epsilon
            elif r < 5 * sigma:
                return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
            else:
                return 0
                
        if interaction_type == "mie":
            if r < sigma:
                return 0
            elif r < 5 * sigma:
                return 4 * epsilon * ((sigma / r)**48 - (sigma / r)**24)
            else:
                return 0        
                
        elif interaction_type == "gs":
            return epsilon * np.exp(-((r / sigma)**2))
            
        elif interaction_type == "yk":
            kappa = 1.0 / sigma
            return epsilon * np.exp(-kappa * r) / r if r != 0 else 0
        
        elif interaction_type == "hc":
            return 0
        
        else:
            return 0
        
        
    
    if (grand_rosenfeld_flag == 1):
        for i in range(len(species)):
            if( sigmai_p[i] <0.0001):
            
                print(" ------------ ERROR FLAG --------------")
                print("... not all the parameters are defined ... if you are using one hard core then it will be hardcore also for the rest of the particles so you need to put some finite size even to the point particles so there will always be some hard core interaction with each particles...\n\n\n")
                print("... so at least put some ghc potential")
                exit(0)
                
    
    
    for i in range(len (epsilonij)):
        for j in range(i+1, len (epsilonij)):
            epsilonij[j][i] = epsilonij[i][j]
            sigmaij[j][i] = sigmaij[i][j]
            interaction_type_ij[j][i] = interaction_type_ij[i][j]
            
            
   
               
    densities = [sp.symbols(f"rho_{i}") for i in range(len(epsilonij))]
    
    
    # Step 1: Mean-field contributions
    mue_mf, partial_pressure_mf, ftotal = free_energy_mean_field(epsilonij, sigmaij, interaction_type_ij, densities)

    # Initial values
    mue = mue_mf.copy()  # Make sure to copy to avoid modifying mue_mf directly
    pressure = partial_pressure_mf



    sigmai = sigmai_p
    # Step 2: Hard-core corrections
    if grand_rosenfeld_flag == 1:
        mue_hc, partial_pressure_hc, fhc = hard_core_approach(sigmai, flag, densities)

        # Add hard-core corrections to each component
        for i in range(len(mue)):
            mue[i] += mue_hc[i]

        ftotal += fhc
        pressure += partial_pressure_hc


    # Step 3: Return final quantities
    return mue, pressure, ftotal, densities

    
    
    

species, epsilonij, sigmaij, interaction_type_ij, sigmai_p, grand_rosenfeld_flag, flag =  data_reader()
mue, pressure, ftotal, densities = canonical_coexistence(species, epsilonij, sigmaij, interaction_type_ij, sigmai_p, grand_rosenfeld_flag, flag)


densities_1 = [sp.symbols(f"rhoa_{i}") for i in range(len(sigmaij))]
densities_2 = [sp.symbols(f"rhob_{i}") for i in range(len(sigmaij))]
mue_a = []
mue_b = []
for j in range(len(mue)):
    updated_mue_1=mue[j]
    updated_mue_2=mue[j]
    for i in range(len(sigmaij)):
        updated_mue_1=updated_mue_1.subs(densities[i], densities_1[i])
        updated_mue_2=updated_mue_2.subs(densities[i], densities_2[i])
    mue_a.append(updated_mue_1)
    mue_b.append(updated_mue_2)
    
pressure_a = pressure
pressure_b = pressure
for i in range(len(sigmaij)):
    pressure_a=pressure_a.subs(densities[i], densities_1[i])
    pressure_b=pressure_b.subs(densities[i], densities_2[i])




mue_a_value = [sp.lambdify(densities_1, mue_a[i], 'numpy') for i in range(len(sigmai_p))]
mue_b_value = [sp.lambdify(densities_2, mue_b[i], 'numpy') for i in range(len(sigmai_p))]

pressure_a_value = sp.lambdify(densities_1, pressure_a, 'numpy')
pressure_b_value = sp.lambdify(densities_2, pressure_b, 'numpy')



def particle_numbers():

    p1 = 0.13
    p2 = 0.13
    c = 0.09

    return  p1, p2, c



def coexistence_residual(vars):
    
    p1, p2, c = particle_numbers()
    
    rhot_a, rhot_b, x_a, x_b, y_a, y_b, p = vars
    
    rhoa_1 = rhot_a * (1- y_a) * (1.0 - x_a)
    rhoa_2 = rhot_a * (1- y_a) * (x_a)
    rhoa_3 = rhot_a * y_a 
    
    rhob_1 = rhot_b * (1- y_b) * (1.0 - x_b)
    rhob_2 = rhot_b * (1- y_b) * (x_b)
    rhob_3 = rhot_b * y_b 
    
    mu1_a = mue_a_value[0](rhoa_1, rhoa_2, rhoa_3)
    mu2_a = mue_a_value[1](rhoa_1, rhoa_2, rhoa_3)
    mu3_a = mue_a_value[2](rhoa_1, rhoa_2, rhoa_3)
    
    pressure_a = pressure_a_value(rhoa_1, rhoa_2, rhoa_3)
    
    mu1_b = mue_b_value[0](rhob_1, rhob_2, rhob_3)
    mu2_b = mue_b_value[1](rhob_1, rhob_2, rhob_3)
    mu3_b = mue_b_value[2](rhob_1, rhob_2, rhob_3)
    
    pressure_b = pressure_b_value(rhob_1, rhob_2, rhob_3)

    eq1 = mu1_a - mu1_b
    eq2 = mu2_a - mu2_b
    eq3 = mu3_a - mu3_b
    eq4 = pressure_a - pressure_b
    eq5 = p * rhoa_1 + (1-p) * rhob_1 - p1
    eq6 = p * rhoa_2 + (1-p) * rhob_2 - p2
    eq7 = p * rhoa_3 + (1-p) * rhoa_3 - c
    
    return eq1, eq2, eq3, eq4, eq5, eq6, eq7
    
    
import numpy as np
from scipy.optimize import root

def random_initial_guess():
    rhot_a = np.random.uniform(0.0, 1.0)
    rhot_b = np.random.uniform(0.0, 1.0)
    x_a = np.random.uniform(0.0, 0.4)
    x_b = np.random.uniform(0.6, 1.0)
    y_a = np.random.uniform(0.0, 1.0)
    y_b = np.random.uniform(0.0, 1.0)
    p = np.random.uniform(0.0, 1.0)
    return [rhot_a, rhot_b, x_a, x_b, y_a, y_b, p]

def solve_phase_coexistence(n_attempts=50):
    for i in range(n_attempts):
        guess = random_initial_guess()
        sol = root(coexistence_residual, guess, method='hybr')

        if sol.success:
            rhot_a, rhot_b, x_a, x_b, y_a, y_b, p = sol.x

            # Check if solution is within bounds and disjoint regions
            if (
                0 <= rhot_a <= 1 and 0 <= rhot_b <= 1 and
                0.0 <= x_a <= 0.4 and 0.6 <= x_b <= 1.0 and
                0 <= y_a <= 1 and 0 <= y_b <= 1 and
                0 <= p <= 1
            ):
                print("✅ Coexistence solution found:")
                print(f"rhot_a = {rhot_a:.4f}, x_a = {x_a:.4f}, y_a = {y_a:.4f}")
                print(f"rhot_b = {rhot_b:.4f}, x_b = {x_b:.4f}, y_b = {y_b:.4f}")
                print(f"p = {p:.4f}")
                return sol.x  # Return the successful solution

    print("❌ No valid coexistence solution found after multiple attempts.")
    return None



solution = solve_phase_coexistence(n_attempts=100)


print(solution)


rhoa_1 = solution[0] * (1- solution[4]) * (1.0 - solution[2])
rhoa_2 = solution[0] * (1- solution[4]) * (solution[2])
rhoa_3 = solution[0] * solution[4]


rhob_1 = solution[1] * (1- solution[5]) * (1.0 - solution[3])
rhob_2 = solution[1] * (1- solution[5]) * (solution[3])
rhob_3 = solution[1] * solution[5]


p1 = solution[6] * rhoa_1 + (1-solution[6]) * rhob_1 
p2 = solution[6] * rhoa_2 + (1-solution[6]) * rhob_2
c  = solution[6] * rhoa_3 + (1-solution[6]) * rhoa_3


print(rhoa_1, rhoa_2, rhoa_3, rhob_1, rhob_2, rhob_3, p1, p2, c)

