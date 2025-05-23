import numpy as np
import json



def calculator_mf_weight_1d():
    import numpy as np
    import json
    from scipy.integrate import quad
    import calculator_pair_potential_custom
    
    
    # File paths
    
    json_file_particles_interactions = "input_data_particles_interactions_parameters.json"
    json_file_simulation_thermodynamics = "input_data_simulation_thermodynamic_parameters.json"
    r_space_file = "supplied_data_r_space.txt"
    output_file = "supplied_data_bulk_mue_rho_r_space.txt"



    # Load thermodynamic properties
    with open(json_file_simulation_thermodynamics, "r") as file:
        data_thermodynamic = json.load(file)
    total_rho = data_thermodynamic["simulation_thermodynamic_parameters"]["rho"]
    temperature =  data_thermodynamic["simulation_thermodynamic_parameters"]["temperature"]


    # Load interaction properties
    with open(json_file_particles_interactions, 'r') as file:
        data_interactions = json.load(file)
    interactions = data_interactions["particles_interactions_parameters"]["interactions"]
    data_species = data_interactions["particles_interactions_parameters"]
    species = {k: v["rho_frac"] * total_rho for k, v in data_species["species"].items()}  # Calculate rho for each species

    
    
    import numpy as np
    import sympy as sp
    from sympy import log, diff
    from scipy import integrate
    from scipy.special import j0


    def hankel_transform(interaction_type, k, r_min, r_max, epsilon, sigma):
        integral, _ = integrate.quad(lambda r: 4 * np.pi * r**2 * interaction_potential(r, epsilon, sigma, interaction_type) * j0(k * r), r_min, r_max, limit=10000)
        return integral


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
                return -4 * epsilon * ((sigma / r)**48 - (sigma / r)**24)
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

    def interaction_potential_r_1d(r, epsilon, sigma, interaction_type):
        
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
                return -4 * epsilon * ((sigma / r)**48 - (sigma / r)**24)
            else:
                return 0        
                
        elif interaction_type == "gs":
            return epsilon * sigma**2 * np.pi * np.exp(-((r /sigma)**2))
            
        elif interaction_type == "yk":
            kappa = 1.0 / sigma
            return epsilon * np.exp(-kappa * r) / r if r != 0 else 0
        
        elif interaction_type == "hc":
            return 0
        
        else:
            return 0
        
    
    
    def interaction_potential_k_1d(k, epsilon, sigma, interaction_type):
        
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
                return -4 * epsilon * ((sigma / r)**48 - (sigma / r)**24)
            else:
                return 0        
                
        elif interaction_type == "gs":
            return epsilon * sigma**3 * np.pi**(3/2) * np.exp(-((2*np.pi* kx *sigma)**2/4.0))
            
        elif interaction_type == "yk":
            kappa = 1.0 / sigma
            return epsilon * np.exp(-kappa * r) / r if r != 0 else 0
        
        elif interaction_type == "hc":
            return 0
        
        else:
            return 0
        
    
    
    
    
    
    
    
    


    
    # for mean field treatment
    epsilonij = [[0.0]*len(species) for _ in range(len(species))]
    sigmaij = [[1.0]*len(species) for _ in range(len(species))]
    interaction_type_ij = [["gs"]*len(species) for _ in range(len(species))]

   
    # for hard core treatment
    sigmai_p = [0.0]*len(species)
    flag = [0]*len(species)

    
    rhos = []
    
    i = 0 
    j = 0
    grand_rosenfeld_flag = 0
    
    for species_type, rho_value in species.items():
    
        rhos.append(species[species_type])
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
            
    
    
   
    k_space_file_path = 'supplied_data_k_space.txt'
    k_space = np.loadtxt(k_space_file_path)
    k_space = np.array(k_space)
    
    r_space_file_path = 'supplied_data_r_space.txt'
    r_space = np.loadtxt(r_space_file_path)
    r_space = np.array(r_space)
    
    r = r_space[:, 0]
    
    k_values = k_space[:, 0]
    
    # Function to calculate weight function in k-space
    
    weight_functions = []
    
    
    i = 0
    for p_species, rho in species.items():

        weight_function_in = []
        j = 0
        for other_species, other_rho in species.items():
            
            if (j>=i):
                weight_function = []
                
                '''
                
                '''
                epsilon = epsilonij[i][j]
                sigma =  sigmaij[i][j]
                interaction_type = interaction_type_ij[i][j]

                
                if interaction_type == "gs" :
                    for kx, ky, kz in k_space:
                        weight_vector = [kx, ky, kz]
                        
                        
                        
                        
                        k_value = kx
                        mod_k = np.sqrt(k_value*k_value)
                        
                        
                       
                        weight  =  interaction_potential_k_1d(kx, epsilon, sigma, interaction_type)
                        
                        
                        value = complex (weight) 
                        
                        weight_function.append(value)
                else :
                
                    
                    r_min = 0.0
                    r_max = 5*sigma
                    Vk = [hankel_transform(interaction_type, k, r_min, r_max, epsilon, sigma) for k in k_values]
                    
                    
                    weight_function = [complex(Vk[i]) if (epsilon * Vk[i] > 0) else complex(0) for i in range(len(Vk))]
                   
                        
                    for it in range (len(weight_function)):
                        if abs(k_values[it])> 10 :
                            weight_function[it] = 0  


               
            
            j = j + 1
            
            weight_function_in.append(np.array(weight_function))
            
        weight_functions.append(weight_function_in)
        i = i + 1 
                
                
        

    # Function to export weight functions to files

    
    i = 0
    for p_species, rho in species.items():
        
        j = 0
        for other_species, other_rho in species.items():
            if (j >= i):
            
                file_name = f"supplied_data_weight_MF_k_space_{p_species}{other_species}.txt"
                np.savetxt(file_name, weight_functions[i][j])
        
            j = j+1
        i = i + 1

    print("\n\n\n ... Mean field weight have been calculated and exported in the file ...")
#calculator_mf_weight_1d()
