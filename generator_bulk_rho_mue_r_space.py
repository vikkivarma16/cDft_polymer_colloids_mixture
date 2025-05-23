# density functional minimizer/bulk rho and mue generator...

# this part of the code generates the bulk values of the density and the chemical potential in uniform way for each r space point supplied by the r-space files for all kind of the interactions specified in the interaction data json data profile... 

def bulk_rho_mue_r_space():
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
        
    
    temperature =  data_thermodynamic["simulation_thermodynamic_parameters"]["temperature"]
    
    species = {k: v["rho_frac"]  for k, v in data_thermodynamic["simulation_thermodynamic_parameters"]["species"].items()}  # Calculate rho for each species
    secondary_species = {k: v["secondary_rho_frac"] for k, v in data_thermodynamic["simulation_thermodynamic_parameters"]["species"].items()}  # Calculate rho for each species

    
    # Load interaction properties
    with open(json_file_particles_interactions, 'r') as file:
        data_interactions = json.load(file)
    interactions = data_interactions["particles_interactions_parameters"]["interactions"]
    data_species = data_interactions["particles_interactions_parameters"]
   

    
    
    import numpy as np
    import sympy as sp
    from sympy import log, diff
    from scipy import integrate
    from scipy.special import j0


    def hankel_transform(interaction_type, k, r_min, r_max, epsilon, sigma):
        integral, _ = integrate.quad(lambda r: 4 * np.pi * r**2 * interaction_potential(r, epsilon, sigma, interaction_type) * j0(k * r), r_min, r_max, limit=10000)
        return integral


    def free_energy_mean_field( epsilonij, sigmaij, interaction_type_ij, rhos):
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
                print(Vk[0])
                
            vij.append(temp)

        densities = [sp.symbols(f"rho{i}") for i in range(len(epsilonij))]

        
        
        fideal = densities[0] * sp.log(densities[0]) - densities[0]
        
        fsub = densities[0]
        for i in range(1, len(epsilonij)):
            fideal += densities[i] * sp.log(densities[i]) - densities[i]
            fsub += densities[i]
            
      

        ftotal = fideal
        for i in range(len(epsilonij)):
            for j in range(len(epsilonij)):
                ftotal += 0.5 * vij[i][j] * densities[i] * densities[j]
                
       
        mue = [diff(ftotal, densities[i]) for i in range(len(epsilonij))]

        muef = [sp.lambdify(densities, mue[i], 'numpy') for i in range(len(epsilonij))]

        mue_values = [f(*rhos) for f in muef]

        return mue_values


    def hard_core_approach(sigmai, rhos, flag):
        densities = [sp.symbols(f"rho_{i}") for i in range(len(sigmai))]
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
        
        mue = [diff(fhc, densities[i]) for i in range(len(sigmai))]
        muef = [sp.lambdify(densities, mue[i], 'numpy') for i in range(len(sigmai))]

        mue_values = [f(*rhos) for f in muef]

        return mue_values


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
        
    
    
    
    
    
    
    
    
    
    
    


    
    # for mean field treatment
    epsilonij = [[0.0]*len(species) for _ in range(len(species))]
    sigmaij = [[1.0]*len(species) for _ in range(len(species))]
    interaction_type_ij = [["gs"]*len(species) for _ in range(len(species))]

   
    # for hard core treatment
    sigmai_p = [0.0]*len(species)
    flag = [0]*len(species)

    
    rhos = []
    secondary_rhos = []
    
    i = 0 
    j = 0
    grand_rosenfeld_flag = 0
    
    for species_type, rho_value in species.items():
    
        rhos.append(species[species_type])
        secondary_rhos.append(secondary_species[species_type])
        
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
            
            
   
               
    valuemf = free_energy_mean_field(epsilonij, sigmaij, interaction_type_ij, rhos)
    
    
    print(valuemf) 
 
    
    sigmai = sigmai_p
    if (grand_rosenfeld_flag == 1):
        valuehc = hard_core_approach(sigmai, rhos, flag)
        
    
    
    print(valuehc)
    
    bulk_mue =[]
    
    
    for i in range(len(species)):
        bulk_mue.append(valuemf[i] + valuehc[i])
        
        
    secondary_valuemf = free_energy_mean_field(epsilonij, sigmaij, interaction_type_ij, secondary_rhos)
    
    print(secondary_valuemf)
    
    
    sigmai = sigmai_p
    if (grand_rosenfeld_flag == 1):
        secondary_valuehc = hard_core_approach(sigmai, secondary_rhos, flag)
    
    print(secondary_valuehc)
    
    bulk_mue =[]
    
    
    for i in range(len(species)):
        bulk_mue.append(valuemf[i] + valuehc[i])
    
    
    secondary_bulk_mue =[]
    
    
    for i in range(len(species)):
        secondary_bulk_mue.append(secondary_valuemf[i] + secondary_valuehc[i])
    
    
    # Load r-space data
    r_space_data = np.loadtxt(r_space_file)
    
    
    
    
    warning  = 0.0
    
    for i in range(len(bulk_mue)):
        warning  =  abs(bulk_mue[i]- secondary_bulk_mue[i]) + warning
        
        print(bulk_mue[i], secondary_bulk_mue[i], "\n\n")
        
    if (warning >1.5):
        print("Warning flag values: ", warning, "... exiting the computation, too much off densities")
        exit(0)
    else:
        
        '''    
        p_rho  = []
        for species_type, rho_value in species.items():
            p_rho.append(rho_value)
            
        s_rho  = []
        for species_type, rho_value in secondary_species.items():
            s_rho.append(rho_value)
        
        for i in range(len(bulk_mue)):
            if (p_rho[i]>s_rho[i]):
                bulk_mue[i] =  bulk_mue[i] + 0.01* bulk_mue[i]
                secondary_bulk_mue[i] =  secondary_bulk_mue[i] - 0.01* secondary_bulk_mue[i]
            else:
                bulk_mue[i] =  bulk_mue[i] - 0.01* bulk_mue[i]
                secondary_bulk_mue[i] =  secondary_bulk_mue[i] + 0.01* secondary_bulk_mue[i]
        
        for i in range(len(bulk_mue)):
            if (bulk_mue[i]>secondary_bulk_mue[i]):
                bulk_mue[i] = secondary_bulk_mue[i]
            else: 
                secondary_bulk_mue[i]= bulk_mue[i]
        
        for i in range(len(bulk_mue)):
            if (p_rho[i]>s_rho[i]):
                bulk_mue[i] =  bulk_mue[i] + 0.01* bulk_mue[i]
                secondary_bulk_mue[i] =  secondary_bulk_mue[i] - 0.01* secondary_bulk_mue[i]
            else:
                bulk_mue[i] =  bulk_mue[i] - 0.01* bulk_mue[i]
                secondary_bulk_mue[i] =  secondary_bulk_mue[i] + 0.01* secondary_bulk_mue[i]
        '''
        
        
        
        print ("...with this differnece", warning, " I am proceeding with the computation..." )
        
        
        
        '''
        for i in range(len(bulk_mue)):
            if (bulk_mue[i]>secondary_bulk_mue[i]):
                bulk_mue[i] = secondary_bulk_mue[i]
            else: 
                secondary_bulk_mue[i]= bulk_mue[i]
       '''

    # Create output data with r positions, rho, and chemical potential values for each point
    output_data = []
    it  = 0
    for r_point in r_space_data:
        row = list(r_point)  # Initial columns are the position (x, y, z) in r-space
        
        if (it < int(0.5*len(r_space_data))):
            i=0
            for species_type, rho_value in species.items():
            
                row.append(rho_value)  # Append rho value for the species
                row.append(bulk_mue[i])  # Append the calculated chemical potential for the species    
                i = i+1
                
        else:
            i=0
            for species_type, rho_value in secondary_species.items():
            
                row.append(rho_value)  # Append rho value for the species
                row.append(secondary_bulk_mue[i])  # Append the calculated chemical potential for the species    
                i = i+1
            
        output_data.append(row)
        
        it  = it + 1 

    # Save output to text file
    np.savetxt(output_file, output_data, fmt="%.6f", header="x y z rho mue")
    
   

    print(f"\n\n... uniform rho and mue with a phase interface has been assigned to each r space points and exported to the supplied data file section ...\n\n\n")
    
# Run the function
#bulk_rho_mue_r_space()
