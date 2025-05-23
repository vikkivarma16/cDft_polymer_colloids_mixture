# density functional minimizer/bulk rho and mue generator...

# this part of the code generates the bulk values of the density and the chemical potential in uniform way for each r space point supplied by the r-space files for all kind of the interactions specified in the interaction data json data profile... 
def calculator_coexistence_densities():
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
        total_rho = data_thermodynamic["simulation_thermodynamic_parameters"]["rho"]
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




    import numpy as np
    from scipy.optimize import root


    co_densities_a = [sp.symbols(f"xa_{i}") for i in range(len(sigmaij))]
    co_densities_b = [sp.symbols(f"xb_{i}") for i in range(len(sigmaij))]

    variables = []

    for i in range(len(co_densities_a)):
        variables.append(co_densities_a[i])
        variables.append(co_densities_b[i])



    def particle_numbers():
        pvec = []
        for species_type, rho_value in species.items():
            pvec.append(rho_value)
        return pvec






    def reduced_to_densities(reduced_vars):
        """
        Convert reduced variables to full species densities for a or b phase.

        Input: reduced_vars = [rhot, x1, x2, ..., xM] with M = N-1
        Output: rho = [rho_1, rho_2, ..., rho_N]
        """
        rhot = reduced_vars[0]
        fractions = reduced_vars[1:]

        # Build nested product chain:
        rho = []
        prod = 1.0
        
        for i in range(len(fractions)):
            prod = prod * (1-fractions[i])
            
        factor = 1.0
        for i in range(len(fractions)+1):
            val = rhot * prod *factor 
            rho.append(val)
            if (i<len(fractions)):
                prod *= 1/(1-fractions[i])
                factor = fractions[i]
            

        return rho  # List of N densities


    def coexistence_residual(vars):
        pvec = particle_numbers()  # should return [p1, p2, ..., pN]
        N = len(species)

        # Split input vars: assume order is [xa0, xa1, ..., xaM, xb0, xb1, ..., xbM, p]
        M = N - 1
        reduced_a = vars[0 : M + 1]         # rhot_a and its M variables
        reduced_b = vars[M + 1 : 2*M + 2]   # rhot_b and its M variables
        
        #print(vars)
        #print(reduced_a)
        #print(reduced_b)
        
        
        p = vars[-1]                        # phase fraction

        # Compute full densities
        rhoa = reduced_to_densities(reduced_a)
        rhob = reduced_to_densities(reduced_b)
        
        #print(rhoa)
        
        #print(rhob)
        
        #exit(0)

        # Compute chemical potentials and pressure
        mu_a = [mue_a_value[i](*rhoa) for i in range(N)]
        mu_b = [mue_b_value[i](*rhob) for i in range(N)]
        pressure_a = pressure_a_value(*rhoa)
        pressure_b = pressure_b_value(*rhob)

        # Coexistence equations
        eqs = []

        # µ_i^a = µ_i^b
        for i in range(N):
            eqs.append(mu_a[i] - mu_b[i])

        # Pressure equality
        eqs.append(pressure_a - pressure_b)

        # Mass conservation (weighted average = given total density)
        for i in range(N):
            eqs.append(p * rhoa[i] + (1 - p) * rhob[i] - pvec[i])

        return eqs


    def random_initial_guess():
        """
        Generates initial guess for:
        [rhot_a, x1_a, ..., xM_a, rhot_b, x1_b, ..., xM_b, p]
        """
        N = len(sigmaij)
        M = N - 1
        guess = []

        # Reduced variables for phase a
        guess.append(np.random.uniform(0.0, 1.0))  # rhot_a
        for i in range(M):
            if i == 0:
                guess.append(np.random.uniform(0.0, 0.4))  # for example: disjoint split
            else:
                guess.append(np.random.uniform(0.0, 1.0))

        # Reduced variables for phase b
        guess.append(np.random.uniform(0.0, 1.0))  # rhot_b
        for i in range(M):
            if i == 0:
                guess.append(np.random.uniform(0.6, 1.0))  # disjoint region
            else:
                guess.append(np.random.uniform(0.0, 1.0))

        # Phase fraction
        guess.append(np.random.uniform(0.0, 1.0))

        return guess



    def solve_general_coexistence(n_attempts=50):
        species_names = []
        for key, value in species.items():
            species_names.append(key)
        
        N = len(species_names)
        M = N - 1  # Number of reduced species

        for _ in range(n_attempts):
            guess = random_initial_guess()
            sol = root(coexistence_residual, guess, method='hybr')
            
            if sol.success:
                # Extract reduced variable blocks
                reduced_a = sol.x[0 : M + 1]          # [rhot_a, x1_a, ..., xM_a]
                reduced_b = sol.x[M + 1 : 2*M + 2]    # [rhot_b, x1_b, ..., xM_b]
                p = sol.x[-1]
                
                rhot_a, *x_a = reduced_a
                rhot_b, *x_b = reduced_b

                # General constraint check (matches random_initial_guess())
                def in_bounds(val, lo, hi):
                    return lo <= val <= hi

                valid = in_bounds(rhot_a, 0.0, 1.0) and in_bounds(rhot_b, 0.0, 1.0)
                valid = valid and in_bounds(x_a[0], 0.0, 0.4) and in_bounds(x_b[0], 0.6, 1.0)

                for i in range(1, M):
                    valid = valid and in_bounds(x_a[i], 0.0, 1.0) and in_bounds(x_b[i], 0.0, 1.0)

                valid = valid and in_bounds(p, 0.0, 1.0)

                if valid:
                    rhoa = reduced_to_densities(reduced_a)
                    rhob = reduced_to_densities(reduced_b)
                    pressure_a = pressure_a_value(*rhoa)
                    pressure_b = pressure_b_value(*rhob)
                    # Build result dictionary
                    result = {"species": {}}
                    for idx, name in enumerate(species_names):
                        result["species"][name] = {
                            "rho_frac": float(rhoa[idx]),
                            "secondary_rho_frac": float(rhob[idx])
                        }

                    

                    print("✅ Coexistence solution found and exported to the json file please open and check: and the pressure is given as", pressure_a, pressure_b)
                    print(f"Phase A: rhot = {rhot_a:.4f}, x = {x_a}")
                    print(f"Phase B: rhot = {rhot_b:.4f}, x = {x_b}")
                    print(f"p = {p:.4f}")
                    
                    import json
                    output_file = "Solution.json"
                    with open(output_file, "w") as f:
                        json.dump(result, f, indent=4)
                    return result

        print("❌ No valid coexistence solution found after multiple attempts.")
        
        rhoa = pvec
        rhob = pvec

        # Build result dictionary
        result = {"species": {}}
        for idx, name in enumerate(species_names):
            result["species"][name] = {
                "rho_frac": float(rhoa[idx]),
                "secondary_rho_frac": float(rhob[idx])
            }
        import json   
        output_file = "Solution.json"
        with open(output_file, "w") as f:
            json.dump(result, f, indent=4)
        
        return result




    result = solve_general_coexistence(n_attempts=100)


    
    return result


