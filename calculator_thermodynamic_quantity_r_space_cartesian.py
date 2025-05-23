# density functional minimizer/thermodynamic quantities generator...


# this part of the code calculates the thermodynamic quantities and generates the variables with the reduced dimension from the MKS unit...



def thermodynamic_variables():

    import numpy as np
    from scipy.optimize import minimize

    # Constants and Parameters
    T = 300  # Temperature in K
    Pressure = 1   # Pressure in bar
    Sigma_ff = 3.305  # Lennard-Jones sigma [Å]
    distance_unit=Sigma_ff
    Eps_kb = 118.05   # Lennard-Jones epsilon/k_B [K]
    kb = 1.3806488e-23  # Boltzmann constant [J/K]

    # Non-dimensionalized Temperature and Pressure
    T_star = T / Eps_kb  # Reduced temperature



    Eps_kb = Eps_kb / Eps_kb

    # Hard-sphere diameter from Lennard-Jones sigma
    d_HS = 3.384 / Sigma_ff  # Diameter in reduced units
      # Cut-off for WCA potential

    Sigma_ff=Sigma_ff/Sigma_ff

    rc_MF = 5.0 * Sigma_ff


    # Weeks-Chandler-Andersen (WCA) potential mean-field contribution
    def wca_potential(r, eps, sigma, rc):
        if r > rc:
            return 0
        rmin = 2**(1/6) * sigma
        if r >= rmin:
            return 4 * eps * ((sigma / r)**12 - (sigma / r)**6) - 4 * eps * ((sigma / rc)**12 - (sigma / rc)**6)
        else:
            return 4 * eps * (1 - (r / rmin)**2)

    # Carnahan-Starling Equation of State for Hard Spheres
    def pressure_cs(rho, T_star, Sigma_ff):
        eta = np.pi * rho * Sigma_ff**3 / 6.0
        return T_star * rho * (1 + eta + eta**2 - eta**3) / (1 - eta)**3

    # Mean-field WCA potential contribution to pressure
    def pressure_wca(rho, Sigma_ff, Eps_kb, rc_MF):
        value=(-np.sqrt(2) * (32.0 / 9.0) * np.pi * Eps_kb * Sigma_ff**3 +
                               (16.0 / 3.0) * np.pi * Eps_kb * Sigma_ff * ((Sigma_ff / (5.0*rc_MF))**3 - 
                               (1.0 / 3.0) * (Sigma_ff / (5.0*rc_MF))**9))
                               
     
        return 0.5 * rho**2 * value

    # Objective function to minimize the difference between calculated and target pressure
    def objective_function(rho, P_Star):
        P_cs = pressure_cs(rho, T_star, Sigma_ff)
        P_wca = pressure_wca(rho, Sigma_ff, Eps_kb, rc_MF)
        P_total = P_cs + P_wca
        
        if (rho>4.0) :
            print(P_Star, P_cs, rho)
        
       
        return (P_total - P_star)**2  # Objective: minimize the difference between calculated and target pressure


        
    pressure_range = np.linspace(0.01, 10.0, 1000)  # 100 points between 0.1 atm and 10 atm

    # Open a text file to write the results
    with open('eos_rho_vs_pressure.txt', 'w') as file:
        # Write the header
       
        
        # Loop over pressure range and calculate corresponding density
        initial_guess = 1.e-5  # Initial guess for density
            
        previous_row= initial_guess
        for P in pressure_range:
            P_star = (P * 1.e5) / (kb * Eps_kb) * (distance_unit * 1.e-10)**3  # Reduced pressure in Lennard-Jones units
            
            # Minimize the objective function to find the optimal bulk density (rho_H)
            #result = minimize(objective_function, initial_guess, args=(P_star,), method='BFGS')
            
            bounds = [(1e-10, 2)]
            
            result = minimize(objective_function, previous_row, args=(P_star,), method='L-BFGS-B', bounds=bounds)
            rho_H = result.x[0]
            
            previous_row=rho_H
            
            eta = np.pi * rho_H * Sigma_ff**3 / 6.0
            
            # Write the pressure and corresponding bulk density to the file
            file.write(f"{eta:>30.12f}  {P:>20.6f}\n")

    print("EOS data written to 'eos_rho_vs_pressure.txt'")


    P_star = (Pressure * 1.e5) / (kb * Eps_kb) * (distance_unit * 1.e-10)**3  # Reduced pressure in Lennard-Jones units
            
    # Minimize the objective function to find the optimal bulk density (rho_H)
    #result = minimize(objective_function, initial_guess, args=(P_star,), method='BFGS')

    bounds = [(1e-10, 2)]

    result = minimize(objective_function, previous_row, args=(P_star,), method='L-BFGS-B', bounds=bounds)
    rho_H = result.x[0]




    # Ideal-Gas contribution
    mu_ID = T * np.log(rho_H)

    # Hard-Sphere Contribution
    EtaH = np.pi * rho_H * Sigma_ff**3 / 6.0
    mu_HS = T * (EtaH * (8.0 - 9.0 * EtaH + 3.0 * EtaH**2) / (1.0 - EtaH)**3)

    # Dispersive contribution (WCA)
    mu_MF_WCA = rho_H * (-np.sqrt(2) * (32.0 / 9.0) * np.pi * Eps_kb * Sigma_ff**3 + 
                         (16.0 / 3.0) * np.pi * Eps_kb * Sigma_ff**3 * ((Sigma_ff / (5.0 * Sigma_ff))**3 - 
                         (1.0 / 3.0) * (Sigma_ff / (5.0 * Sigma_ff))**9))

    # Total chemical potential
    mu_H = mu_ID + mu_HS + mu_MF_WCA

    # Output the result
    print(f"Ideal Gas Contribution (mu_ID): {mu_ID}")
    print(f"Hard-Sphere Contribution (mu_HS): {mu_HS}")
    print(f"Dispersive Contribution (mu_MF_WCA): {mu_MF_WCA}")
    print(f"Total Chemical Potential (mu_H): {mu_H}")

    return 0; 
thermodynamic_variables()
