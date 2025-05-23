# density functional minimizer/potential generator in r and k space

# this code generates potential values for the distance measured from the center of a particle for the visualization and nature of the potential purpose both in the r and k space within the cartesian framework....
# in K space its a repetitive values from all the particles so its values are for the set of frequencies from zero to maximum relevant values... 

import numpy as np
def pair_potential_custom_1(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array

    # Calculate the WCA potential based on the distance
    v_r[r < r_cutoff] = -epsilon  # Region 1: r < 2^(1/6) * sigma
    v_r[(r >= r_cutoff) & (r < r_upper_bound)] = 4 * epsilon * ((sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**12 - 2 * (sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**6)  # Region 2: 2^(1/6) * sigma < r < 5*sigma
    v_r[r >= r_upper_bound] = 0  # Region 3: r > 5 * sigma
    return 0
    
    
    
def pair_potential_integrant_custom_1(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array

    # Calculate the WCA potential based on the distance
    
    if (r<r_cutoff):
        v_r = -epsilon
    elif r >= r_cutoff and r < r_upper_bound:
        v_r = 4 * epsilon * ((sigma / r)**12 - 2 * (sigma / r)**6)
    elif r >= r_upper_bound:
        v_r = 0
    
    return 0
    
    
    
def pair_potential_custom_2(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array
    
    # Calculate the WCA potential based on the distance
    v_r[r < r_cutoff] = -epsilon  # Region 1: r < 2^(1/6) * sigma
    v_r[(r >= r_cutoff) & (r < r_upper_bound)] = 4 * epsilon * ((sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**12 - 2 * (sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**6)  # Region 2: 2^(1/6) * sigma < r < 5*sigma
    v_r[r >= r_upper_bound] = 0  # Region 3: r > 5 * sigma
    
    return 0
    
    
    
'''
import calculator_pair_potential_custom
for i in range(1, 4):  # Adjust the range based on the number of functions
    func_name = f'pair_potential_custom_{i}'
    func = getattr(calculator_pair_potential_custom, func_name)
    result = func()  # Call the function
    print(f"Result from {func_name}: {result}")
'''

def pair_potential_integrant_custom_2(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array

    # Calculate the WCA potential based on the distance
    
    print("... here the problem persist ...\n\n")
    if (r<r_cutoff):
        v_r = -epsilon
    elif r >= r_cutoff and r < r_upper_bound:
        v_r = 4 * epsilon * ((sigma / r)**12 - 2 * (sigma / r)**6)
    elif r >= r_upper_bound:
        v_r = 0
    
    return 0
    
    

def pair_potential_custom_3(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array
    
    # Calculate the WCA potential based on the distance
    v_r[r <= EPSILON] = 2000000000
    v_r[(r < r_upper_bound) & (r > EPSILON)] = epsilon * ((2.0/5.0)*(sigma / r[(r < r_upper_bound) & (r > EPSILON)])**10 -   (sigma / r[(r < r_upper_bound) & (r > EPSILON)])**4)  # Region 2: 2^(1/6) * sigma < r < 5*sigma
    v_r[r >= r_upper_bound] = 0  # Region 3: r > 5 * sigma
    
    return 0
    
    
    
def pair_potential_integrant_custom_3(r, v_r, epsilon, sigma):

    EPSILON = 0.0001
    
    import numpy as np
    
    
    r_cutoff = 2**(1/6) * sigma
    r_upper_bound = 5 * sigma

    # Initialize the potential array
    
    # Calculate the WCA potential based on the distance
    if (r < r_upper_bound) and (r>EPSILON):
        v_r = epsilon * ((2.0/5.0)*(sigma / r)**10 -   (sigma / r)**4)  # Region 2: 2^(1/6) * sigma < r < 5*sigma
    elif (r >= r_upper_bound):
        v_r = 0  # Region 3: r > 5 * sigma
    if (r<EPSILON):
        v_r = 20000000
    
    return 0    

'''
import calculator_pair_potential_custom
for i in range(1, 4):  # Adjust the range based on the number of functions
    func_name = f'pair_potential_custom_{i}'
    func = getattr(calculator_pair_potential_custom, func_name)
    result = func()  # Call the function
    print(f"Result from {func_name}: {result}")
'''

