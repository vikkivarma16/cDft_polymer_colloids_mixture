# density functional minimizer/executor

# this is the main code which read and write the data in an executable format and then also run the program for the calculation...






# this part read the input file and print the json file for the further processing of the input file in the selected format...







# Density Functional Minimizer/Executor
# Main code that reads, writes, and processes input data to execute the program

# Density Functional Minimizer/Executor
import numpy as np
import json
from generator_k_and_r_space_box import r_k_space as rks
from generator_input_data import data_exporter
from generator_uniform_rho_mue_r_space import bulk_rho_mue_generator as rmg
from generator_pair_potential import potential_r_k_space as potentials


# File paths
input_file_name = 'executor_input.in'
json_file_space_properties = "input_space_properties.json"
json_file_species = "input_species_properties.json"
json_file_interaction = "input_interactions_properties.json"
json_file_thermodynamics = "input_thermodynamic_properties.json"

# Read and export configuration data
try:
    data_exporter(input_file_name)
    print("Input data exported successfully.")
except Exception as e:
    print("Error exporting input data:", e)

# Load space properties
try:
    with open(json_file_space_properties, "r") as file:
        space_properties = json.load(file)["space_properties"]
    space_dimension = space_properties["dimension"]
    space_confinement = space_properties["confinement"]
except FileNotFoundError:
    print("Space properties file not found.")
    exit()

# Load thermodynamic properties
try:
    with open(json_file_thermodynamics, "r") as file:
        thermodynamics = json.load(file)["thermodynamic_properties"]
    temperature, rho, iteration = thermodynamics["temperature"], thermodynamics["rho"], thermodynamics["iteration_max"]
    print(f"\n\n... thermodynamic properties loaded: Temperature = {temperature}, Rho = {rho}, Max Iterations = {iteration} ...\n\n")
except FileNotFoundError:
    print("Thermodynamic properties file not found.")
    exit()

# Load species properties
try:
    with open(json_file_species, 'r') as file:
        species = {k: v["rho_frac"] for k, v in json.load(file)["species"].items()}
    print(f"\n... loaded {len(species)} species for simulation ...\n\n")
except FileNotFoundError:
    print("Species properties file not found.")
    exit()

# Load interaction properties
try:
    with open(json_file_interaction, 'r') as file:
        data_interactions = json.load(file)["interactions"]
        interaction_types = {}
        closest_distances = {}
        interaction_strength = {}
        cutoff_ranges = {}
        
        interaction_types["primary"] = {k: v["type"] for k, v in data_interactions["primary"].items()}
        closest_distances["primary"] = {k: v["sigma"] for k, v in data_interactions["primary"].items()}
        interaction_strength["primary"] = {k: v["epsilon"] for k, v in data_interactions["primary"].items()}
        cutoff_ranges["primary"] = {k: v["cutoff"] for k, v in data_interactions["primary"].items()}
        
        interaction_types["secondary"] = {k: v["type"] for k, v in data_interactions["secondary"].items()}
        closest_distances["secondary"] = {k: v["sigma"] for k, v in data_interactions["secondary"].items()}
        interaction_strength["secondary"] = {k: v["epsilon"] for k, v in data_interactions["secondary"].items()}
        cutoff_ranges["secondary"] = {k: v["cutoff"] for k, v in data_interactions["secondary"].items()}
        
        interaction_types["tertiary"] = {k: v["type"] for k, v in data_interactions["tertiary"].items()}
        closest_distances["tertiary"] = {k: v["sigma"] for k, v in data_interactions["tertiary"].items()}
        interaction_strength["tertiary"] = {k: v["epsilon"] for k, v in data_interactions["tertiary"].items()}
        cutoff_ranges["tertiary"] = {k: v["cutoff"] for k, v in data_interactions["tertiary"].items()}
        
        
    print("\n... interactions loaded successfully and are given as ...")
except FileNotFoundError:
    print("Interaction properties file not found.")
    exit()


print("\n") 
interaction_level = "primary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... primary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")          
interaction_level = "secondary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... secondary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")       
interaction_level = "tertiary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... tertiary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}  \n")
          
print("\n\n\n")

# Check for confinement and execute dependent functions
if space_confinement == "pbox":
    print(f"...preparing for the simulation in {space_dimension}D space with periodic boundary conditions...\n\n\n")
    try:
        
        
        
        # preliminary calculation before minimization....
        print("...preliminary calculation are being done...\n\n")
        rks()  # Generate r and k space...
        potentials()  # Visualize potential functions...
        # if there is a hard core interaction then the FMT weight will be required...
        if (space_dimension==1):
            from calculator_FMT_weights_1d_cartesian import fmt_weights_1d as cf1
            cf1()
            
        elif (space_dimension==2):
            # this part is waiting for the editing purpose...
            from calculator_FMT_weights_2d_cartesian import fmt_weights_1d as cf2
            cf2()
        
        
            
       
        rmg()  # Generate rho and chemical potential data
        
        
    except Exception as e:
        print("Error in generating r/k space, plotting potentials, or rho-mue data:", e)
else:
    print("Selected confinement not supported.")
    exit()

