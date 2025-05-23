# density functional minimizer/input data generator

# this is part of the code prints all the input data in the JSON format so that it could be read and executed by the different modules


# this code is also a gateway to block the code and stop it from running for a given configuration which is not developed yet ...

# this part is perfect and requires no changes ... 





def data_exporter_particles_interactions_parameters(input_file):

    import json
    # Initialize dictionaries to hold the data
    # number of dimension and confinement types
    
    
    
    # a submodule to check whether a value is convertible or not 
    
    def is_convertible_to_float(value):
        try:
            float(value)  # Try to convert the string to a float
            return True
        except ValueError:
            return False
    
    
    
    
    
    
    
    # type of particles in the system
    species_data = {}
    
    
    # interaction among the particles in the system
    interactions_data = {}
    interactions_data ["primary"] = {}
    interactions_data ["secondary"] = {}
    interactions_data ["tertiary"] = {}
    
    
    # thermodynamic parameters defined for the system
    thermodynamic_properties = {}

    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Parsing the input file
    # Default species data
    species = "a"
    species_data[species] = {"rho_frac": 1.0}
    
    
   
    
    
    
    
    for key in ["temperature", "rho", "iteration_max"]:
                thermodynamic_properties[key] = "NA"
    
    
    for line in lines:
        line = line.strip()  # Remove leading and trailing whitespace
        if not line or line.startswith("#"):
            continue  # Skip empty lines and comments

        # split the line into key and value
        if "=" in line and "interaction" not in line and "wall" not in line:
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()  # Remove quotes around strings
            if key == "species":
                species_names = value.split(",")
                
                rho_frac = 1.0 / len(species_names)
                for name in species_names:
                    species_data[name.strip()] = {"rho_frac": rho_frac}  # Initialize species with rho_frac
                
                # Initialize default interaction data for species
                i=0
                j=0
                flag={}
                for name1 in species_names:
                    j=0
                    for name2 in species_names:
                        if j>=i :
                            interaction = name1.strip() + name2.strip()
                            flag[interaction]=[0,  0,  0]
                            interactions_data["primary"][interaction] = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 0.0
                            }
                        j = j+1
                    i = i+1
                interactions_data["secondary"] = {}
                interactions_data["tertiary"] = {}
                
            elif key == "species_canonical_fraction":
                rho_values = value.split(",")
                for i, name in enumerate(species_data.keys()):
                    species_data[name]["rho_frac"] = float(rho_values[i].strip()) 
            elif key in ["temperature", "rho", "iteration_max"]:
                thermodynamic_properties[key] = float(value)

        # Handle interactions
        elif "interaction" in line and "wall" not in line:
            
            interaction, dummy = line.split("=", 1)
            
            
            dummy, properties_str = line.split(":", 1)
            
            dummy, pair= interaction.split(":", 1) 
            
            # Remove leading/trailing spaces from interaction type
            pair = pair.strip()

            # Step 2: Split the properties part (after the '=') by commas to get each key-value pair
            properties_list = properties_str.split(",")

            # Initialize a dictionary to store the property values
            properties_dict = {}

            # Step 3: Process each key-value pair and store them in a dictionary
            
            if pair in flag:
    
                if flag[pair][0]==0:   
                    flag[pair][0]=1
                    temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 0.0
                            }
                    
                    for prop in properties_list:
                    
                        param_key, param_value = prop.split("=", 1)
                        param_key = param_key.strip()
                        param_value = param_value.strip()
                        
                        # Only replace default values with input values
                    
                        if param_key != pair:
                            if is_convertible_to_float(param_value)==True:
                                temp[param_key] = float(param_value)
                                
                            else:
                                temp[param_key] = param_value
                                
                            if "type" == "hc":
                                temp["cutoff"]= temp["sigma"]
                                
                        elif param_key == pair :

                            temp["type"] = param_value
                    
                    interactions_data["primary"][pair]=temp
                    
                    
                elif flag[pair][1]==0:
                    
                    flag[pair][1]=1
                    temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 0.0
                            }
                    
                    for prop in properties_list:
                    
                        param_key, param_value = prop.split("=", 1)
                        param_key = param_key.strip()
                        param_value = param_value.strip()
                        
                        # Only replace default values with input values
                    
                        if param_key != pair:
                            if is_convertible_to_float(param_value)==True:
                                temp[param_key] = float(param_value)
                                
                            else:
                                temp[param_key] = param_value
                                
                            if "type" == "hc":
                                temp["cutoff"]= temp["sigma"]
                                
                        elif param_key == pair :

                            temp["type"] = param_value
                    
                    interactions_data["secondary"][pair]=temp
                    
                    
                elif flag[pair][2]==0:
                    
                    flag[pair][2]=1
                    temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 0.0
                            }
                    
                    for prop in properties_list:
                    
                       
                        param_key, param_value = prop.split("=", 1)
                        param_key = param_key.strip()
                        param_value = param_value.strip()
                        

                        # Only replace default values with input values
                    
                        
                        if param_key != pair:
                            if is_convertible_to_float(param_value)==True:
                                temp[param_key] = float(param_value)
                                
                            else:
                                temp[param_key] = param_value
                                
                            if "type" == "hc":
                                temp["cutoff"]= temp["sigma"]
                                
                        elif param_key == pair :

                            temp["type"] = param_value
                    
                    interactions_data["tertiary"][pair]=temp
                    
                else:
                
                    print("Error: too many parameters have been defined for the interaction between the pair", pair, "which is not appropriate... please reduce and adjust the maximum number of interactions to 3 ")
                    exit(0)
                    
                
        
       
            
            
            # Save to interactions_data

    # Save the parsed data to JSON files
    
    input_data_particles_interactions = {}
    input_data_particles_interactions["species"] = species_data
    input_data_particles_interactions["interactions"] = interactions_data

    with open('input_data_particles_interactions_parameters.json', 'w') as f:
        json.dump({"particles_interactions_parameters":input_data_particles_interactions}, f, indent=4)
        
    print("\n\n\n... particles properties have been exported ...\n\n\n")
    
    return 0
# Specify the input file path
#input_file_path = 'executor_input.in'  # Replace with your actual input file name
#dummy = data_exporter_particles_interactions_parameters(input_file_path)


input_file_name = 'executor_input.in'
data_exporter_particles_interactions_parameters(input_file_name)
# this part of the code run the various segments of the code before running the actual minimizer codes.....



