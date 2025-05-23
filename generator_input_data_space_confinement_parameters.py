# density functional minimizer/input data generator

# this is part of the code prints all the input data in the JSON format so that it could be read and executed by the different modules


# this code is also a gateway to block the code and stop it from running for a given configuration which is not developed yet ...





def data_exporter_space_confinement_parameters(input_file):

    import json

    
    
    # a submodule to check whether a value is convertible or not 
    
    def is_convertible_to_float(value):
        try:
            float(value)  # Try to convert the string to a float
            return True
        except ValueError:
            return False
            
            
            
            
            
    # species data is being fed in case if there is an interaction between particles and the aperiodic boundary conditions ...
    try:
        with open('input_data_particles_interactions_parameters.json', 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        print("... Error: The file 'input_data_particles_interactions.json' was not found please first run the code: 'generator_input_data_particles_interactions.py' ...")
        exit(0)
    
    prop = data["particles_interactions_parameters"]
    # Loop through all species and print their rho_frac values
    species_data={}
    for key, attributes in prop["species"].items():
        rho_frac={}
        rho_frac["rho_frac"] = attributes["rho_frac"]
        species_data[key] = rho_frac
    file.close()
    
    
    
    
    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Setting up the default box and space properties
    
    space_properties={}
    
    # confinement type specific position space simulation properties
    
    box_properties = {}
    sphere_properties={}
    cylinder_properties={}
    
    
    space_properties["dimension"]=1
    space_properties["confinement"]="pbox"
   
    box_properties["box_length"] = [5.0, 5.0, 5.0]
    box_properties["box_points"] = [1000, 1000, 1000]
    
    
    for line in lines:
        line = line.strip()  # Remove leading and trailing whitespace
        if not line or line.startswith("#"):
            continue  # Skip empty lines and comments

        # split the line into key and value
        if "=" in line and "interaction" not in line and "wall" not in line:
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()  # Remove quotes around strings

            if key == "space_dimension":
                space_properties["dimension"]= int(value)
                
        # setting up default confinement properties depending upon the dimension of the system
                
                space_properties["confinement"]="pbox"
                box_properties["box_length"]=[10.0, 0.0, 0.0]
                box_properties["box_points"]=[0.0, 0.0, 0.0]
                for i in range(int(value)):
                    box_properties["box_length"][i]=10.0
            
            elif key == "space_confinement":
                space_properties["confinement"] = value
                
                if (space_properties["confinement"] == "abox"):
                    #setting up the default properties for the aperiodic confinement
                    boundary_type = {}
                    boundary_type["aperiodicity_blocker"] = "NA"
                    
                
                
                
            elif key == "box_extension":
                box_properties["box_length"]=[10, 10, 10]
                box_properties["box_length"] = [float(x) for x in value.split(",")]
            elif key == "box_points":
                box_properties["box_points"]=[100, 100, 100]
                
                box_properties["box_points"]= [float(x) for x in value.split(",")]
                
            elif key ==  "sphere_extension":  
                sphere_properties["radius_of_spherical_confinement"]= [5.0]
                sphere_properties["radius_of_spherical_confinement"]= value
            elif key ==  "sphere_points":
                sphere_properties["sphere_points"]= [100]
                sphere_properties["sphere_points"]= value
                

         
        
        elif "aperiodicity_blocker" in line:
        
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            
            boundary_type["aperiodicity_blocker"] = value
            if boundary_type["aperiodicity_blocker"] == "wall" or   boundary_type["aperiodicity_blocker"] == "walls":
                walls_properties = {}
                walls_properties["walls_particles_type"] = ["a"]
                walls_properties["walls_position"] = [[0.0, 0.0, 0.0]]
                walls_properties["walls_normal"] = [[1.0, 0.0, 0.0 ]]
                walls_interactions_data = {}
                walls_interactions_data["primary"] = {}
                walls_interactions_data["secondary"] = {}
                walls_interactions_data["tertiary"] = {}
        
        
        elif "wall" in line and "aperiodicity_blocker" not in line  and  boundary_type["aperiodicity_blocker"] != "NA" and boundary_type["aperiodicity_blocker"] != "na":
            if space_properties["confinement"] == "abox":
                if boundary_type["aperiodicity_blocker"] == "wall" or   boundary_type["aperiodicity_blocker"] == "walls":
                
                    
                    dummy, wall_prop = line.split(":", 1)
                    
                    
                    key, value = wall_prop.split("=", 1)
                
                    key= key.strip()
                    value =value.strip()
                    
                    if (key=="particles"):
                        particles_name = value.split(",")
                        li=[]
                        for name in particles_name:
                            p_name = name.strip()
                            li.append(str(p_name))
                        walls_properties["walls_particles_type"] = li
                        
                        # setting up default values for the interactions of the fluid particles with the walls particles ...
                        
                        walls_flag={}
                       
                                                
                        for name1 in species_data:
                            for name2 in walls_properties["walls_particles_type"]:    
                                    interaction = name1.strip() + name2.strip()
                                    walls_flag[interaction]=[0,  0,  0]
                                    
                        
                      
                        
                    elif (key == "position"):
                        position_list = "".join(value.split())
                        position_list = position_list.split("),(")
                        li= []
                        for position in position_list:
                            position = position.strip("(")
                            position = position.strip(")")
                           
                            li.append([float (cord) for cord in position.split(",")])
                        walls_properties["walls_position"] = li 
                        
                        
                        
                    elif (key == "orientation" or key == "normal"):
                        orientation_list = "".join(value.split())
                        orientation_list = orientation_list.split("),(")
                        li= []
                        for orientation in orientation_list:
                            orientation = orientation.strip("(")
                            orientation = orientation.strip(")")
                            li.append([float (cord) for cord in orientation.split(",")])
                        walls_properties["walls_normal"] = li 
                        
                        
                        
                    elif "wall_interaction" in line:
                      
                       
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
                        
                        if pair in walls_flag:
                            
                            if walls_flag[pair][0]==0: 
                                
                                temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 2.0
                            }
                                walls_flag[pair][0]=1
                                
                                
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
                                
                                walls_interactions_data["primary"][pair]=temp
                            elif walls_flag[pair][1]==0:
                                
                                temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 2.0
                            }
                                walls_flag[pair][1]=1
                                
                                
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
                                
                                walls_interactions_data["secondary"][pair]=temp
                                
                                
                            elif walls_flag[pair][2]==0:
                                temp = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 2.0
                            }
                                walls_flag[pair][2]=1
                                for prop in properties_list:
                                
                                   
                                    param_key, param_value = prop.split("=", 1)
                                    param_key = param_key.sinput_data_space_confinement=trip()
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
                                
                                walls_interactions_data["tertiary"][pair]=temp
                                
                            else:
                            
                                print("Error: too many parameters have been defined for the interaction between the pair", pair, "which is not appropriate... please reduce and adjust the maximum number of interactions to 3 ")
                                exit(0)
                            
            # Save to interactions_data
    # Save the parsed data to JSON files
    
    
    
    
    input_data_space_confinement={}
    input_data_space_confinement["space_properties"] = space_properties
    
    
    
    
    if space_properties["confinement"]=="pbox":
        input_data_space_confinement["box_properties"] = box_properties
        with open('input_data_space_confinement_parameters.json', 'w') as f:
            json.dump({"space_confinement_parameters": input_data_space_confinement}, f , indent=4)
        f.close()     
        
        print("\n\n... space and confinement properties has been exported ... \n\n\n")
           
            
            
            
    elif space_properties["confinement"]=="abox":
        if boundary_type["aperiodicity_blocker"] == "wall" or boundary_type["aperiodicity_blocker"] == "walls":
            walls_properties["walls_interactions"] = walls_interactions_data 
            input_data_space_confinement["box_properties"] = box_properties
            input_data_space_confinement["walls_properties"] =  walls_properties
            with open('input_data_space_confinement_parameters.json', 'w') as f:
                json.dump({"space_confinement_parameters": input_data_space_confinement}, f , indent=4)
            f.close()
            print("\n\n... space and confinement properties has been exported ... \n\n\n")
            print("\n\n... periodicity has been removed with the introduction of walls ... the second simplest aperiodicity is being used in the system, the first is aperiodicity without any wall or confinement ... \n\n")
        elif (boundary_type["aperiodicity_blocker"] == "NA"  or boundary_type["aperiodicity_blocker"] == "na"):
            
            input_data_space_confinement["box_properties"] = box_properties
            with open('input_data_space_confinement_parameters.json', 'w') as f:
                json.dump({"space_confinement_parameters": input_data_space_confinement}, f , indent=4)
            f.close()
            print("\n\n... space and confinement properties has been exported ... \n\n\n")
            print("\n\n... periodicity has been removed without introducing any walls ... the simplest aperiodicity is being used in the system ... \n\n")
        else:
            print("\n\nerror: an unknown kind of aperiodicity has been introduced...\n\n")
            
            
            
        
    elif space_properties["confinement"]=="asphere":
        print("\n\n properties has not been specified for the asphere.. the code is still in development... Thank you so much...\n\n\n")
        exit(0)
        
        
        
        
    elif space_properties["confinement"]=="psphere":
        print("\n\n periodicity is not applicable in the spherical simulation confinement... please chose the physical options... Thank you so much...\n\n\n")
        exit(0)



    
    return 0
# Specify the input file path
#input_file_path = 'executor_input.in'  # Replace with your actual input file name
#dummy = data_exporter_space_confinement(input_file_path)




# this part of the code run the various segments of the code before running the actual minimizer codes.....



