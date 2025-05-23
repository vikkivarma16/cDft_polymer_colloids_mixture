# density functional minimizer/input data generator

# this code no need any interaction once edited in its full functional formate...




def data_exporter_simulation_thermodynamic_parameters(input_file, result):

    import json
   
    def is_convertible_to_float(value):
        try:
            float(value)  # Try to convert the string to a float
            return True
        except ValueError:
            return False
    

    # thermodynamic parameters defined for the system
    
    thermodynamic_properties = {}

    # setting up default thermodynamic properties ... 
    
    
    with open(input_file, 'r') as file:
        lines = file.readlines()
    # Parsing the input file

    for key in ["temperature", "rho", "iteration_max", "pressure"]:
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

            
            if key in ["temperature", "rho", "iteration_max", "pressure"]:
                
                if key== "iteration_max":
                    thermodynamic_properties[key] = int(value)
                else :
                    thermodynamic_properties[key] = float(value)
    
    thermodynamic_properties["species"] = result["species"] 

    if thermodynamic_properties["temperature"]== "NA" or thermodynamic_properties["rho"]== "NA" or thermodynamic_properties["iteration_max"]== "NA":
        print ("error: specify thermodynamic properties like rho=value, temperature= value, iteration_max= value, in row wise stacking format otherwise simulation will not run..")
        exit(0)
    else:
        with open('input_data_simulation_thermodynamic_parameters.json', 'w') as f:
            json.dump({"simulation_thermodynamic_parameters": thermodynamic_properties}, f, indent=4)
            
        print("\n\n\n... thermodynamic parameters have been exported ...\n\n\n")
    
    return 0
# Specify the input file path
#input_file_path = 'executor_input.in'  # Replace with your actual input file name
#dummy = data_exporter_simulation_thermodynamic_parameters(input_file_path)


# this part of the code run the various segments of the code before running the actual minimizer codes.....



