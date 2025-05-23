# density functional minimizer/r and k space points generator within the cartesian

# this part of the code generates the values of the r space points and the corresponding k space points within the cartesian framework 




import json 
import numpy as np 
def r_k_space():
    
    def read_input(filename):
        """Reads input file and returns parameters as a dictionary."""
        with open(filename, 'r') as f:
            params = json.load(f)
        return params

    def generate_r_space_1d(length_x, num_points_x):
        """Generates real space for 1D Cartesian grid."""
        x_space = np.linspace(0, length_x, int(num_points_x))
        return x_space, np.zeros_like(x_space), np.zeros_like(x_space)  # y and z are zeros

    def generate_r_space_2d(length_x, length_y, num_points_x, num_points_y):
        """Generates real space for 2D Cartesian grid."""
        x_space = np.linspace(0, length_x, int(num_points_x))
        y_space = np.linspace(0, length_y, int(num_points_y))
        x_grid, y_grid = np.meshgrid(x_space, y_space)
        
        # Assign z values to zero
        z_grid = np.zeros_like(x_grid)
        
        return x_grid, y_grid, z_grid

    def generate_r_space_3d(length_x, length_y, length_z, num_points_x, num_points_y, num_points_z):
        """Generates real space for 3D Cartesian grid."""
        x_space = np.linspace(0, length_x, int(num_points_x))
        y_space = np.linspace(0, length_y, int(num_points_y))
        z_space = np.linspace(0, length_z, int(num_points_z))
        x_grid, y_grid, z_grid = np.meshgrid(x_space, y_space, z_space)
        return x_grid, y_grid, z_grid

    def generate_k_space_1d(length_x, num_points_x):
        """Generates k-space for 1D Cartesian grid using fftfreq."""
        
        print("\n\n\n\n\nlength x", length_x, "\n\n\n\n\n")
        
        
        dkx = length_x / (int(num_points_x) - 1)
        kx = np.fft.fftfreq(int(num_points_x), d=dkx)
        return kx, np.zeros_like(kx), np.zeros_like(kx)  # ky and kz are zeros

    def generate_k_space_2d(length_x, length_y, num_points_x, num_points_y):
        """Generates k-space for 2D Cartesian grid using fftfreq."""
        dkx = length_x / (int(num_points_x) - 1)
        dky = length_y / (int(num_points_y) - 1)
        
        kx = np.fft.fftfreq(int(num_points_x), d=dkx)
        ky = np.fft.fftfreq(int(num_points_y), d=dky) 
        
        # Create grids for kx and ky and assign kz values to zero
        kx_grid, ky_grid = np.meshgrid(kx, ky)
        kz_grid = np.zeros_like(kx_grid)  # Create a zero array with the same shape as kx_grid

        return kx_grid, ky_grid, kz_grid

    def generate_k_space_3d(length_x, length_y, length_z, num_points_x, num_points_y, num_points_z):
        """Generates k-space for 3D Cartesian grid using fftfreq."""
        dkx = length_x / (int(num_points_x) - 1)
        dky = length_y / (int(num_points_y) - 1)
        dkz = length_z / (int(num_points_z) - 1)
        
        kx = np.fft.fftfreq(int(num_points_x), d=dkx) 
        ky = np.fft.fftfreq(int(num_points_y), d=dky) 
        kz = np.fft.fftfreq(int(num_points_z), d=dkz) 
        
        kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz)
        
        return kx_grid, ky_grid, kz_grid
    """Main function to read input and calculate r_space and k_space."""
    # Read input parameters from 'input_box_properties.json'
    
    
    
    params_box = read_input('input_box_properties.json')
    
    # Extract box properties
    box_length = params_box['box_properties']['box_length']
    num_points = params_box['box_properties']['box_points']

    # Read space properties from 'input_space_properties.json'
    params_space = read_input('input_space_properties.json')
    
    # Extract space properties
    dimension = int(params_space['space_properties']['dimension'])
    
    if dimension == 1:
        # Generate real space and k-space for 1D
        x, y, z = generate_r_space_1d(box_length[0], num_points[0])
        kx, ky, kz = generate_k_space_1d(box_length[0], num_points[0])

        # Save results to text files
        np.savetxt('r_space.txt', np.column_stack((x, y, z)))
        np.savetxt('k_space.txt', np.column_stack((kx, ky, kz)))
        print("\n\n\n... the corresponding 1D real space and k space data has been generated ...\n\n\n")

    elif dimension == 2:
        # Generate real space and k-space for 2D
        x, y, z = generate_r_space_2d(box_length[0], box_length[1], num_points[0], num_points[1])
        kx, ky, kz = generate_k_space_2d(box_length[0], box_length[1], num_points[0], num_points[1])

        # Save results to text files
        np.savetxt('r_space.txt', np.column_stack((x.flatten(), y.flatten(), z.flatten())))
        np.savetxt('k_space.txt', np.column_stack((kx.flatten(), ky.flatten(), kz.flatten())))
        print("\n\n\n... the corresponding 2D real space and k space data has been generated ...\n\n\n")

    elif dimension == 3:
        # Generate real space and k-space for 3D
        x, y, z = generate_r_space_3d(box_length[0], box_length[1], box_length[2], num_points[0], num_points[1], num_points[2])
        kx, ky, kz = generate_k_space_3d(box_length[0], box_length[1], box_length[2], num_points[0], num_points[1], num_points[2])

        # Save results to text files
        np.savetxt('r_space.txt', np.column_stack((x.flatten(), y.flatten(), z.flatten())))
        np.savetxt('k_space.txt', np.column_stack((kx.flatten(), ky.flatten(), kz.flatten())))
        print("\n\n\n... the corresponding 3D real space and k space data has been generated ...\n\n\n")

    else:
        print("Currently, only 1D, 2D, and 3D spaces are supported.")
        

    return 0

# Call the function here
#dummy = r_k_space()  # Call the main function

