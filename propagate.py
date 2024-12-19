import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from numba import njit, prange
from pyproj import Transformer

@njit
def knife_edge_diffraction_loss(d, h, wavelength):
    v = h * np.sqrt(2 / (wavelength * d))
    if v < -0.78:
        return 0
    else:
        return 6.9 + 20 * np.log10(np.sqrt((v - 0.1)**2 + 1) + v - 0.1)

def antenna_gain_pattern(antenna_type, theta, phi, hpbw, fnbw):
    if antenna_type == 'dipole':
        return 2.15
    elif antenna_type == 'directional':
        return 8 - 0.1 * theta
    elif antenna_type == 'parabolic':
        return 20 - 0.2 * theta
    else:
        return 0

def simulate_rf_coverage(dem_data, tx_location, tx_power_dbm, frequency_MHz, tx_height, antenna_type, angle_of_attack, direction, hpbw, fnbw):
    """
    Simulates the RF coverage over a terrain using a digital elevation model (DEM).

    Parameters:
    dem_data (numpy.ndarray): The digital elevation model data.
    tx_location (tuple): The (row, col) location of the transmitter in the DEM.
    tx_power_dbm (float): The transmitter power in dBm.
    frequency_MHz (float): The frequency of the transmission in MHz.
    tx_height (float): The height of the transmitter above the ground in meters.
    antenna_type (str): The type of antenna ('dipole', 'directional', 'parabolic').
    angle_of_attack (float): The angle of attack of the antenna in degrees.
    direction (float): The direction the antenna is facing in degrees.
    hpbw (float): The half-power beamwidth of the antenna in degrees.
    fnbw (float): The first null beamwidth of the antenna in degrees.

    Returns:
    numpy.ndarray: The simulated RF coverage map.
    """
    
    # Speed of light in meters per second, fundamental constant for RF calculations
    c = 3e8
    
    # Convert frequency from MHz to Hz for accurate wavelength calculation
    frequency_hz = frequency_MHz * 1e6
    
    # Calculate the wavelength in meters, essential for understanding wave propagation characteristics
    wavelength = c / frequency_hz

    # Extract the row and column indices of the transmitter location from the DEM
    tx_row, tx_col = tx_location
    
    # Calculate the elevation of the transmitter above sea level by adding the transmitter height to the DEM elevation
    tx_elevation = dem_data[tx_row, tx_col] + tx_height

    # Create a grid of row and column indices for the DEM, representing the entire terrain
    rows, cols = np.indices(dem_data.shape)
    
    # Calculate the distance from the transmitter to each cell in the DEM using Euclidean distance formula
    distances = np.sqrt((rows - tx_row)**2 + (cols - tx_col)**2)
    
    # set the signal sample height about 6 feet from the ground
    altitudes = dem_data + 2
    
    # Calculate the free-space path loss in dB using the Friis transmission equation
    # This accounts for the loss of signal strength over distance in free space
    path_loss_db = 20 * np.log10(distances + 1e-6) + 20 * np.log10(frequency_hz) - 147.55
    
    # Calculate the atmospheric loss in dB, which increases with altitude due to factors like air density and humidity
    atmospheric_loss_db = 0.02 * altitudes

    # Calculate the diffraction loss using terrain data, which accounts for signal bending around obstacles
    diffraction_loss_db = calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength)
    
    # Calculate the elevation angles from the transmitter to each cell in degrees
    # These angles are crucial for understanding the vertical spread of the signal
    elevation_angles = np.degrees(np.arctan2(rows - tx_row, distances))
    
    # Calculate the azimuth angles from the transmitter to each cell in degrees
    # These angles help in understanding the horizontal spread and directionality of the signal
    azimuth_angles = (np.degrees(np.arctan2(cols - tx_col, rows - tx_row)) - direction) % 360
    
    # Calculate the antenna gain for each cell based on its elevation and azimuth angles
    # Antenna gain patterns vary with angle and type, affecting signal strength in different directions
    antenna_gain_db = np.vectorize(antenna_gain_pattern)(antenna_type, elevation_angles - angle_of_attack, azimuth_angles, hpbw, fnbw)

    # Sum all losses (path, atmospheric, diffraction) and subtract antenna gain to get total loss in dB
    total_loss_db = path_loss_db + atmospheric_loss_db + diffraction_loss_db - antenna_gain_db
    
    # Calculate RF coverage by subtracting total loss from transmitter power
    # This gives the signal strength at each point in the DEM
    rf_coverage = tx_power_dbm - total_loss_db

    return rf_coverage

@njit(parallel=True)
def calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength):
    # Get the number of rows and columns in the DEM data
    rows, cols = dem_data.shape
    
    # Initialize an array to store the diffraction loss values
    diffraction_loss_db = np.zeros((rows, cols), dtype=np.float32)

    # Loop over each row in parallel
    for row in prange(rows):
        # Loop over each column in parallel
        for col in prange(cols):
            # Calculate the distance from the transmitter to the current cell
            distance = np.sqrt((row - tx_row)**2 + (col - tx_col)**2)
            
            # If the distance is greater than 0, calculate the diffraction loss
            if distance > 0:
                # Number of samples along the path for interpolation
                num_samples = int(distance) * 10
                
                # Interpolate the heights along the path from the transmitter to the current cell
                path_heights = np.linspace(tx_elevation, dem_data[row, col], num=num_samples)
                
                # Interpolate the row and column indices along the path
                row_indices = np.linspace(tx_row, row, num=num_samples).astype(np.int32)
                col_indices = np.linspace(tx_col, col, num=num_samples).astype(np.int32)
                
                # Initialize the maximum obstacle height to negative infinity (for dBm value)
                max_obstacle_height = -np.inf
                
                # Loop over each sample along the path from transmitter to cell
                for i in range(num_samples):
                    r_idx = row_indices[i]
                    c_idx = col_indices[i]
                    
                    # Update the maximum obstacle height if the current cell is higher
                    if dem_data[r_idx, c_idx] > max_obstacle_height:
                        max_obstacle_height = dem_data[r_idx, c_idx]
                
                # Calculate the obstacle height relative to the path heights
                obstacle_height = max_obstacle_height - np.min(path_heights)
                
                # Calculate the diffraction loss using the knife-edge model
                diffraction_loss_db[row, col] = knife_edge_diffraction_loss(distance, obstacle_height, wavelength)
    
    # Return the array of diffraction loss values
    return diffraction_loss_db

def latlon_to_rowcol(dem_dataset, lat, lon):
    transform = dem_dataset.transform
    crs = dem_dataset.crs
    
    if crs.to_string() != 'EPSG:4326':
        transformer = Transformer.from_crs('EPSG:4326', crs.to_string(), always_xy=True)
        lon, lat = transformer.transform(lon, lat)
    
    col, row = ~transform * (lon, lat)
    row = int(row)
    col = int(col)
    
    return row, col

def main(dem_path):

    try:
        with rasterio.open(dem_path) as src:
            meta = src.meta
            dem = src.read(1)
            dem_transform = src.transform
            dem_crs = src.crs
            print('Coord Ref System:', dem_crs.to_string())
            
            #38.0512799,-78.5483868

            latitude = 38.0512799                    
            longitude = -78.5483868

            row, col = latlon_to_rowcol(src, latitude, longitude)
                
            tx_location = (row, col)
            tx_power_dbm = 1.0
            tx_height = 10
            frequency_MHz = 700
            antenna_type = 'parabolic'
            angle_of_attack = 0.0
            direction = 315
            hpbw = 60
            fnbw = 120

            coverage = simulate_rf_coverage(dem, tx_location, tx_power_dbm, frequency_MHz, tx_height, antenna_type, angle_of_attack, direction, hpbw, fnbw)

            plt.figure(figsize=(10, 8))
            plt.imshow(coverage, cmap='viridis', extent=(0, dem.shape[1], 0, dem.shape[0]))
            plt.colorbar(label='Signal Strength (dBm)')
            plt.title('Simulated RF Coverage Map')
            plt.xlabel('Column Index')
            plt.ylabel('Row Index')
            plt.savefig('simulated_rf_coverage_map.png')            

            output_path = 'coverage_map.tif'
            if os.path.exists(output_path):
                os.remove(output_path)
            with rasterio.open(
                output_path,
                'w',
                driver='GTiff',
                height=coverage.shape[0],
                width=coverage.shape[1],
                count=1,
                dtype=coverage.dtype,
                crs=dem_crs,
                transform=dem_transform,
            ) as dst:
                dst.write(coverage, 1)

            print("Done")
    except Exception as e:
        print(f"An error occurred: {e}")

# Generate a random terrain (DEM data) using the diamond-square algorithm for testing
def diamond_square(size=129, scale=100):
    def diamond_step(arr, step_size, scale):
        half_step = step_size // 2
        for y in range(half_step, size-1, step_size):
            for x in range(half_step, size-1, step_size):
                avg = (arr[y-half_step][x-half_step] +
                       arr[y-half_step][x+half_step] +
                       arr[y+half_step][x-half_step] +
                       arr[y+half_step][x+half_step]) / 4.0
                arr[y][x] = avg + (np.random.rand() * 2 * scale) - scale

    def square_step(arr, step_size, scale):
        half_step = step_size // 2
        for y in range(0, size-1+half_step, half_step):
            for x in range((y + half_step) % step_size, size-1+half_step, step_size):
                avg = (arr[(y-half_step) % (size-1)][x] +
                       arr[(y+half_step) % (size-1)][x] +
                       arr[y][(x-half_step) % (size-1)] +
                       arr[y][(x+half_step) % (size-1)]) / 4.0
                arr[y][x] = avg + (np.random.rand() * 2 * scale) - scale

    arr = np.zeros((size, size))
    arr[0][0] = np.random.rand() * scale
    arr[0][size-1] = np.random.rand() * scale
    arr[size-1][0] = np.random.rand() * scale
    arr[size-1][size-1] = np.random.rand() * scale

    step_size = size - 1
    while step_size > 1:
        diamond_step(arr, step_size, scale)
        square_step(arr, step_size, scale)
        step_size //= 2
        scale /= 2.0

    return arr

# Test the simulation function with the generated DEM data
def test_simulate_rf_coverage():
    #dem_data = generate_random_dem()
    dem_data = diamond_square()
    
    # Define test parameters
    tx_location = (50, 50) # center of the DEM
    tx_power_dbm = 1.0
    frequency_MHz = 700
    tx_height = 10
    antenna_type = 'parabolic'
    angle_of_attack = 0.0
    direction = 315
    hpbw = 60
    fnbw = 120
    
    # Run the simulation function with the test parameters and DEM data
    coverage_map = simulate_rf_coverage(dem_data, tx_location, tx_power_dbm,
                                        frequency_MHz, tx_height,
                                        antenna_type,
                                        angle_of_attack,
                                        direction,
                                        hpbw,
                                        fnbw)
    
    # Plot the coverage map for visual inspection
    plt.figure(figsize=(10, 8))
    plt.imshow(coverage_map, cmap='viridis', extent=(0, dem_data.shape[1], 0, dem_data.shape[0]))
    plt.colorbar(label='Signal Strength (dBm)')
    plt.title('Simulated RF Coverage Map (Test)')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.savefig('simulated_rf_coverage_map.png')
    

test_simulate_rf_coverage()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate RF coverage using a DEM file or generate a random terrain.")
    parser.add_argument('filename', nargs='?', default=None, help="Path to the DEM file (must end with .tif). If not provided, a random terrain will be generated.")
    args = parser.parse_args()

    if args.filename and args.filename.endswith('.tif'):
        main(args.filename)
    else:
        test_simulate_rf_coverage()
