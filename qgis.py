import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os

def knife_edge_diffraction_loss(d, h, wavelength):
    """
    Calculate the knife-edge diffraction loss.
    :param d: Distance between transmitter and receiver (meters)
    :param h: Height of the obstacle (meters)
    :param wavelength: Wavelength of the signal (meters)
    :return: Diffraction loss in dB
    """
    v = h * np.sqrt(2 / (wavelength * d))
    if v < -0.78:
        return 0
    else:
        return 6.9 + 20 * np.log10(np.sqrt((v - 0.1)**2 + 1) + v - 0.1)

def simulate_rf_coverage(dem_data, tx_location, tx_power_dbm, frequency_MHz, tx_height):
    # Constants
    c = 3e8  # Speed of light in m/s
    frequency_hz = frequency_MHz * 1e6  # Convert MHz to Hz
    wavelength = c / frequency_hz  # Wavelength in meters

    # Initialize coverage map
    rf_coverage = np.zeros_like(dem_data, dtype=float)

    # Transmitter location in DEM coordinates
    tx_row, tx_col = tx_location

    # Transmitter altitude from DEM
    tx_elevation = dem_data[tx_row, tx_col] + tx_height

    atmospheric_absorption_coefficient = 0.02  # Example coefficient in dB/m

    # Calculate RF coverage
    rows, cols = np.indices(dem_data.shape)
    distances = np.sqrt((rows - tx_row)**2 + (cols - tx_col)**2)
    altitudes = dem_data + tx_height

    path_loss_db = 20 * np.log10(distances) + 20 * np.log10(frequency_hz) - 147.55
    atmospheric_loss_db = atmospheric_absorption_coefficient * altitudes

    # Calculate diffraction loss using knife-edge model
    diffraction_loss_db = np.zeros_like(dem_data, dtype=float)
    
    for row in range(dem_data.shape[0]):
        for col in range(dem_data.shape[1]):
            if distances[row, col] > 0:
                los_path = np.linspace(tx_elevation, dem_data[row, col], num=int(distances[row, col]))
                max_obstacle_height = np.max(los_path)
                obstacle_height = max_obstacle_height - min(tx_elevation, dem_data[row, col])
                diffraction_loss_db[row, col] = knife_edge_diffraction_loss(distances[row, col], obstacle_height, wavelength)

    total_loss_db = path_loss_db + atmospheric_loss_db + diffraction_loss_db
    rf_coverage = tx_power_dbm - total_loss_db

    return rf_coverage

def latlon_to_rowcol(dem_path, lat, lon):
    # Open the DEM file
    with rasterio.open(dem_path) as dem_dataset:
        # Get the transform and metadata
        transform = dem_dataset.transform
        crs = dem_dataset.crs
        
        # Check if the DEM is in WGS84, if not, reproject the coordinates
        if crs.to_string() != 'EPSG:4326':
            from pyproj import Transformer
            transformer = Transformer.from_crs('EPSG:4326', crs.to_string(), always_xy=True)
            lon, lat = transformer.transform(lon, lat)
        
        # Convert latitude and longitude to row and column
        row, col = ~transform * (lon, lat)
        
        # Convert to integer row and column indices
        row = int(row)
        col = int(col)
        
        return row, col

# Open the TIF file
dem_path = 'file.tif'

with rasterio.open(dem_path) as src:
    # Read the metadata
    meta = src.meta

    # Extract the elevation data from the DEM
    dem = src.read(1)
    dem_transform = src.transform
    dem_crs = src.crs
    print('Coord Ref System:', dem_crs.to_string())
    
    # Example coordinates
    latitude = 37.9899576
    longitude = -78.4940777

    row, col = latlon_to_rowcol(dem_path, latitude, longitude)
        
    # Transmitter parameters
    tx_location = (row, col)  # Example location (row, col)
    tx_power = 1  # Transmitter power in dBm
    tx_height = 1  # Transmitter height in meters above ground from DEM
    frequency_MHz = 97.5  # Frequency in MHz

    # Simulate coverage
    coverage = simulate_rf_coverage(dem, tx_location, tx_power, frequency_MHz, tx_height)

    # Visualize the coverage map and save to file
    plt.figure(figsize=(10, 8))
    plt.imshow(coverage, cmap='viridis', extent=(0, dem.shape[1], 0, dem.shape[0]))
    plt.colorbar(label='Signal Strength (dBm)')
    plt.title('Simulated RF Coverage Map')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')

    # Save the plot to a file instead of showing it interactively
    plt.savefig('simulated_rf_coverage_map.png')
    print("The simulated RF coverage map has been saved to 'simulated_rf_coverage_map.png'.")

    # Write coverage to disk
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

    # Perform analysis on the elevation data
    print("Minimum elevation:", np.min(dem))
    print("Maximum elevation:", np.max(dem))
