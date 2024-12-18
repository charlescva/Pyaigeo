import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os
from numba import njit, prange

@njit
def knife_edge_diffraction_loss(d, h, wavelength):
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

    # Transmitter location in DEM coordinates
    tx_row, tx_col = tx_location
    tx_elevation = dem_data[tx_row, tx_col] + tx_height

    rows, cols = np.indices(dem_data.shape)
    distances = np.sqrt((rows - tx_row)**2 + (cols - tx_col)**2)
    altitudes = dem_data + tx_height

    path_loss_db = 20 * np.log10(distances + 1e-6) + 20 * np.log10(frequency_hz) - 147.55
    atmospheric_loss_db = 0.02 * altitudes  # Assuming absorption coefficient of 0.02 dB/m

    diffraction_loss_db = calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength)
    
    total_loss_db = path_loss_db + atmospheric_loss_db + diffraction_loss_db
    rf_coverage = tx_power_dbm - total_loss_db

    return rf_coverage

@njit(parallel=True)
def calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength):
    diffraction_loss_db = np.zeros_like(dem_data, dtype=np.float32)
    rows, cols = dem_data.shape
    
    for row in prange(rows):
        for col in range(cols):
            distance = np.sqrt((row - tx_row)**2 + (col - tx_col)**2)
            if distance > 0:
                los_path = np.linspace(tx_elevation, dem_data[row, col], num=10)
                max_obstacle_height = np.max(los_path)
                obstacle_height = max_obstacle_height - min(tx_elevation, dem_data[row, col])
                diffraction_loss_db[row, col] = knife_edge_diffraction_loss(distance, obstacle_height, wavelength)
    return diffraction_loss_db

def latlon_to_rowcol(dem_path, lat, lon):
    with rasterio.open(dem_path) as dem_dataset:
        transform = dem_dataset.transform
        crs = dem_dataset.crs

        if crs.to_string() != 'EPSG:4326':
            from pyproj import Transformer
            transformer = Transformer.from_crs('EPSG:4326', crs.to_string(), always_xy=True)
            lon, lat = transformer.transform(lon, lat)
        
        row, col = ~transform * (lon, lat)
        return int(row), int(col)

# Open the TIF file
dem_path = 'file.tif'

with rasterio.open(dem_path) as src:
    dem = src.read(1)
    row, col = latlon_to_rowcol(dem_path, 37.9899576, -78.4940777)
        
    # Transmitter parameters
    tx_location = (row, col)
    tx_power = 1  # Transmitter power in dBm
    tx_height = 1  # Transmitter height in meters
    frequency_MHz = 97.5

    coverage = simulate_rf_coverage(dem, tx_location, tx_power, frequency_MHz, tx_height)

    plt.figure(figsize=(10, 8))
    plt.imshow(coverage, cmap='viridis', extent=(0, dem.shape[1], 0, dem.shape[0]))
    plt.colorbar(label='Signal Strength (dBm)')
    plt.title('Simulated RF Coverage Map')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')

    plt.savefig('simulated_rf_coverage_map.png')
    print("The simulated RF coverage map has been saved to 'simulated_rf_coverage_map.png'.")

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
        crs=src.crs,
        transform=src.transform,
    ) as dst:
        dst.write(coverage, 1)

    print("Minimum elevation:", np.min(dem))
    print("Maximum elevation:", np.max(dem))
