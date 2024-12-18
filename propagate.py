import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os
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
    print("Simulating the RF Propogation...")
    c = 3e8
    frequency_hz = frequency_MHz * 1e6
    wavelength = c / frequency_hz

    tx_row, tx_col = tx_location
    tx_elevation = dem_data[tx_row, tx_col] + tx_height

    rows, cols = np.indices(dem_data.shape)
    distances = np.sqrt((rows - tx_row)**2 + (cols - tx_col)**2)
    altitudes = dem_data + tx_height
    
    path_loss_db = 20 * np.log10(distances + 1e-6) + 20 * np.log10(frequency_hz) - 147.55
    atmospheric_loss_db = 0.02 * altitudes

    diffraction_loss_db = calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength)
    
    elevation_angles = np.degrees(np.arctan2(rows - tx_row, distances))
    azimuth_angles = (np.degrees(np.arctan2(cols - tx_col, rows - tx_row)) - direction) % 360
    print("Applying Gain Pattern for " + antenna_type + " antenna...")
    antenna_gain_db = np.vectorize(antenna_gain_pattern)(antenna_type, elevation_angles - angle_of_attack, azimuth_angles, hpbw, fnbw)

    total_loss_db = path_loss_db + atmospheric_loss_db + diffraction_loss_db - antenna_gain_db
    rf_coverage = tx_power_dbm - total_loss_db

    return rf_coverage

@njit(parallel=True)
def calculate_diffraction_loss(dem_data, tx_elevation, tx_row, tx_col, wavelength):
    print("Calculating Diffraction Terrain based Loss...")
    rows, cols = dem_data.shape
    diffraction_loss_db = np.zeros((rows, cols), dtype=np.float32)

    for row in prange(rows):
        for col in prange(cols):
            distance = np.sqrt((row - tx_row)**2 + (col - tx_col)**2)
            if distance > 0:
                num_samples = int(distance) * 10
                path_heights = np.linspace(tx_elevation, dem_data[row, col], num=num_samples)
                
                row_indices = np.linspace(tx_row, row, num=num_samples).astype(np.int32)
                col_indices = np.linspace(tx_col, col, num=num_samples).astype(np.int32)
                
                max_obstacle_height = -np.inf
                for i in range(num_samples):
                    r_idx = row_indices[i]
                    c_idx = col_indices[i]
                    if dem_data[r_idx, c_idx] > max_obstacle_height:
                        max_obstacle_height = dem_data[r_idx, c_idx]
                
                obstacle_height = max_obstacle_height - np.min(path_heights)
                diffraction_loss_db[row, col] = knife_edge_diffraction_loss(distance, obstacle_height, wavelength)    
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

def main():
    dem_path = 'file2x.tif'
    
    try:
        with rasterio.open(dem_path) as src:
            meta = src.meta
            dem = src.read(1)
            dem_transform = src.transform
            dem_crs = src.crs
            print('Coord Ref System:', dem_crs.to_string())
            
            latitude = 37.985041083363242                    
            longitude = -78.480348470078056

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
                crs=dem_crs,
                transform=dem_transform,
            ) as dst:
                dst.write(coverage, 1)

            print("Done")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
