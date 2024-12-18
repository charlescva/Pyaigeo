# Pyaigeo
AI, Python and Geospatial

This code simulates RF (radio frequency) coverage over a given terrain using a digital elevation model (DEM).

**Imports and Setup:**
Libraries like rasterio, numpy, matplotlib, os, numba, and pyproj are imported for handling geospatial data, numerical computations, plotting, and performance optimization.

**Knife-Edge Diffraction Loss Calculation:**
The knife_edge_diffraction_loss function calculates the diffraction loss based on the distance, height of the obstacle, and wavelength.

**Antenna Gain Pattern:**
The antenna_gain_pattern function returns the gain of different types of antennas (dipole, directional, parabolic) based on the angle.

**RF Coverage Simulation:**
The simulate_rf_coverage function calculates the RF coverage by considering path loss, atmospheric loss, diffraction loss, and antenna gain. It uses the DEM data and transmitter parameters to compute the signal strength over the terrain.

**Diffraction Loss Calculation:**
The calculate_diffraction_loss function computes the terrain-based diffraction loss using a parallelized approach with numba.

**Coordinate Transformation:**
The latlon_to_rowcol function converts latitude and longitude to row and column indices in the DEM.

**Main Function:**
The main function reads the DEM data, sets up the transmitter parameters, and calls the simulate_rf_coverage function. It then plots and saves the RF coverage map.

**Key Points:**
Performance Optimization: The use of numba for JIT compilation and parallel processing speeds up the diffraction loss calculations.
Geospatial Handling: rasterio and pyproj are used for reading DEM data and transforming coordinates.
Visualization: matplotlib is used to plot the simulated RF coverage map.
