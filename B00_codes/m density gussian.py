#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:41:19 2024

@author: zui
"""

import magpylib as magpy
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()   #output start time
starttime=now
print("start time =", now)

# Constants
Bohr_magneton = 9.274 * 1e-24
M0 = Bohr_magneton * 10  # Peak magnetic moment magnitude in SI units
sigma = 0.5e-6  # Standard deviation of the Gaussian profile (in meters)

# Define grid for magnetic moments on z=0 plane
Nx, Ny = 400, 200  # Number of dipoles in x and y directions

x_positions = np.linspace(-2, 2, Nx) * 1e-6  # Change unit to um
y_positions = np.linspace(-1, 1, Ny) * 1e-6  # Change unit to um

X, Y = np.meshgrid(x_positions, y_positions)

# Gaussian profile for magnetization along the x-direction
def gaussian_profile(x, M0, sigma):
    return M0 * np.exp(-x**2 / (2 * sigma**2))

# Initialize dipoles with a Gaussian diffusion profile in the x direction
dipoles = []
dipole_map = np.zeros((Nx, Ny))

for i in range(len(x_positions)):
    magnetization = gaussian_profile(x_positions[i], M0, sigma)
    
    for j in range(len(y_positions)):
        dipole = magpy.misc.Dipole(moment=(0, 0, magnetization), position=(x_positions[i], y_positions[j], 0))
        dipole_map[i, j] = magnetization
        dipoles.append(dipole)

# Step 2: Create a collection of dipoles
collection = magpy.Collection(*dipoles)

# Define observation points on a plane at height z
z_obs = 100 * 1e-9  # Height of observation plane, change unit to nm, 50nm
Nx_obs, Ny_obs = 100, 2  # Observation grid size

x_obs = np.linspace(-4, 4, Nx_obs) * 1e-6  # Change unit to um
y_obs = np.linspace(-0.1, 0.1, Ny_obs) * 1e-6  # Change unit to um
X_obs, Y_obs = np.meshgrid(x_obs, y_obs)
observation_points = np.array([X_obs.flatten(), Y_obs.flatten(), np.full(X_obs.size, z_obs)]).T

# Calculate the magnetic field at these observation points
print('Calculate B field...')
B_field = collection.getB(observation_points)

# Reshape the magnetic field components for plotting
Bx = B_field[:, 0].reshape(X_obs.shape)  # Bx[y,x]
By = B_field[:, 1].reshape(X_obs.shape)
Bz = B_field[:, 2].reshape(X_obs.shape)

# Plot dipole density along x direction
plt.figure()
plt.contourf(X_obs*1e6, Y_obs*1e6, np.linalg.norm(B_field, axis=1).reshape(X_obs.shape), cmap='viridis')
plt.colorbar(label='Magnetic Field Magnitude (T)')
plt.xlabel('X (um)')
plt.ylabel('Y (um)')
plt.title(f'Magnetic Field Magnitude at z={z_obs*1e9:.0f} nm')
plt.show()



Z=dipole_map.T/Bohr_magneton


plt.figure(figsize=(5, 4))
plt.plot(x_positions * 1e6, Z[Ny//2,:], label='Dipole density')  # Convert x_obs to micrometers for plotting
plt.xlabel('X (um)')
plt.ylabel('Dipole Magnitude (Bohr magneton)')
plt.title('Dipole Magnitude Map')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



# Plot Bz as a function of x at fixed y and z

plt.figure(figsize=(5, 4))
plt.plot(x_obs * 1e6, Bz[1,:], label='Bz')  # Convert x_obs to micrometers for plotting
plt.xlabel('X (um)')
plt.ylabel('Bz (T)')
plt.title('Bz as a Function of X')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

