#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:09:27 2024

- Define the 2D plane of magnetic moments pointing in z
- Define the observation points on a plane at height 
- Use Magpylib to calculate the magnetic field:

    in SI unit for all inputs and output

@author: zui
"""

import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt

# Define grid for magnetic moments on z=0 plane
Nx, Ny = 100, 100  # Number of points in x and y directions

x = np.linspace(-50, 50, Nx) *1e-6 # change unit to um


y = np.linspace(-50, 50, Ny) * 1e-6 # change unit to um

X, Y = np.meshgrid(x, y)

# Define magnetic moments for each dipole (e.g., all pointing in z-direction)
M0 = 9.274 * 1e-24 *10 # Magnetic moment magnitude in SI unit
deltaM = M0/100 # linear decreasing profile

sources = []


x_index=0
for x_pos, y_pos in zip(X.flatten(), Y.flatten()):
    
    M = M0 - deltaM * x_index # M linearly decrease
 
    dipole = magpy.misc.Dipole(moment=(0, 0, M))
    dipole.move((x_pos, y_pos, 0))
    sources.append(dipole)
    
    x_index = x_index+1


# Combine all dipoles into a single collection
collection = magpy.Collection(*sources)


# Define observation points on a plane at height z
z_obs = 50 * 1e-9 # Height of observation plane, change unit to nm
Nx_obs, Ny_obs = 40, 40  # Observation grid size
obs_num = Nx_obs * Ny_obs
x_obs = np.linspace(-50, 50, Nx_obs) * 1e-6 # change unit to um

y_obs = np.linspace(-50, 50, Ny_obs) * 1e-6 # change unit to um
X_obs, Y_obs = np.meshgrid(x_obs, y_obs)
observation_points = np.array([X_obs.ravel(), Y_obs.ravel(), np.full(X_obs.size, z_obs)]).T

# observers=[]
# for i in range(obs_num):
#     observers.append(observation_points[i])
    

# Calculate magnetic field at observation points
B_field = collection.getB(observation_points)
B_field = B_field.reshape(Nx_obs, Ny_obs, 3)
Bx, By, Bz = B_field.T

# Visualize the magnetic field magnitude
B_magnitude = np.linalg.norm(B_field, axis=2)


plt.figure()
plt.contourf(X_obs, Y_obs, Bz, cmap='viridis')
plt.colorbar(label='Magnetic Field Magnitude (T)')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'Magnetic Field Magnitude at z= {z_obs * 1e9: .1f} nm')
plt.show()
