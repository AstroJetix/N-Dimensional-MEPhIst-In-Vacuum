import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

q = 1.6e-19
epsilon_0 = 8.854e-12
lambda_D = 0.01
v_z = 1e4

rho = np.linspace(-0.05, 0.05, 100)
z = np.linspace(-0.1, 0.1, 100)
rho_grid, z_grid = np.meshgrid(rho, z)

def debye_potential(rho_grid, z_grid, t):
    r_prime = np.sqrt(rho_grid**2 + (z_grid - v_z * t)**2)
    phi = (q / (4 * np.pi * epsilon_0)) * np.exp(-r_prime / lambda_D) / r_prime
    return phi

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_title("Dynamic Equipotential Lines of Debye Potential Over Time")
ax.set_xlabel("ρ (meters)")
ax.set_ylabel("z (meters)")

time_steps = np.linspace(0, 5e-6, 50)

def update(t):
    ax.clear()
    ax.set_title("Dynamic Equipotential Lines of Debye Potential Over Time")
    ax.set_xlabel("ρ (meters)")
    ax.set_ylabel("z (meters)")
    phi = debye_potential(rho_grid, z_grid, t)
    contour = ax.contourf(rho_grid, z_grid, phi, levels=20, cmap="viridis")

ani = FuncAnimation(fig, update, frames=time_steps, interval=100)

animation_path = "debye_potential_animation.gif"
ani.save(animation_path, writer="imagemagick", fps=10)

animation_path
