import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, e, k, m_e

# Constants
T = 1e6  # Температура плазмы в Кельвинах
n_e = 1e18  # Концентрация электронов (электроны на куб. метр)
m_ch = m_e * 1000  # Масса заряженной частицы (1000 масс электрона)
v_ch = 1e5  # Скорость заряженной частицы (м/с)
lambda_D = np.sqrt(epsilon_0 * k * T / (n_e * e**2)) #Дебаевская длина в метрах

# Скорости плазменных частиц
v_th_e = np.sqrt(k * T / m_e)  # Тепловая скорость электронов
v_th_p = np.sqrt(k * T / m_e)  # Тепловая скорость позитронов

# Сетка для симуляции 2D (в метрах)
x_grid = np.linspace(-0.1, 0.1, 100)
y_grid = np.linspace(-0.1, 0.1, 100)
X, Y = np.meshgrid(x_grid, y_grid)

# Начальные позиции частиц плазмы (электронов и позитронов)
n_particles = 2000  # Число частиц плазмы
electrons_positions = np.random.rand(n_particles, 2) * 0.2 - 0.1
positrons_positions = np.random.rand(n_particles, 2) * 0.2 - 0.1

# Начальные скорости частиц плазмы
electrons_velocities = np.random.randn(n_particles, 2) * v_th_e
positrons_velocities = np.random.randn(n_particles, 2) * v_th_p

# Функция расчёта кулоновского потенциала
def coulomb_potential(x, y, x_p, y_p, charge):
    r = np.sqrt((x - x_p)**2 + (y - y_p)**2)
    return np.where(r != 0, charge / (4 * np.pi * epsilon_0 * r), 0)

# Функция расчёта дебаевского потенциала
def debye_potential(x, y, x_ch, y_ch, charge, v_ch, t):
    r = np.sqrt((x - x_ch - v_ch * t)**2 + (y - y_ch)**2)
    return np.where(r != 0, 
           charge * np.exp(-r / lambda_D) / (4 * np.pi * epsilon_0 * r), 0)

# Шаги времени для симуляции
time_steps = np.linspace(0, 1e-8, 1000)
charged_particle_position = np.array([0.0, 0.0])
charged_particle_velocity = np.array([v_ch, 0.0])

# Создание анимации взаимодействия
fig, ax = plt.subplots(figsize=(8, 8))

def update(frame, electrons_positions, positrons_positions, 
           charged_particle_position):
    ax.clear()
    ax.set_xlim(-0.1, 0.1)
    ax.set_ylim(-0.1, 0.1)
    ax.set_title(f"Plasma Interaction at t = {time_steps[frame]:.1e} s")
    
    # Обновление положения заряженной частицы с учётом временного шага
    delta_t = time_steps[1] - time_steps[0]
    charged_particle_position += charged_particle_velocity * delta_t
    
    # Граничные условия: отражение от границ
    charged_particle_position = np.where(np.abs(charged_particle_position) > 0.1, 
                                -charged_particle_position, 
                                charged_particle_position)
    
    # Обновление позиций частиц плазмы
    electrons_positions += electrons_velocities * delta_t
    positrons_positions += positrons_velocities * delta_t
    
    # Граничные условия для плазмы
    electrons_positions = np.where(np.abs(electrons_positions) > 0.1, 
                                   -electrons_positions, electrons_positions)
    positrons_positions = np.where(np.abs(positrons_positions) > 0.1, 
                                   -positrons_positions, positrons_positions)
    
    # Расчёт потенциала на сетке
    potential = np.zeros_like(X)
    for i in range(len(X)):
        for j in range(len(Y)):
            potential[i, j] += debye_potential(X[i, j], Y[i, j], 
            charged_particle_position[0], charged_particle_position[1], e, 
            v_ch, time_steps[frame])
    
    for i in range(len(electrons_positions)):
        potential += coulomb_potential(X, Y, electrons_positions[i, 0], 
                                       electrons_positions[i, 1], -e) \
                    + debye_potential(X, Y, electrons_positions[i, 0], 
                    electrons_positions[i, 1], -e, v_ch, time_steps[frame])
    for i in range(len(positrons_positions)):
        potential += coulomb_potential(X, Y, positrons_positions[i, 0], 
                                       positrons_positions[i, 1], e) \
                    + debye_potential(X, Y, positrons_positions[i, 0], 
                    positrons_positions[i, 1], e, v_ch, time_steps[frame])
    
    # Визуализация линий эквипотенциала
    ax.contourf(X, Y, potential, levels=20, cmap='viridis')
    ax.scatter(electrons_positions[:, 0], electrons_positions[:, 1], 
               color='blue', s=1, label='Electrons')
    ax.scatter(positrons_positions[:, 0], positrons_positions[:, 1], 
               color='red', s=1, label='Positrons')
    
    ax.legend(loc='upper right')

# Создание анимации
from matplotlib.animation import FuncAnimation
ani = FuncAnimation(fig, update, frames=len(time_steps), 
                    fargs=(electrons_positions, positrons_positions, 
                           charged_particle_position), interval=10)

# Сохранение анимации
animation_path = 'debye_plasma_interaction.gif'
ani.save(animation_path, writer='imagemagick', fps=10)

plt.show()
