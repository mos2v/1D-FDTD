import numpy as np
import matplotlib.pyplot as plt

def signal(time, pulse_width, pulse_delay, omega0, amplitude = 1.0):
    return amplitude * np.exp(- ((time - pulse_delay) / pulse_width) ** 2) * np.sin(omega0 * time)

eps0 = 8.8541878128e-12
mu0 = 1.256637062e-6
c0 = 2.99792458e8
imp0 = np.sqrt(mu0 / eps0)

simulation_size = 10e-6
step_size = 5e-9
N_space_cells = int(simulation_size/ step_size)
print(f"there are {N_space_cells} FDTD cells")


dt = step_size / c0
simulation_time = 10e-12
N_time_steps = int(simulation_time / dt)
print(f"there are {N_time_steps} FDTD time steps")


Ex = np.zeros(N_space_cells)
Hz = np.zeros(N_space_cells)
eps = np.ones(N_space_cells)
refractive_index = np.sqrt(eps)

h_coeff = dt / (mu0 + step_size) 
e_coeff = dt / (eps0 * eps * step_size)

c = c0 / refractive_index[0]
c_ = c0 / refractive_index[-1]
a = (c * dt - step_size) / (c * dt + step_size)
a_ = (c_ * dt - step_size) / (c_ * dt + step_size)


center_wavelength = 550e-9
omega0 = 2 * np.pi * c0 / center_wavelength
pulse_width = 20e-15
pulse_delay = 4 * pulse_width

time = np.linspace(0, simulation_time)


for n in range(N_time_steps):
    Hz_prev = Hz.copy()
    Ex_prev = Ex.copy()
    

    # for j in range(0, N_space_cells - 1):
    #     Hz[j] = Hz_prev[j] + h_coeff * (Ex[j+1] - Ex[j])

    Hz[:N_space_cells-1] = Hz_prev[:N_space_cells-1] + h_coeff * (Ex[1:] - Ex[:N_space_cells-1])
    #Hz = [Hz_prev[j] + dt / (mu0 + step_size) * (Ex[j+1] - Ex[j]) for j in range(0, N_space_cells - 1)]    

    # for j in range(1, N_space_cells - 1):
    #     Ex[j] = Ex_prev[j] +  e_coeff[j] * (Hz[j] - Hz[j - 1])

    Ex[1:N_space_cells-1] = Ex_prev[1:N_space_cells-1] +  e_coeff[1:N_space_cells-1] * (Hz[1:N_space_cells-1] - Hz[:N_space_cells - 2])

    #Ex = [Ex_prev[j] + dt / (eps0 * eps[j] * step_size) * (Hz[j] - Hz[j - 1]) for j in range(1, N_space_cells - 1)]  

    Ex[0] = Ex_prev[1] + a * (Ex[1] - Ex_prev[0])
    Ex[-1] = Ex_prev[-2] + a_ * (Ex[-2] - Ex_prev[-1])


    if n % 100 == 0:
        print(n)
        print(np.min(Ex), np.max(Ex))    
