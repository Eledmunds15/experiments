# Taken from...
# https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.265501

# Import Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants
T = 600  # Temperature in Kelvin
k = 8.617E-5  # Boltzmann constant in eV/K
beta = (1/(k*T))  # Inverse temperature factor

# Material Constants
Em_v = 0.61  # Vacancy migration energy in eV
Em_c = 0.35  # Core migration energy in eV
Ef_v = 0.67  # Vacancy formation energy in eV
Ef_c = 0.5  # Core formation energy in eV
D_0 = 1.51E-5  # Pre-exponential factor for diffusion (m^2/s)
vv = 9.3E13  # Frequency factor (Hz)
a = 0.285E-9  # Lattice constant in meters
r_c = 1E-9  # Core radius in meters
R_inf = 1E-3  # Assuming R_inf is equal to core radius for simplicity
omega = a**3  # Atomic volume (m^3)
sigma_a = 1E4  # Applied stress in Pascals (reduced)
c_inf = 1E6  # Vacancy concentration at infinity (vacancies/m^3)
b = 2.5E-11  # Burgers vector in meters (reduced)

# Adjust c0_v (normalized to a reasonable value)
c0_v = 1E-6  # Set c0_v to a reasonable value for vacancy concentration (vacancies/m^3)

## Functions for Climb Rate Calculation

def D_c():
    """Diffusion coefficient in the core."""
    return D_0 * np.exp(-beta * Em_c)

def D_v():
    """Diffusion coefficient in the vacancy region."""
    return D_0 * np.exp(-beta * Em_v)

def l_c(E_c_v):
    """Characteristic length scale for core diffusion."""
    term = (D_c() * r_c) * np.exp(beta * E_c_v) / (a * vv)
    return np.sqrt(term)

def l_v(E_v_c):
    """Characteristic length scale for vacancy diffusion."""
    term = (D_v() * r_c) * np.exp(beta * E_v_c) / (a * vv)
    return np.sqrt(term)

def alpha():
    """Alpha factor from Eq. (4) in the paper."""
    term = ((l_v(Em_c)**2 + (r_c**2) * np.log(R_inf/r_c))) / (2 * (l_v(Em_c)**2))
    return np.sqrt(term)

def coth(x):
    """Hyperbolic cotangent function."""
    return 1/np.tanh(x)  # coth(x) = 1/tanh(x), which avoids overflow

def climb_rate(d_j, E_c_v, E_v_c):
    """Calculate climb rate using the given equation."""
    term1 = (2 * np.pi * D_v() * c0_v) / (b * np.log(R_inf/r_c))
    term2 = (c_inf / c0_v) - np.exp(beta * sigma_a * omega)
    
    # Calculate terms involving d_j, l_c, and l_v
    alpha_value = alpha()
    argument = d_j / (2 * alpha_value * l_c(E_c_v))
    
    # Adjust scaling for term3 (modifying d_j's impact)
    term3 = (l_v(E_v_c)**2 / r_c**2) * (1 + 2 * alpha_value**2 * (argument * coth(argument) - 1))
    
    return term1 * term2 * term3

def equilibrium_climb_velocity():
    """Calculate equilibrium climb velocity (v_eq)"""
    term1 = (2 * np.pi * D_v() * c0_v) / (b * np.log(R_inf/r_c))
    term2 = (c_inf / c0_v) - np.exp(beta * sigma_a * omega)
    
    # The equilibrium climb velocity
    v_eq = term1 * term2
    
    return v_eq

# Calculate climb rates for both values of E_v_c
v_eq = equilibrium_climb_velocity()

# Simplified test for climb rate calculation
def test_climb_rate():
    E_c_v = 0.4  # Core-vacancy exchange energy barrier in eV
    E_v_c = 0.4  # Vacancy-core exchange energy barrier in eV

    # Test 2: Test equilibrium climb velocity
    v_eq = equilibrium_climb_velocity()
    print(f"Test 1 Equilibrium Climb Velocity: {v_eq} m/s")

    # Test 3: Check climb rate for extreme small and large values of d_j
    d_j_small = a*10  # Small inter-jog distance in meters (1 pm)
    d_j_large = a*10E3   # Large inter-jog distance in meters (1 µm)

    v_same = climb_rate(a, E_c_v, E_v_c)
    v_small = climb_rate(d_j_small, E_c_v, E_v_c)
    v_large = climb_rate(d_j_large, E_c_v, E_v_c)

    print(f"Test 2 Climb Rate (d_j == a): {v_same} m/s")
    print(f"Test 3 Climb Rate (small d_j > a): {v_small} m/s")
    print(f"Test 4 Climb Rate (large d_j >> a): {v_large} m/s")

test_climb_rate()
