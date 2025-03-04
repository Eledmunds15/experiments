import numpy as np
import matplotlib.pyplot as plt
import math

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
a_rad = a/(math.sqrt(2)*2) # Atomic radius of aluminium
omega = (4/3)*math.pi*a**3  # Atomic volume (m^3)
sigma_a = 1E6  # Applied stress in Pascals (reduced)
c_inf = 1E-5  # Vacancy concentration at infinity (vacancies/m^3)
b = 2.5E-10  # Burgers vector in meters (reduced)

# Adjust c0_v (normalized to a reasonable value)
c0_v = math.exp(-beta*Ef_v)  # Set c0_v to a reasonable value for vacancy concentration (vacancies/m^3)

## Functions for Climb Rate Calculation
def D_c():
    """Diffusion coefficient in the core."""
    return D_0 * math.exp(-beta * Em_c)

def D_v():
    """Diffusion coefficient in the vacancy region."""
    return D_0 * math.exp(-beta * Em_v)

def l_c(E_c_v):
    """Characteristic length scale for core diffusion."""
    term = (D_c() * r_c) * math.exp(beta * E_c_v) / (a * vv)
    return math.sqrt(term)

def l_v(E_v_c):
    """Characteristic length scale for vacancy diffusion."""
    term = (D_v() * r_c) * math.exp(beta * E_v_c) / (a * vv)
    return math.sqrt(term)

def alpha():
    """Alpha factor from Eq. (4) in the paper."""
    term = ((l_v(Em_c)**2 + (r_c**2) * math.log(R_inf/r_c))) / (2 * (l_v(Em_c)**2))
    return math.sqrt(term)

def coth(x):
    """Hyperbolic cotangent function."""
    try:
        return (math.exp(x) + math.exp(-x))/(math.exp(x) - math.exp(-x))
    except:
        return 1

def climb_rate(d_j, E_c_v, E_v_c):
    """Calculate climb rate based on the scaled inter-jog distance."""
    
    alpha_term = alpha()
    l_v_val = l_v(E_v_c)
    l_c_val = l_c(E_c_v)

    argument = d_j/(2*alpha_term*l_c_val)

    term1 = (2*math.pi*D_v()*c0_v)/(b)
    term2 = ((c_inf/c0_v) - math.exp(beta*sigma_a*omega))
    term3 = math.log(R_inf / r_c)  # Logarithmic term
    term4 = ((l_v_val**2) / (r_c**2)) * (1 + 2 * alpha_term**2 * (argument * coth(argument) - 1))

    # Calculate the climb rate using the simplified formula
    climb = (term1 * term2) / (term3 + term4)

    return climb

def climb_rate_eq():

    term1 = 2*math.pi*D_v()*c0_v
    term2 = b*math.log(R_inf / r_c)
    term3 = (c_inf/c0_v) - math.exp(beta*sigma_a*omega)

    return term1*term3/term2

def test_climb_rate():
    E_v_c = 0.6
    
    E_c_v = 0.4

    d_j_a_log = np.logspace(0, 6, 10, a)

    d_j_values_log = d_j_a_log*a

    # Calculate the climb rate for each value of d_j/a
    climb_rate_results = []
    
    for d_j_scaled in d_j_values_log:
        climb = climb_rate(d_j_scaled, E_c_v, E_v_c)
        climb_rate_results.append((d_j_scaled, climb))
    
    v_eq = climb_rate_eq()

    print(f"Equilibrium Climb Rate (m/s): {v_eq}")

    # Output the results
    print("d_j | Climb Rate (v)")
    print("----------------------------------")
    for result in climb_rate_results:
        print(f"{result[0]:.3e} | {result[1]:.3e}")

def plot_results_unscaled():

    E_v_c_1 = 0.2
    E_v_c_2 = 0.4
    E_v_c_3 = 0.6
    E_v_c_4 = 0.8
    E_v_c_5 = 1.0

    E_c_v = 0.6

    d_j_a_log = np.logspace(0, 9, 100, a)

    d_j_values_log = d_j_a_log*a

    # Calculate the climb rate for each value of d_j
    climb_rate_results_log_1 = [climb_rate(d_j, E_c_v, E_v_c_1) for d_j in d_j_values_log]
    climb_rate_results_log_2 = [climb_rate(d_j, E_c_v, E_v_c_2) for d_j in d_j_values_log]
    climb_rate_results_log_3 = [climb_rate(d_j, E_c_v, E_v_c_3) for d_j in d_j_values_log]
    climb_rate_results_log_4 = [climb_rate(d_j, E_c_v, E_v_c_4) for d_j in d_j_values_log]
    climb_rate_results_log_5 = [climb_rate(d_j, E_c_v, E_v_c_5) for d_j in d_j_values_log]

    # Plot the results
    plt.figure(figsize=(8, 6))
    plt.plot(d_j_values_log, climb_rate_results_log_1, label='E_c_v = 0.2eV', color='red', linestyle='-', marker='o', markersize=2)
    plt.plot(d_j_values_log, climb_rate_results_log_2, label='E_c_v = 0.4eV', color='blue', linestyle='-', marker='o', markersize=2)
    plt.plot(d_j_values_log, climb_rate_results_log_3, label='E_c_v = 0.6eV', color='green', linestyle='-', marker='o', markersize=2)
    plt.plot(d_j_values_log, climb_rate_results_log_4, label='E_c_v = 0.8eV', color='purple', linestyle='-', marker='o', markersize=2)
    plt.plot(d_j_values_log, climb_rate_results_log_5, label='E_c_v = 1.0eV', color='orange', linestyle='-', marker='o', markersize=2)

    # Set the x-axis to logarithmic scale
    plt.xscale('log')
    plt.xlim([min(d_j_values_log), max(d_j_values_log)])

    # Set minor grid and major grid
    plt.grid(True, which='both', axis='both', linestyle='--', linewidth=1, alpha=0.8)
    plt.minorticks_on()  # Enable minor ticks for a finer grid
    plt.grid(True, which='minor', axis='both', linestyle=':', linewidth=0.7, alpha=0.5)  # Mini grid style

    # Labels and title
    plt.xlabel('Inter-jog Distance (d_j) [m]', fontsize=12)
    plt.ylabel('Climb Rate (v) [m/s]', fontsize=12)
    plt.title('Climb Rate vs Inter-jog Distance (Unscaled)', fontsize=14)
    plt.legend()
    
    # Show the plot
    plt.show()

def plot_results_scaled():

    E_v_c_1 = 0.2
    E_v_c_2 = 0.4
    E_v_c_3 = 0.6
    E_v_c_4 = 0.8
    E_c_v = 0.9

    d_j_a_log = np.logspace(0, 9, 100, a)

    d_j_values_log = d_j_a_log*a

    v_eq = climb_rate_eq()

    # Calculate the climb rate for each value of d_j
    climb_rate_results_log_1 = [climb_rate(d_j, E_c_v, E_v_c_1)/v_eq for d_j in d_j_values_log]
    climb_rate_results_log_2 = [climb_rate(d_j, E_c_v, E_v_c_2)/v_eq for d_j in d_j_values_log]
    climb_rate_results_log_3 = [climb_rate(d_j, E_c_v, E_v_c_3)/v_eq for d_j in d_j_values_log]
    climb_rate_results_log_4 = [climb_rate(d_j, E_c_v, E_v_c_4)/v_eq for d_j in d_j_values_log]

    # Plot the results
    plt.figure(figsize=(8, 6))
    plt.plot(d_j_a_log, climb_rate_results_log_1, label='E_c_v = 0.2eV', color='red', linestyle='-', marker='o', markersize=3)
    plt.plot(d_j_a_log, climb_rate_results_log_2, label='E_c_v = 0.4eV', color='green', linestyle='-', marker='o', markersize=3)
    plt.plot(d_j_a_log, climb_rate_results_log_3, label='E_c_v = 0.6eV', color='blue', linestyle='-', marker='o', markersize=3)
    plt.plot(d_j_a_log, climb_rate_results_log_4, label='E_c_v = 0.8eV', color='purple', linestyle='-', marker='o', markersize=3)

    # Set the x-axis to logarithmic scale
    plt.xscale('log')

    # Set minor grid and major grid
    plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.7)
    plt.minorticks_on()  # Enable minor ticks for a finer grid
    plt.grid(True, which='minor', axis='both', linestyle=':', linewidth=0.5)  # Mini grid style

    # Labels and title
    plt.xlabel('d_j/a', fontsize=12)
    plt.ylabel('v/v_eq', fontsize=12)
    plt.title('Climb Rate vs Inter-jog Distance (Scaled)', fontsize=14)
    plt.legend()
    
    # Show the plot
    plt.show()

plot_results_unscaled()