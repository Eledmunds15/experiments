import numpy as np
import matplotlib.pyplot as plt
import math

# Constants (these will now be passed as parameters)
def get_constants():
    return {
        'k': 8.617E-5,  # Boltzmann constant in eV/K
        'D_0': 1.2E-5,  # Pre-exponential factor for diffusion (m^2/s) for Fe
        'vv': 1.0E13,  # Frequency factor (Hz) for Fe
        'a': 2.866E-10,  # Lattice constant in meters (for Fe)
        'Em_v': 0.74,  # Vacancy migration energy in eV (for Fe)
        'Em_c': 0.55,  # Core migration energy in eV (average value for Fe)
        'Ef_v': 1.05,  # Vacancy formation energy in eV (for Fe)
        'Ef_c': 0.65,  # Core formation energy in eV (for Fe)
        'sigma_a': 1E6,  # Applied stress in Pascals (reduced)
        'r_c': 1E-9,  # Core radius in meters
        'R_inf': 1E-3,  # Assuming R_inf is equal to core radius for simplicity
        'c_inf': 1E-5,  # Vacancy concentration at infinity (vacancies/m^3)
        'b': 2.5E-10,  # Burgers vector in meters (reduced)
    }

def climb_rate(d_j, E_c_v, E_v_c, T, constants):
    """Calculate climb rate based on the scaled inter-jog distance."""
    
    k = constants['k']
    D_0 = constants['D_0']
    vv = constants['vv']
    a = constants['a']
    sigma_a = constants['sigma_a']
    c_inf = constants['c_inf']
    r_c = constants['r_c']
    R_inf = constants['R_inf']
    b = constants['b']
    
    beta = 1 / (k * T)  # Inverse temperature factor
    c0_v = math.exp(-beta * constants['Ef_v'])  # Vacancy concentration at zero stress
    
    def D_c():
        """Diffusion coefficient in the core."""
        return D_0 * math.exp(-beta * constants['Em_c'])

    def D_v():
        """Diffusion coefficient in the vacancy region."""
        return D_0 * math.exp(-beta * constants['Em_v'])

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
        term = ((l_v(constants['Em_c'])**2 + (r_c**2) * math.log(R_inf / r_c))) / (2 * (l_v(constants['Em_c'])**2))
        return math.sqrt(term)

    def coth(x):
        """Hyperbolic cotangent function."""
        try:
            return (math.exp(x) + math.exp(-x)) / (math.exp(x) - math.exp(-x))
        except:
            return 1

    # Now calculate the climb rate
    alpha_term = alpha()
    l_v_val = l_v(E_v_c)
    l_c_val = l_c(E_c_v)
    argument = d_j / (2 * alpha_term * l_c_val)

    term1 = (2 * math.pi * D_v() * c0_v) / b
    term2 = ((c_inf / c0_v) - math.exp(beta * sigma_a * constants['omega']))
    term3 = math.log(R_inf / r_c)  # Logarithmic term
    term4 = ((l_v_val**2) / (r_c**2)) * (1 + 2 * alpha_term**2 * (argument * coth(argument) - 1))

    climb = (term1 * term2) / (term3 + term4)

    return climb

def plot_results_temperature():
    constants = get_constants()  # Retrieve all constants from a single function

    # Parametric changes
    temperature_values = [300, 400, 500, 600, 700]  # Different temperature values to explore (in Kelvin)

    E_v_c_values = [0.2, 0.4, 0.6, 0.8, 1.0]  # Different E_v_c values for comparison

    d_j_a_log = np.logspace(0, 9, 100, constants['a'])
    d_j_values_log = d_j_a_log * constants['a']

    for T in temperature_values:
        for E_v_c in E_v_c_values:
            # Calculate the climb rate for each combination of parameters
            climb_rate_results = [climb_rate(d_j, E_c_v=0.6, E_v_c=E_v_c, T=T, constants=constants) for d_j in d_j_values_log]

            # Plot results for each temperature
            plt.figure(figsize=(8, 6))
            plt.plot(d_j_values_log, climb_rate_results, label=f'E_v_c = {E_v_c}eV, T = {T}K', marker='o', markersize=2)

            # Set the x-axis to logarithmic scale
            plt.xscale('log')
            plt.xlim([min(d_j_values_log), max(d_j_values_log)])

            # Add a horizontal line for equilibrium climb rate (optional)
            plt.axhline(y=climb_rate_eq(), color='black', linestyle='-', label='Equilibrium Climb Velocity')

            # Set minor grid and major grid
            plt.grid(True, which='both', axis='both', linestyle='--', linewidth=1, alpha=0.8)
            plt.minorticks_on()
            plt.grid(True, which='minor', axis='both', linestyle=':', linewidth=0.7, alpha=0.5)

            # Labels and title
            plt.xlabel('Inter-jog Distance (d_j) [m]', fontsize=12)
            plt.ylabel('Climb Rate (v) [m/s]', fontsize=12)
            plt.title(f'Climb Rate vs Inter-jog Distance (T = {T}K)', fontsize=14)
            plt.legend()

            # Show the plot
            plt.show()

# Call the function to explore the effect of temperature
plot_results_temperature()
