import matplotlib.pyplot as plt

D = 0.085                       #[m]
C = 0.09                        #[m]
L = 0.15                        #[m]
mpiston = 0.45                  #[kg]
mbielle = 0.65                  #[kg]
tau = 15                        #[-]
Mair_carb = 14.5                #[kg_air/kg_fuel]
import scipy.integrate as i
import numpy as np

q = 43e6                        # [J / kg_essence]
Qtot = q / Mair_carb
Vmax = (np.pi * (D**2) *C) / 4  # [V]
gamma = 1.3                     # [-]
beta = (2 * L) / C

def dQ(theta, thetaC, deltaThetaC, m_air):
    return (Qtot * m_air * 0.5 * np.pi * np.sin((np.pi * (theta - thetaC)) / deltaThetaC)) 

def F_pied(theta, p, omega):
    return (np.pi * (D/2) ** 2) * p + mpiston * (D / 2) * (omega ** 2 ) * np.cos(theta)

def F_tete(theta, p, omega):
    return - (np.pi * (D/2) ** 2) * p + (mpiston + mbielle) * (D / 2) * (omega ** 2 ) * np.cos(theta)

def F_crit(fpied,ftete):
    fpied_max = max(fpied)
    ftete_max = max(-ftete)
    fmax = max(fpied_max, ftete_max)
    return fmax

def t_crit(F_crit, alpha):
    L_flamb = L             # [m]
    A = np.pi * (D/2) ** 2  # [m^2]
    sigma_c = 450e6         # [Pa]
    E = 2e11                # [Pa]
    if alpha == (419 / 12):
        K = 1               # [-]
    elif alpha == (131 / 12):
        K = 0.5             # [-]
    sol = np.roots([-1 / F_crit, 1 / (11 * sigma_c), ((K*L_flamb)**2) / (((np.pi)**2) * E * alpha)])
    for i in range(len(sol)) :
        if (sol[i].imag != 0 ):
            sol[i] = 0
    sol_f = np.max(sol)
    sol_f = np.sqrt(sol_f)
    return abs(sol_f)


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    s *= 10e4
    theta = np.radians(theta)
    thetaC =  - 1 * np.radians(thetaC)
    deltaThetaC = np.radians(deltaThetaC)

    n = (s*Vmax)/(8.314*303.15)     # [mol]
    m_gaz = n * 0.02897             # [kg_air]
    Q_output = 0.5 * Qtot * m_gaz * (1 - np.cos(np.pi * ((theta - thetaC) / deltaThetaC)))
    Q_output[theta < thetaC] = 0
    Q_output[theta > thetaC + deltaThetaC] = 0

    V_output = Vmax * (0.5 * (1 - np.cos(theta) + beta - np.sqrt(beta ** 2 - np.sin(theta)*np.sin(theta))) + 1 / (tau - 1))

    def dp(theta, p_):
        dQ = 0
        dV = (Vmax / 2) * np.sin(theta) * (1 + np.cos(theta) / np.sqrt(beta **2 - np.sin(theta) * np.sin(theta)))
        V_output = Vmax * (0.5 * (1 - np.cos(theta) + beta - np.sqrt(beta ** 2 - np.sin(theta)*np.sin(theta))) + 1 / (tau - 1))
        if (theta > thetaC and theta < thetaC + deltaThetaC):
            dQ = ( (Qtot * np.pi) / (2 * deltaThetaC) * (np.sin(np.pi * ((theta - thetaC) / (deltaThetaC) ))) ) * m_gaz
        return - gamma * (p_ / V_output) * dV + (gamma - 1) * dQ / V_output

    p = i.solve_ivp(dp, (theta[0], theta[-1]), [s], t_eval = theta)
    p_output = p.y[0]

    omega = (2 * np.pi * rpm) / 60                  # [rad / s] 
    F_pied_output = F_pied(theta, p_output, omega)  # [N]
    F_tete_output = F_tete(theta, p_output, omega)  # [N]

    fcrit = F_crit(F_pied_output, F_tete_output)
    t = max(t_crit(fcrit, 419 / 12), t_crit(fcrit, 131 / 12))
    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)

theta = np.linspace(-180,180,1000)
V_out,Q_out,Fp,Ft,p_out,t = myfunc(2000,1,theta,24,40)

print(t)

# Créer une figure et des axes pour les sous-graphiques
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Sous-graphique 1: V_out
axs[0, 0].plot(theta, V_out, 'g-')
axs[0, 0].set_title('V_out')
axs[0, 0].grid(True)

# Sous-graphique 2: Q_out
axs[0, 1].plot(theta, Q_out)
axs[0, 1].set_title('Q_out')
axs[0, 1].grid(True)

# Sous-graphique 3: p_out
axs[1, 0].plot(theta, p_out, 'r')
axs[1, 0].set_title('p_out')
axs[1, 0].grid(True)

# Sous-graphique 4: Fp et Ft
axs[1, 1].plot(theta, Fp, 'r', label='Fp')
axs[1, 1].plot(theta, Ft, 'g', label='Ft')
axs[1, 1].set_title('Fp and Ft')
axs[1, 1].legend()
axs[1, 1].grid(True)

# Ajuster les sous-graphiques pour éviter les chevauchements
plt.tight_layout()

# Afficher les sous-graphiques
plt.show()
