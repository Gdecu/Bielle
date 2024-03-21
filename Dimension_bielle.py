import scipy as sc
from scipy.integrate import odeint
import numpy as np

D = 0.085                       #[m]
C = 0.09                        #[m]
L = 0.15                        #[m]
mpiston = 0.45                  #[kg]
mbielle = 0.65                  #[kg]
tau = 15                        #[-]

Mair_carb = 14.5                # [kg_air/kg_fuel]
q = 43e6                        # [J / kg_air]

Vmax = (np.pi * (D**2) *C) / 4  # [V]
gamma = 1.3                     # [-]

beta = (2 * L) / C # test


def Qtot(s):
    n = (s*Vmax)/(8.314*303) # PK T = 303 K ??!
    m_gaz = n * 0.029
    return m_gaz * q

def Volume (theta):
    return Vmax * (0.5 * (1 - np.cos(theta) + beta - np.sqrt(beta ** 2 - np.sin(theta)*np.sin(theta))) + 1 / (tau - 1))

def dV(theta):
    return Vmax * (np.sin(theta) + np.sin(theta) * np.cos(theta) * 1 / ( np.sqrt(beta**2 - (np.sin(theta))**2)))
# Faut peut-être diviser par 2

def dQ(theta, thetaC, deltaThetaC, s):
    return (Qtot(s) * np.pi * np.sin((np.pi * (theta - thetaC)) / deltaThetaC)) 

def dp(p, theta, thetaC, deltaThetaC, s):
    dQ_ = dQ (theta, thetaC, deltaThetaC, s)
    #dp = - gamma * (p / Vmax) * dV(theta) + (gamma - 1) * (dQ_ / Vmax)
    # A verif avec d'autres pour voir comment ils ont fait
    if (theta < -2*np.pi  and theta > 2*np.pi):
        dp = 0
    elif (theta < thetaC):
        dp = - gamma * (p / Vmax) * dV(theta)
    elif (theta < thetaC + deltaThetaC):
        dp = - gamma * (p / Vmax) * dV(theta) + (gamma - 1) * (dQ_ / Vmax)
    elif (theta < np.pi):                        
        dp = - gamma * (p / Vmax) * dV(theta)
    elif (theta < 2*np.pi):                      
        dp = s-p
    return dp



def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    p0 = s * 10e5
    p = sc.integrate.solve_ivp(dp, (theta[0], theta[-1]), [p0], t_eval =theta)
    p_output = p.y[0] # Pas trop capté mais vsy pas lchoix



    (V_output, Q_output, F_pied_output, F_tete_output, p_output, t) = (0,0,0,0,0,0)
    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)