import pandas as pd
import numpy as np
from scipy.integrate import trapz
from hapi import *
import time

# ==
# ....Constants
# ==
T_ref  = 296.              # Referenced Temperature        [K]
P_ref  = 1.                # Referenced Pressure           [atm.]
h      = 6.62607015E-34    # Plank constant                [j * s]
c0     = 2.99792458E8      # speed of light in a vacuum    [m * s-1]
k      = 1.380649E-23      # Boltzmann constant            [j * K-1]
Const1 = 2.*np.pi*h*c0**2  # The first radiation constant  [j * cm2 * s-1]
Const2 = h*c0/k            # The second radiation constant [m * K]
sigma  = 5.670374419E-8    # Stefan-Boltzmann constant     [W * m-2 * K-4]
# ==
#
nu_begin = 0               # Wavenumber at fetching start  [cm-1]
nu_end   = 100000          # Wavenumber at fetching end    [cm-1]
# ==
#  species -> moleID
#    CO    ->   5
#    CO2   ->   2
#    H2O   ->   1
#    CH4   ->   6
#    H2O2  ->   25
#    H2CO  ->   20
# ==

def main():
    # .... Target species data (HERE ARE THE VALUES YOU HAVE TO CHANGE)
    species_name   = 'CO'
    moleculeID     = 5   # Molecule ID of HITRAN
    isotopologueID = 1   # Isotopologue ID of HITRAN

    # .... Data fetch
    db_begin('data')
    fetch(species_name, moleculeID, isotopologueID, nu_begin, nu_end)

    # .... Estimation temperature
    Tstt = 300.   # [K]
    Tend = 3000.  # [K]
    T_arr = np.linspace(Tstt, Tend, 100)

    kp = np.zeros(len(T_arr))
    for i in range(len(T_arr)):
        nu, abs_coeff = absorptionCoefficient_SDVoigt(((moleculeID, isotopologueID),),
                                                      species_name,
                                                      Environment={
                                                          'p': P_ref,
                                                          'T': T_arr[i]
                                                      },
                                                      HITRAN_units=False)
        kp[i] = trapz(abs_coeff * Ib_nu(nu, T_arr[i]), nu) / Ib(T_arr[i])
        print(f'kp: {kp[i]:.4E} @ {T_arr[i]:.3E} [K]')
    df = pd.DataFrame({'T[K]': T_arr, 'kp[cm-1 * atm-1]': kp})
    df.to_csv(f'PMAC_{species_name}.csv', index=False)
# ==
#


# ==
# .... Functions
# ==
def Ib(T):  # Total blackbody intensity [W * m-2]
    """
    Calculate the total blackbody intensity at temperature T.
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
    Returns
    -------
    float
        Total blackbody intensity in W * m-2.
    """
    return sigma*T**4/np.pi

def Ib_nu(nu, T): #  Black body intensity at nu [W * m-2 * cm]
    """
    Calculate the blackbody intensity at wavenumber nu and temperature T.
    Parameters
    ----------
    nu : float or np.ndarray
        Wavenumber in cm-1.
    T : float
        Temperature in Kelvin.
    Returns
    -------
    float or np.ndarray
        Blackbody intensity at wavenumber nu and temperature T in W * m-2 * cm.
    """
    c1nu = Const1 * 1E8  # [W * m2 * cm-4]
    c2nu = Const2 * 1E2  # [K * cm-1]
    return c1nu * nu**3 / (np.exp(c2nu*nu/T)-1.) / np.pi
## ==
##


if __name__ == '__main__':
    main()
