import pandas as pd
import numpy as np
from scipy.integrate import trapz, simps
from hapi import *
import time

##==
##    This code calculate Plank mean absorption length [cm]
##    defined in the article below.
##    $$L_p\equiv\int_{0}^{\infty} \kappa(u)G(u) du$$
##    $$G(u)\equiv15u^3/\pi^4(e^u-1)$$
##    $$u\equiv hc\omega/kT$$
##
##    Y. Ju et al., PCI 27 (1998), 2691.
##==

##==
##....Constants in CGS unit
##==
h      = 6.62607015E-27    # Plank's constant              [erg * s]
c0     = 2.99792458E+10     # speed of light in a vacuum   [cm * s-1]
k      = 1.380649E-16      # Boltzmann's constant          [erg * K-1]
Const1 = 2.*np.pi*h*c0**2  # The first radiation constant  [j * cm2 * s-1]
Const2 = h*c0/k            # The second radiation constant [m * K]
##==
##
nu_begin = 0               # Wavenumber at fetching start  [cm-1]
nu_end   = 10000           # Wavenumber at fetching end    [cm-1]
#==
#  Here is the table for HITRAN database
#
#  species -> moleID
#    CO    ->   5
#    CO2   ->   2
#    H2O   ->   1
#    CH4   ->   6
#    H2O2  ->   25
#    H2CO  ->   20
#     SF6  ->   30
#==

def main():
    #.... Target species data (HERE ARE THE VALUES YOU SHOULD CHANGE)
    Temperature = 1300 # K
    Pressure    = 40   # atm

    species_name = "CO2"
    moleculeID = 2
    isotopologueID = 1   # Isotopologue ID of HITRAN

    #.... Data fetch
    db_begin('data')
    if not os.path.isfile(f"data/{species_name}.header"):
        fetch(species_name, moleculeID, isotopologueID, nu_begin, nu_end)

    #.... Estimate the absorption coeficcient
    #nu, absorption_coefficient = absorptionCoefficient_Voigt(((moleculeID, isotopologueID),), species_name,
    #                                                                      Environment = {'p': Pressure, 'T': Temperature}, HITRAN_units=False)
    nu, absorption_coefficient = absorptionCoefficient_Voigt(SourceTables=species_name,
                                                             Environment = {'T': Temperature, "p": Pressure}, HITRAN_units=False)
    u = calc_u(nu, Temperature)
    G = calc_G(u)
    invLP = calc_invLp(absorption_coefficient, G, u)
    df = pd.DataFrame({"wave number [cm-1]": nu,
                       "absorption coefficient [cm-1]": absorption_coefficient,
                       "u": u,
                       "G": G})
    df.to_csv(f"CSV/{species_name}_T{Temperature}_P{Pressure}.csv", index=False)
    Lp = 1/invLP
    print(f"Optical thickness of {species_name} @ T = {Temperature:.1f} [K], P = {Pressure:.1f} [atm]: {Lp:.4e} [cm]")
##==
##


##==
#.... Here is the functions
##==
def calc_u(nu, T): # [erg**2/K**2]
    return h*c0/k * nu/T

def calc_G(u):
    return 15./np.pi**4 * u**3. / ( np.exp(u)-1 )

def calc_invLp(kappa, G, u):
    return trapz( kappa*G, u )
##==
##

if __name__ == '__main__':
    main()
