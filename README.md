# PMCA Calculator

A scientific computing tool for calculating Planck Mean Absorption Coefficients (PMAC) and absorption lengths for molecular species using the HITRAN spectroscopic database.

## Overview

This project provides three main calculation tools for spectroscopic analysis:

1. **Planck Mean Absorption Coefficients**: Temperature-dependent absorption coefficients for gas species
2. **Planck Mean Absorption Length**: Optical thickness calculations based on spectral integration
3. **Polynomial Fitting**: Curve fitting tools for temperature-dependent absorption data

## File Structure

- `main.py`: Calculate Planck mean absorption coefficients vs temperature
- `calc_Lp.py`: Calculate Planck mean absorption length (optical thickness)
- `fit.py`: Polynomial fitting tool for absorption coefficient data

## Prerequisites

### Required Python Packages
- `pandas`: Data manipulation and analysis
- `numpy`: Numerical computing
- `scipy`: Scientific computing (integration, curve fitting)
- `matplotlib`: Plotting and visualization
- `hapi`: HITRAN database interface

### HITRAN Database Setup
1. Create a free [HITRAN]( https://hitran.org/) account
2. Download the `hapi.py` library from HITRAN
3. Place `hapi.py` in your project directory

## Usage

### 1. Calculate Planck Mean Absorption Coefficients

```bash
python main.py
```

**Configuration**: Edit the following variables in `main.py`:
- `species_name`: Target molecular species (e.g., 'CO', 'CO2', 'H2O')
- `moleculeID`: HITRAN molecule ID
- `isotopologueID`: HITRAN isotopologue ID
- `T_arr`: Temperature range for calculations

**Output**: `PMAC_{species}.csv` with temperature-dependent absorption coefficients

### 2. Calculate Planck Mean Absorption Length

```bash
python calc_Lp.py
```

**Configuration**: Edit the following variables in `calc_Lp.py`:
- `Temperature`: Target temperature [K]
- `Pressure`: Target pressure [atm]
- `species_name`: Target molecular species
- `moleculeID`: HITRAN molecule ID

**Output**:
- Console output with calculated absorption length
- `CSV/{species}_T{temp}_P{pressure}.csv` with spectral data

### 3. Polynomial Fitting

```bash
python fit.py
```

**Prerequisites**: Run `main.py` first to generate `PMAC_{species}.csv`

**Configuration**: Edit the following variables in `fit.py`:
- `target_speces`: Target species name
- `ORDER`: Polynomial order (default: 6)
- `T_THRESH`: Temperature threshold for fitting [K]

**Output**:
- `coeffs_{species}.dat`: Polynomial coefficients
- Visualization plot

## Physical Constants

The calculations use precise physical constants:
- Planck constant: $h = 6.62607015 \times 10^{-34}$ J·s
- Speed of light: $c_0 = 2.99792458 \times 10^8$ m/s
- Boltzmann constant: $k = 1.380649 \times 10^{-23}$ J/K
- Stefan-Boltzmann constant: $\sigma = 5.670374419 \times 10^{-8}$ W·m⁻²·K⁻⁴

## Mathematical Background

The calculations are based on:
- Planck mean absorption coefficient: $\kappa_P = \frac{\int_0^{\infty} \kappa(\nu) I_b(\nu,T) d\nu}{\int_0^{\infty} I_b(\nu,T) d\nu}$
- Planck mean absorption length: $L_p = \int_0^{\infty} \kappa(u) G(u) du$ where $u = \frac{hc\omega}{kT}$
- Weighting function: $G(u) = \frac{15u^3}{\pi^4(e^u-1)}$
- Reference: Y. Ju et al., PCI 27 (1998), 2691

## Notes

- Data is cached in the `data/` directory after first fetch
- CSV outputs are saved in `CSV/` directory (for calc_Lp.py)
- Ensure sufficient disk space for HITRAN database cache
- Temperature ranges should be physically reasonable for the species

## References
You can learn how to calculate PMCA.

1. H. Zhang, M. F. Modest, Evaluation of the Planck-mean absorption coefficients from HITRAN and HITEMP databases, Journal of Quantitative Spectroscopy and Radiative Transfer, 73 (2002) 649--653.
https://doi.org/10.1016/S0022-4073(01)00178-9
