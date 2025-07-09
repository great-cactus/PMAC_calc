# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a scientific computing project for calculating Planck Mean Absorption Coefficients (PMAC) for molecular species using the HITRAN database. The codebase performs spectroscopic calculations related to radiation absorption by gases at various temperatures and pressures.

## Core Architecture

The project consists of three main Python modules:

1. **main.py**: Calculates Planck mean absorption coefficients as a function of temperature
   - Fetches spectroscopic data from HITRAN database using `hapi` library
   - Computes absorption coefficients using Voigt profiles
   - Integrates over wavelength spectrum to get Planck mean values
   - Outputs results to CSV files

2. **calc_Lp.py**: Computes Planck mean absorption length (optical thickness)
   - Implements the mathematical formulation from Y. Ju et al., PCI 27 (1998), 2691
   - Calculates dimensionless parameter `u = hcω/kT` and weighting function `G(u)`
   - Uses CGS units for physical constants
   - Outputs spectral data and final absorption length

3. **fit.py**: Polynomial fitting tool for absorption coefficient data
   - Fits 6th-order polynomials to temperature-dependent absorption data
   - Splits data at threshold temperature (800K) for better fitting
   - Generates visualization plots and coefficient files

## Dependencies

The project requires:
- `pandas` for data handling
- `numpy` for numerical computations
- `scipy` for integration and curve fitting
- `matplotlib` for visualization
- `hapi` library for HITRAN database access

## Key Configuration

- **Species selection**: Modify `moleculeID` and `isotopologueID` in main scripts
- **Temperature ranges**: Adjust `T_arr` in main.py or `Temperature` in calc_Lp.py
- **Pressure conditions**: Set `P_ref` (main.py) or `Pressure` (calc_Lp.py)
- **Wavenumber range**: Configure `nu_begin` and `nu_end` for spectral coverage

## Data Flow

1. HITRAN database access requires account setup and `hapi.py` download
2. Spectroscopic data is fetched and cached in `data/` directory
3. Absorption coefficients computed using Voigt line profiles
4. Temperature-dependent calculations produce CSV output files
5. Polynomial fitting generates coefficient files for further use

## Physical Constants

The code uses precise physical constants in both SI and CGS units:
- Planck constant, speed of light, Boltzmann constant
- Radiation constants for blackbody calculations
- Stefan-Boltzmann constant for intensity calculations

## Output Files

- `PMAC_{species}.csv`: Temperature-dependent absorption coefficients
- `{species}_T{temp}_P{pressure}.csv`: Spectral absorption data
- `coeffs_{species}.dat`: Polynomial fitting coefficients