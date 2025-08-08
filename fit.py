import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ORDER = 6       # orede of fittign function
T_THRESH = 800  # Threshold temperature for converting the fitting [K]
COEF = 1E2      # COEF determine the unit,
#                 if COEF == 1E2, [atm-1 * m-1], elif COEF == 1, [atm-1 * cm-1]


def main():
    # .... Set target species
    target_speces = 'H2CO'
    file_name = f'PMAC_{target_speces}.csv'
    df = pd.read_csv(file_name)
    # .... Concatinate data at threshold temperature
    dfLow = df[df['T[K]'] <= T_THRESH*1.1]
    dfHigh = df[df['T[K]'] >= T_THRESH*0.9]
    # .... Estimate polynomials
    poptLow = get_polys(dfLow)
    poptHigh = get_polys(dfHigh)

    # ... Check your polynomials by showing the figure
    T_Low = df[df['T[K]'] <= T_THRESH]['T[K]']
    T_High = df[df['T[K]'] >= T_THRESH]['T[K]']
    kp_Low = func(T_Low, poptLow)
    kp_High = func(T_High, poptHigh)
    y_sci = np.concatenate([kp_Low, kp_High])

    plot(df['T[K]'], df['kp[cm-1 * atm-1]']*COEF, y_sci)
    output('coeffs_{}.dat'.format(target_speces), poptLow, poptHigh)


def get_polys(df):
    T = df['T[K]'].values
    kp = df['kp[cm-1 * atm-1]'].values*COEF

    popt, pcov = curve_fit(func, T, kp, p0=np.ones(ORDER))

    return popt


def func(x, *params):
    y = np.zeros_like(x)
    for i, param in enumerate(params):
        y += np.array(param * x ** i)

    return y


def plot(x, y, y_sci):
    fig, ax = plt.subplots(1,1)
    ax.scatter(x, y, alpha=0.3, label='original')
    ax.plot(x, y_sci, color='red', label='scipy curve_fit')
    plt.legend(loc='best')
    plt.show()


def output(file_name, poptLow, poptHigh):
    with open(file_name, 'w') as fw:
        fw.write('<memo>\n')
        fw.write('{}th order polynomial coefficients for Planck mean absorption coefficients\n'.format(ORDER))
        fw.write('2 polynomials for each temperature range, Higher or Lower than {} [K]\n'.format(T_THRESH))
        fw.write('In <data> section...\n')
        fw.write('T_THRESH[K]\n')
        fw.write('a0Low[cm-1 * atm-1], a0High[cm-1 * atm-1]\n')
        fw.write('a1Low[cm-1 * atm-1], a1High[cm-1 * atm-1]\n')
        fw.write('.............\n')
        fw.write('a5Low[cm-1 * atm-1], a5High[cm-1 * atm-1]\n')
        fw.write('<memo>\n')
        fw.write('<data>\n')
        fw.write('{:.0e}\n'.format(T_THRESH))
        for pl, ph in zip(poptLow, poptHigh):
            fw.write('{:.6e},{:.6e}\n'.format(pl, ph))
        fw.write('<data>\n')


if __name__ == '__main__':
    main()
