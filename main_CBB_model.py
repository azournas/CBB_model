import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import constant_inds  # Assuming constant_inds.py is in the same directory
import read_excel_model
from ODES import ODEs

def main_CBB_model():
    species, S, rate_inds = read_excel_model.read_excel_model()
    k, kidcs = constant_inds.constant_inds()
    kidcs = kidcs.astype(int)
    Rknames = pd.read_excel('CBB_model.xlsx', header=None).iloc[:, 1].values
    yT = pd.read_excel('CBB_Y.xlsx')
    y0 = yT['value'].values
    tspan = [0, 180]
    tvec = np.linspace(tspan[0], tspan[1], 100)
    vcarb = []
    vcef = []
    ATP = []
    ADP = []
    NADPH = []
    NADP = []

    dydt0 = ODEs(0, y0, k[kidcs.astype(int)-1], k, rate_inds, S, Rknames, species, kidcs)

    Sol = solve_ivp(lambda t, y: ODEs(t, y, k[kidcs.astype(int)], k, rate_inds, S, Rknames, species, kidcs), 
                    tspan, y0, 
                    method='LSODA', t_eval=np.linspace(tspan[0], tspan[1], 1000),
                    rtol = 1e-12,
                    atol = 1e-12)

    Ctot = 5 * Sol.y.T[:, 0] + \
           5 * Sol.y.T[:, 2] + \
           3 * Sol.y.T[:, 4] + \
           3 * Sol.y.T[:, 5] + \
           3 * Sol.y.T[:, 7]

    plt.figure()
    plt.subplot(3, 2, 1)
    plt.plot(Sol.t, Sol.y.T[:, [0, 2, 4, 5, 7]], linewidth=2)
    plt.legend([species[i] for i in [0, 2, 4, 5, 7]], loc='best')
    plt.axvline(180)
    plt.title('CBB cycle intermediate concentrations')

    plt.subplot(3, 2, 2)
    plt.plot(Sol.t, Sol.y.T[:, [1, 3, 6, 8]], linewidth=2)
    plt.legend([species[i] for i in [1, 3, 6, 8]], loc='best')
    plt.axvline(180)
    plt.title('Energy components Concentrations')

    vcarb = Sol.y.T[:, rate_inds[1]] * k[kidcs[1]]
    plt.subplot(3, 2, 3)
    plt.plot(Sol.t, vcarb, linewidth=2)
    plt.axvline(180)
    plt.title('carboxylation rate')

    vcef = Sol.y.T[:, rate_inds[9][0]] * Sol.y.T[:, rate_inds[9][1]] * k[kidcs[9]]
    plt.subplot(3, 2, 4)
    plt.plot(Sol.t, vcef, linewidth=2)
    plt.axvline(180)
    plt.title('PSI-CEF rate')

    plt.subplot(3, 2, 5)
    plt.plot(Sol.t, Ctot)
    plt.title('total carbon atoms')

    plt.tight_layout()
    plt.show()
    foo = 1
if __name__ == "__main__":
    main_CBB_model()



