"""
Function for MT 1-D
"""

import numpy as np
import math as m


def FWD_MT1D(rho, thickness, period):
    global w, Zn, dj, wj, ej, rj, re, Zj, Z, apparent_resistivity, phase
    nd = len(period)  # number of datas
    Z = np.zeros(nd, dtype=complex)
    apparent_resistivity = np.zeros(nd)
    phase = np.zeros(nd)

    for i in range(nd):
        mp = 4 * m.pi * 10 ** -7  # magnetic permeability
        w = 2 * m.pi / period[i]  # angular frequency
        nl = len(rho)  # number of layers
        impedances = np.zeros(nl, dtype=complex)

        # calculate basement impedance
        impedances[nl - 1] = np.sqrt(complex(0, w * mp * rho[nl - 1]))

        for j in range(nl - 2, -1, -1):
            R = rho[j]
            H = thickness[j]

            # iteratively calculated the apparent resistivity from bottom to top layer
            dj = np.sqrt(complex(0, w * mp / R))
            wj = dj * R
            ej = np.exp(-2 * H * dj)

            belowImpedance = impedances[j + 1]
            rj = (wj - belowImpedance) / (wj + belowImpedance)
            re = rj * ej
            Zj = wj * ((1 - re) / (1 + re))
            impedances[j] = Zj

        Z[i] = impedances[0]
        apparent_resistivity[i] = (abs(Z[i]) * abs(Z[i])) / (mp * w)
        phase[i] = np.arctan([Z[i].imag / Z[i].real]) * 180 / m.pi

    return apparent_resistivity, phase

def misfit_MT(rho_cal, rho_obs, phase_cal, phase_obs):
    global a, b, d2r, nd, e, misfit
    a = rho_cal
    b = phase_cal
    d2r = m.pi/180 #degree to radian
    nd = len(rho_obs) #length of data
    e = np.zeros(nd)
    for k in range (nd):
        e[k] =abs(m.log10(a[k]/rho_obs[k]))+abs(d2r*(b[k]-phase_obs[k]))
    return (sum(e)/nd)