""" Inversion process using SOS Algorithm """


import numpy as np
import random as r
from MT_1D import FWD_MT1D, misfit_MT
import matplotlib.pyplot as plt

# Create synthetic data

r_layer = [100, 10, 1000]  # resistivity of each layer
t_layer = [600, 1200]  # thickness of each layer
period = np.logspace(-3, 3, 55)  # period from 10^-3 to 10^3 divide into 50 points
ndata = len(period)
rho_data, phase_data = FWD_MT1D(r_layer, t_layer, period)

# Estimated number of layers (predicted)
lyr = 3
npar = 2 * lyr - 1  # number of parameter model

# Model space
npop = 50  # number of populations
ngen = 50  # number of iteration/generation
r_min = 1  # minimum value of resistivity
r_max = 2000  # maximum value of resistivity
h_min = 200  # minimum value of thickness
h_max = 2000  # maximum value of thickness

# Create initial model
res = np.zeros((npop, lyr))
thick = np.zeros((npop, lyr - 1))
model = []
for ip in range(npop):
    for im in range(lyr):
        res[ip][im] = r_min + r.random() * (r_max - r_min)
    for im in range(lyr - 1):
        thick[ip][im] = h_min + r.random() * (h_max - h_min)
    model.append(np.hstack(np.hstack([res[ip], thick[ip]])))

# calculated the apparent resistivities, phases, and misfit of each initial model
rho_mod = np.zeros((npop, ndata))
phase_mod = np.zeros((npop, ndata))
misfit = np.zeros(npop)
for ip in range(npop):
    rho_mod[ip], phase_mod[ip] = FWD_MT1D(res[ip], thick[ip], period)
    misfit[ip] = misfit_MT(rho_mod[ip], rho_data, phase_mod[ip], phase_data)

# inversion process
misfit_inv = np.zeros(ngen)
for it in range(ngen):
    for ip in range(npop):
        idx = np.argmin(misfit)  # index that shows model with minimum misfit
        res_best = res[idx]
        thick_best = thick[idx]
        model_best = model[idx]

        # mutualism phase
        j = r.randint(0, npop - 1)
        k = r.randint(0, npop - 1)
        while j == k | ip == j | ip == k:
            j = r.randint(0, npop - 1)
            k = r.randint(0, npop - 1)

        # model tested and targeted model
        mod_tes = model[ip], model[j]
        misfit_tes = misfit[ip], misfit[j]
        mod_tar = model[k]

        mod_mut = np.zeros((2, npar))
        rho_mut = np.zeros((2, ndata))
        phase_mut = np.zeros((2, ndata))
        misfit_mut = np.zeros(2)
        for i in range(2):
            for im in range(npar):
                mod_mut[i][im] = mod_tes[i][im] + r.random() * (mod_tar[im] - (mod_tes[0][im] + mod_tes[1][im]) / 2)
                if im < lyr:
                    if mod_mut[i][im] < r_min:
                        mod_mut[i][im] = r_min
                    if mod_mut[i][im] > r_max:
                        mod_mut[i][im] = r_max
                else:
                    if mod_mut[i][im] < h_min:
                        mod_mut[i][im] = h_min
                    if mod_mut[i][im] > h_max:
                        mod_mut[i][im] = h_max

            rho_mut[i], phase_mut[i] = FWD_MT1D(mod_mut[i][0:lyr], mod_mut[i][lyr:npar], period)
            misfit_mut[i] = misfit_MT(rho_mut[i], rho_data, phase_mut[i], phase_data)

        if misfit_mut[0] < misfit[ip]:
            model[ip] = mod_mut[0]
            misfit[ip] = misfit_mut[0]
            res[ip] = mod_mut[0][0:lyr]
            thick[ip] = mod_mut[0][lyr:npar]
            rho_mod[ip] = rho_mut[0]
            phase_mod[ip] = phase_mut[0]
        if misfit_mut[1] < misfit[j]:
            model[j] = mod_mut[1]
            misfit[j] = misfit_mut[1]
            res[j] = mod_mut[1][0:lyr]
            thick[j] = mod_mut[1][lyr:npar]
            rho_mod[j] = rho_mut[1]
            phase_mod[j] = phase_mut[1]

        # commensalism phase
        j = r.randint(0, npop - 1)
        while j == ip:
            j = r.randint(0, npop - 1)

        mod_com = np.zeros(npar)
        for im in range(npar):
            mod_com[im] = model[ip][im] + (-1 + 2 * r.random()) * (model_best[im] - model[j][im])
            if im < lyr:
                if mod_com[im] < r_min:
                    mod_com[im] = r_min
                if mod_com[im] > r_max:
                    mod_com[im] = r_max
            else:
                if mod_com[im] < h_min:
                    mod_com[im] = h_min
                if mod_com[im] > h_max:
                    mod_com[im] = h_max
        rho_com, phase_com = FWD_MT1D(mod_com[0:lyr], mod_com[lyr:npar], period)
        misfit_com = misfit_MT(rho_com, rho_data, phase_com, phase_data)

        if misfit_com < misfit[ip]:
            model[ip] = mod_com
            misfit[ip] = misfit_com
            res[ip] = mod_com[0:lyr]
            thick[ip] = mod_com[lyr:npar]
            rho_mod[ip] = rho_com
            phase_mod[ip] = phase_com

        # parasitism phase
        j = r.randint(0, npop - 1)
        while j == ip:
            j = r.randint(0, npop - 1)

        mod_par = model[ip]
        p = r.randint(0, npar - 1)
        if p < lyr:
            mod_par[p] = r_min + r.random() * (r_max - r_min)
        else:
            mod_par[p] = h_min + r.random() * (h_max - h_min)
        rho_par, phase_par = FWD_MT1D(mod_par[0:lyr], mod_par[lyr:npar], period)
        misfit_par = misfit_MT(rho_par, rho_data, phase_par, phase_data)

        if misfit_par < misfit[j]:
            model[j] = mod_par
            misfit[j] = misfit_par
            res[j] = mod_par[0:lyr]
            thick[j] = mod_par[lyr:npar]
            rho_mod[j] = rho_par
            phase_mod[j] = phase_par

    # best model
    inv = np.argmin(misfit)
    model_inv = model[inv]
    res_inv = res[inv]
    thick_inv = thick[inv]
    rho_inv = rho_mod[inv]
    phase_inv = phase_mod[inv]
    misfit_inv[it] = misfit[inv]

# Plotting
plot_res_obs = np.zeros(2 * lyr)
plot_res_inv = np.zeros(2 * lyr)
plot_thi_obs = np.zeros(2 * lyr)
plot_thi_inv = np.zeros(2 * lyr)
for i in range(lyr):
    plot_res_obs[2 * i] = r_layer[i]
    plot_res_obs[2 * i + 1] = r_layer[i]
    plot_res_inv[2 * i] = res_inv[i]
    plot_res_inv[2 * i + 1] = res_inv[i]
plot_thi_obs[0] = 0
plot_thi_obs[2 * lyr - 1] = 10000
plot_thi_inv[0] = 0
plot_thi_inv[2 * lyr - 1] = 10000
for i in range(lyr - 1):
    plot_thi_obs[2 * i + 1] = plot_thi_obs[2 * i] + t_layer[i]
    plot_thi_obs[2 * i + 2] = plot_thi_obs[2 * i] + t_layer[i]
    plot_thi_inv[2 * i + 1] = plot_thi_inv[2 * i] + thick_inv[i]
    plot_thi_inv[2 * i + 2] = plot_thi_inv[2 * i] + thick_inv[i]

# figure 1
# App. resistivity vs Period
plt.figure(figsize=(15, 10))
plt.axes([0.1, 0.45, 0.475, 0.45])
plt.loglog(period, rho_inv, color='red', linewidth=4, label='calculated')
plt.loglog(period, rho_data, '.', color='blue', markersize=14, label='observed')
plt.xlim([0.001, 1000])
plt.ylim([1, 10000])
plt.xlabel('Periode (sec)', fontsize=15)
plt.ylabel('App. Resistivity (Ohm.m)', fontsize=15)
plt.legend(prop={'size': 14})

# Phase vs Period
plt.axes([0.1, 0.1, 0.475, 0.25])
plt.semilogx(period, phase_inv, color='red', linewidth=4, label='calculated')
plt.semilogx(period, phase_data, '.', color='blue', markersize=14, label='observed')
plt.xlim([0.001, 1000])
plt.ylim([0, 90])
plt.xlabel('Periode (sec)', fontsize=15)
plt.ylabel('Phase (deg)', fontsize=15)
plt.legend(prop={'size': 14})

# resist vs depth
plt.axes([0.65, 0.1, 0.275, 0.8])
plt.step(plot_res_obs, plot_thi_obs, color='blue', linestyle='--', linewidth=4, label='Observe')
plt.step(plot_res_inv, plot_thi_inv, color='red', linestyle='-', linewidth=4, label='Inverse')
plt.xscale('log')
plt.xlim([1, 10 ** 4])
plt.ylim([0, 5000])
plt.gca().invert_yaxis()
plt.xlabel('Resistivity (Ohm.m)', fontsize=15)
plt.ylabel('Depth (m)', fontsize=15)
plt.legend(prop={'size': 14})
plt.show()

# figure 2
x = list(range(ngen))
plt.figure(figsize=(8, 6))
plt.plot(x, misfit_inv)
plt.xlim([0, ngen])
plt.title('Misfit Curve', fontsize=25)
plt.xlabel('Iteration', fontsize=15)
plt.ylabel('Misfit', fontsize=15)
