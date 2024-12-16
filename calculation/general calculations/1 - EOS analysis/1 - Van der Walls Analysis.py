# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from tqdm import tqdm
import numpy as np

# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 200
fluids = ["water", "carbonDioxide", "methane"]

# %%-------------------------------------   SHOW VAN DER WALLS                  -------------------------------------> #
v_mods = np.logspace(np.log10(0.4), np.log10(5), n_points)
t_mods = np.logspace(np.log10(0.4), np.log10(5), n_points)

xv, yv = np.meshgrid(v_mods, t_mods, indexing='ij')
zs = 8 * yv / (3 * xv - 1) - 3 / (xv ** 2)

colors = pl.cm.Oranges(np.linspace(1,0,n_points))
plt_every = 2
for i in range(int(np.floor(n_points / plt_every))):

    plt.plot(xv[:, int(i*plt_every)], zs[:, int(i*plt_every)], color=colors[int(i*plt_every)])

plt.ylim((0.1, 100))
plt.xlim((0.4, 5))
plt.xscale("log")
plt.yscale("log")
plt.xlabel("reduced volume")
plt.ylabel("reduced pressure")
plt.show()

# %%-------------------------------------   CALCULATE VAN DER WALLS             -------------------------------------> #
v_mods = np.logspace(np.log10(0.2), 3, n_points)
p_mods = np.logspace(-2, 2, n_points)
t_mods = np.empty((n_points**2, len(fluids) + 1))
t_mods[:, :] = np.nan

xv, yv = np.meshgrid(v_mods, p_mods, indexing='ij')
xv = np.array(xv).flatten()
yv = np.array(yv).flatten()
t_mods[:, 0] = (yv + 3 /(xv ** 2)) * (3 * xv - 1) / 8

# %%-------------------------------------   CALCULATE REFPROP                   -------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=n_points * n_points * len(fluids))

for i in range(len(fluids)):

    fluid = fluids[i]
    tp = ThermodynamicPoint([fluid], [1], unit_system="MOLAR BASE SI")
    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    v_crit = 1 / tp.get_variable("rho")

    for j in range(len(xv)):

        p_real = yv[j] * p_crit
        rho_real = 1 / (xv[j] * v_crit)

        try:

            tp.set_variable("P", p_real)
            tp.set_variable("rho", rho_real)
            t_real = tp.get_variable("T")

            if t_real < 0:
                raise ValueError

            t_rel = t_real / t_crit
            t_mods[j, i + 1] = t_rel

        except Exception:

            pass

        pbar.update(1)

pbar.close()

# %%-------------------------------------   EVALUATE SATURATION REFPROP         -------------------------------------> #
pbar = tqdm(desc="Evaluate Saturation Condition", total=n_points * len(fluids))
p_red_sats = np.logspace(-2, 0, n_points)

v_sat_liq_list = np.empty((n_points, len(fluids)))
v_sat_vap_list = np.empty((n_points, len(fluids)))
v_sat_liq_list[:, :] = np.nan
v_sat_vap_list[:, :] = np.nan

for i in range(len(fluids)):

    fluid = fluids[i]
    tp = ThermodynamicPoint([fluid], [1], unit_system="MOLAR BASE SI")
    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    rho_crit = tp.get_variable("rho")

    for j in range(len(p_red_sats) - 1):

        p_real = p_red_sats[j] * p_crit

        try:

            tp.set_variable("P", p_real)
            tp.set_variable("x", 0)
            rho_real = tp.get_variable("rho")

            if rho_real < 0:

                raise ValueError

            v_rel = rho_crit / rho_real
            v_sat_liq_list[j, i] = v_rel

        except:

            pass

        try:

            tp.set_variable("P", p_real)
            tp.set_variable("x", 1)
            rho_real = tp.get_variable("rho")

            if rho_real < 0:

                raise ValueError

            v_rel = rho_crit / rho_real
            v_sat_vap_list[j, i] = v_rel

        except:

            pass

        pbar.update(1)

    v_sat_liq_list[-1, i] = 1
    v_sat_vap_list[-1, i] = 1

pbar.close()

# %%-------------------------------------   EVALUATE ERRORS                     -------------------------------------> #
abs_errs = np.abs(t_mods[:, 1:] - t_mods[:, 0:1])
rel_errs = np.abs(abs_errs / t_mods[:, 1:])

# %%-------------------------------------   PLOT ERRORS                         -------------------------------------> #
ratios = [4 if i < len(fluids) else 1 for i in range(len(fluids) + 1)]
fig, axis = plt.subplots(2, len(fluids) + 1, gridspec_kw={'width_ratios':ratios})
fig.set_size_inches(15, 10)
axis = axis.T

v_rels, p_rel = np.meshgrid(v_mods, p_mods, indexing='ij')
abs_levels = np.logspace(np.log10(np.nanmin(abs_errs)), np.log10(np.nanmax(abs_errs)), 50)
rel_levels = np.logspace(np.log10(np.nanmin(rel_errs)), np.log10(np.nanmax(rel_errs)), 50)

for i in range(len(fluids)):

    abs_err = abs_errs[:, i].reshape(v_rels.shape)
    rel_err = rel_errs[:, i].reshape(v_rels.shape)

    axis[i, 0].set_title(fluids[i])
    c_bar_abs = axis[i, 0].contourf(v_rels, p_rel, abs_err, abs_levels, norm='log', cmap=pl.cm.plasma)
    c_bar_rel = axis[i, 1].contourf(v_rels, p_rel, rel_err, rel_levels, norm='log')

    for j in range(2):

        axis[i, j].plot(v_sat_liq_list[:, i], p_red_sats, color='black')
        axis[i, j].plot(v_sat_vap_list[:, i], p_red_sats, color='black')

        axis[i, j].set_xscale("log")
        axis[i, j].set_yscale("log")
        axis[i, j].set_xlabel("reduced volume")
        axis[i, j].set_ylabel("reduced pressure")

fig.colorbar(c_bar_abs, cax=axis[-1, 0])
fig.colorbar(c_bar_rel, cax=axis[-1, 1])

plt.tight_layout(pad=1)
plt.show()