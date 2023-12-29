# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
from main_classes.eos import SRKEOS, PREOS, RKEOS
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   T_SAT CALCULATIONS             -------------------------------------> #
acntr_factors = np.linspace(0.3, -0.3, 10)
p_rels = np.logspace(-3, 0, 4000)

if 0 not in acntr_factors:

    acntr_factors = np.append(acntr_factors, 0)
    acntr_factors.sort()

t_sat_arr = np.empty(np.array([len(acntr_factors), len(p_rels)]))
pbar = tqdm(desc="Evaluating Points", total=len(acntr_factors) * len(p_rels))
eos = PREOS

for k in range(len(acntr_factors)):

    if acntr_factors[k] == 0:

        fluid = eos(

            p_crit=1e6,
            t_crit=1e5,
            cp_ideal=850,
            m_molar=0.044,
            acntr=acntr_factors[k]

        )

    else:

        fluid = eos(

            p_crit=1e6,
            t_crit=1e5,
            cp_ideal=850,
            m_molar=0.044,
            acntr=acntr_factors[k]

        )

    for j in range(len(p_rels)):

        t_sat, zl, zv = fluid.t_sat(p_rels[j] * fluid.p_crit)
        t_sat_arr[k, j] = t_sat / fluid.t_crit
        pbar.update(1)

pbar.close()


# %%-------------------------------------   OPTIMIZE RATIO                      -------------------------------------> #
def opt_function(i_curr, minimise_max=False):

    def funct(x):

        a = x[0]
        b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i_curr] + 1)) - 1)
        t_sat_actr = np.exp(b / a * (np.power(p_rels, a) - 1))

        if minimise_max:
            return np.nanmax(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

        return np.nanmean(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

    return funct


ratios = list()
t_sat_approx = np.empty(np.array(t_sat_arr.shape))
t_sat_approx[:, :] = np.nan

pbar = tqdm(desc="Optimizing Ratios", total=len(acntr_factors))
for i in range(len(acntr_factors)):

    res = minimize(opt_function(i), np.array([0.11]))
    ratios.append(res.x[0])

    a = res.x[0]
    b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i] + 1)) - 1)
    t_sat_approx[i, :] = np.exp(b / a * (np.power(p_rels, a) - 1))

    pbar.update(1)

pbar.close()
rel_error = (t_sat_approx - t_sat_arr) / t_sat_arr


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
ax_err = ax.twinx()

for i in range(len(acntr_factors)):

    line = ax.plot(t_sat_arr[i, :], p_rels, label="$\\omega$ = {:0.2f}".format(acntr_factors[i]))[0]
    ax.plot(t_sat_approx[i, :], p_rels, "--", color=line.get_color(), alpha=0.75)
    ax_err.plot(t_sat_approx[i, :], rel_error[i, :] * 100, "-.", color=line.get_color(), alpha=0.50)

ax.set_title("Saturation Line - Approximation")
ax.set_xlabel("$T_{rel}$ [-]")
ax.set_ylabel("$p_{rel}$ [-]")
ax_err.set_ylabel("Error $T_{sat}}$ [%]")
ax.set_yscale("log")
ax.legend()
plt.show()

# %%-------------------------------------   PLOT RATIOS                        -------------------------------------> #
plt.plot(acntr_factors, ratios)
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "6 - Saturation Polinomial.png"))