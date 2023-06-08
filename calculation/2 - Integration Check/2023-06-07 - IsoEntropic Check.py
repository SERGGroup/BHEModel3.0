# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
from scipy.integrate import RK45
import matplotlib.pyplot as plt
import scipy.constants
import pandas as pd
import numpy as np
import os


# %%-------------------------------------   FUNCTION DEFINITION                 -------------------------------------> #
R = scipy.constants.R
g = scipy.constants.g

fluid = "Carbon Dioxide"
tp_in = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
tp_in_mol = ThermodynamicPoint([fluid], [1], unit_system="MOLAR BASE SI")

tp_in.set_variable("T", 300)
tp_in.set_variable("P", 1e5)
tp_in_mol.set_variable("T", 300)
tp_in_mol.set_variable("P", 1e5)

m_fluid = tp_in.get_variable("rho") / tp_in_mol.get_variable("rho")
tp_curr = tp_in.duplicate()

def new_dict():

    return {

        "p": list(),
        "T": list(),
        "s": list(),
        "rho": list(),

        "dT": list(),
        "drho": list()

    }


points_simple = new_dict()
points_overall = new_dict()

def rk_funct_simple(z, y):

    p_curr = z
    rho_curr = y[0]

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("rho", rho_curr)

    dsdp = tp_curr.get_derivative("s", "P", "T")
    dsdt = tp_curr.get_derivative("s", "T", "P")
    drhodp = tp_curr.get_derivative("rho", "P", "T")
    drhodt = tp_curr.get_derivative("rho", "T", "P")

    dt = - dsdp / dsdt
    drho = drhodp - dsdp / dsdt * drhodt

    points_simple["p"].append(p_curr)
    points_simple["rho"].append(rho_curr)
    points_simple["T"].append(tp_curr.get_variable("T"))
    points_simple["s"].append(tp_curr.get_variable("s"))

    points_simple["dT"].append(dt)
    points_simple["drho"].append(drho)

    return [drho]

def rk_overall_der(z, y):

    p_curr = z
    rho_curr = y[0]

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("rho", rho_curr)

    cv = tp_curr.get_variable("cv")
    t_curr = tp_curr.get_variable("T")

    dpdrho = tp_curr.get_derivative("P", "rho", "T")
    dtdrho = tp_curr.get_derivative("T", "rho", "P")
    dpdv = - rho_curr ** 2 * dpdrho
    dtdv = - rho_curr ** 2 * dtdrho

    dpdt = tp_curr.get_derivative("P", "T", "rho")
    dtdp = tp_curr.get_derivative("T", "P", "rho")

    r_dag = - t_curr * dpdt ** 2 / dpdv

    dv = 1 / (dpdv * (1 + r_dag/cv))
    dt = dtdp + dtdv * dv

    drho = -rho_curr**2 * dv

    points_overall["p"].append(p_curr)
    points_overall["rho"].append(rho_curr)
    points_overall["T"].append(tp_curr.get_variable("T"))
    points_overall["s"].append(tp_curr.get_variable("s"))

    points_overall["dT"].append(dt)
    points_overall["drho"].append(drho)

    return [drho]


# %%-------------------------------------   RK INTEGRATION                      -------------------------------------> #
p0 = 0.2 * tp_in.RPHandler.PC
tp_in.set_variable("P", p0)
tp_in.set_variable("Q", 1)
rho0 = tp_in.get_variable("rho")

p_ratio = 2
p1 = p0 * p_ratio

points_simple = new_dict()
points_overall = new_dict()

integrator_simple = RK45(rk_funct_simple, p0, [rho0], p1)
integrator_overall = RK45(rk_overall_der, p0, [rho0], p1)

while integrator_simple.status == 'running':
    integrator_simple.step()

while integrator_overall.status == 'running':
    integrator_overall.step()

output_simple = integrator_simple.y
output_overall = integrator_overall.y

err_abs_overall = (output_overall[0] - output_simple[0])
err_rel_overall = err_abs_overall / output_simple[0]
print("the overall error was {:0.2f}%".format(err_rel_overall*100))


# %%-------------------------------------   INIT PLOT                           -------------------------------------> #
h_fig = 8
n_sub_plots = 2
fig, super_axs = plt.subplots(1, n_sub_plots)
fig.set_size_inches(n_sub_plots * h_fig * 1.2, h_fig)
i = 0


# %%-------------------------------------   PLOT CURVES                         -------------------------------------> #
ax = super_axs[i]
i += 1
ax_twn = ax.twinx()
ax.set_title(fluid, fontsize=20)

df_simple = pd.DataFrame.from_dict(points_simple)
df_overall = pd.DataFrame.from_dict(points_overall)

df_simple = df_simple.sort_values("p", ascending=True)
df_overall = df_overall.sort_values("p", ascending=True)

df_simple = df_simple.drop_duplicates(["p"])
base_x = np.linspace(min(df_simple["p"]), max(df_simple["p"]), 100)

axs = [ax, ax_twn]
keys = ["T", "rho"]
labels = ["Temperature [C]", "Density [kg/$m^3$]"]
colors = ["tab:blue", "orange"]
scales = ["log", "std"]
corrs = [1, 1]
lines = list()
overall_scatters = list()

for i in range(len(axs)):

    ax = axs[i]
    key = keys[i]
    label = labels[i].split(" [")[0]

    cs = scipy.interpolate.CubicSpline(df_simple["p"], df_simple[key])

    new_line = ax.plot(

        base_x / 1e5, cs(base_x)*corrs[i],
        label=label, color=colors[i], zorder=0

    )

    new_ovr_scatter = ax.scatter(

        df_overall["p"] / 1e5, df_overall[key] * corrs[i],
        label="$R^{\dagger}$ based", color=colors[i], marker="+", zorder=200,
        s=150

    )

    ax.set_xlabel("Pressure [bar]", fontsize=17)
    ax.set_ylabel(labels[i], fontsize=17, labelpad=12)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)

    lines.append(new_line[0])
    overall_scatters.append(new_ovr_scatter)


lines.append(overall_scatters[0])
labs = [l.get_label() for l in lines]
ax.legend(lines, labs, fontsize=15)


# %%-------------------------------------   SHOW PLOT                           -------------------------------------> #
plt.tight_layout(pad=5)
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "2 - Integration Check", "output", "isoentropic check.png"))
