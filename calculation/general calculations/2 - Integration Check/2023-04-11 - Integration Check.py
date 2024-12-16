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

        "z":list(),

        "p": list(),
        "T": list(),
        "h": list(),
        "s": list(),
        "rho": list(),

        "dp": list(),
        "dT": list(),
        "drho": list(),

        "dhdp": list(),
        "dhdt": list(),

        "C0": list(),

    }


points_simple = new_dict()
points_enthalpy = new_dict()
points_overall = new_dict()

def rk_funct_simple(z, y):

    p_curr = y[0]
    h_curr = tp_in.get_variable("h") - g * z

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("h", h_curr)

    t_curr = tp_curr.get_variable("T")
    rho_curr = tp_curr.get_variable("rho")

    dp = - tp_curr.get_variable("rho") * g

    dhdp = tp_curr.get_derivative("h", "P", "T")
    dhdt = tp_curr.get_derivative("h", "T", "P")

    drhodp = tp_curr.get_derivative("rho", "P", "T")
    drhodt = tp_curr.get_derivative("rho", "T", "P")

    dt = (1 / rho_curr - dhdp) / dhdt * dp
    drho = drhodt * dt + drhodp * dp

    points_simple["z"].append(z)

    points_simple["p"].append(p_curr)
    points_simple["h"].append(h_curr)
    points_simple["s"].append(tp_curr.get_variable("s"))
    points_simple["T"].append(t_curr)
    points_simple["rho"].append(rho_curr)

    points_simple["dp"].append(dp)
    points_simple["dT"].append(dt)
    points_simple["drho"].append(drho)

    points_simple["dhdp"].append(dhdp)
    points_simple["dhdt"].append(dhdt)

    points_simple["C0"].append(tp_curr.get_variable("CP0"))

    return [dp]

def rk_enthalpy_der(z, y):

    p_curr = y[0]
    rho_curr = y[1]
    v_curr = 1/rho_curr

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("rho", rho_curr)

    h_curr = tp_curr.get_variable("h")
    t_curr = tp_curr.get_variable("T")
    cp_0 = tp_curr.get_variable("CP0")

    dpdrho = tp_curr.get_derivative("P", "rho", "T")
    dtdrho = tp_curr.get_derivative("T", "rho", "P")
    dpdv = - rho_curr ** 2 * dpdrho
    dtdv = - rho_curr ** 2 * dtdrho

    dpdt = tp_curr.get_derivative("P", "T", "rho")
    dtdp = tp_curr.get_derivative("T", "P", "rho")
    dhdp = tp_curr.get_derivative("h", "P", "T")
    dhdt = tp_curr.get_derivative("h", "T", "P")

    dp = - rho_curr * g
    dv = - (dpdt * (v_curr - dhdp) / dhdt - 1) / dpdv * dp
    dt = dtdp * dp + dtdv * dv

    drho = dp/dpdrho + dt/dtdrho

    points_enthalpy["z"].append(z)

    points_enthalpy["p"].append(p_curr)
    points_enthalpy["h"].append(h_curr)
    points_enthalpy["s"].append(tp_curr.get_variable("s"))
    points_enthalpy["T"].append(t_curr)
    points_enthalpy["rho"].append(rho_curr)

    points_enthalpy["dp"].append(dp)
    points_enthalpy["dT"].append(dt)
    points_enthalpy["drho"].append(drho)

    points_enthalpy["dhdp"].append(dhdp)
    points_enthalpy["dhdt"].append(dhdt)

    points_enthalpy["C0"].append(cp_0)

    return [dp, drho]

def rk_overall_der(z, y):

    p_curr = y[0]
    rho_curr = y[1]
    v_curr = 1/rho_curr

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("rho", rho_curr)

    h_curr = tp_curr.get_variable("h")
    t_curr = tp_curr.get_variable("T")
    cp_0 = tp_curr.get_variable("CP0")
    cp = tp_curr.get_variable("CP")

    dpdrho = tp_curr.get_derivative("P", "rho", "T")
    dtdrho = tp_curr.get_derivative("T", "rho", "P")
    dpdv = - rho_curr ** 2 * dpdrho
    dtdv = - rho_curr ** 2 * dtdrho

    dpdt = tp_curr.get_derivative("P", "T", "rho")
    dtdp = tp_curr.get_derivative("T", "P", "rho")
    dhdp = tp_curr.get_derivative("h", "P", "T")
    dhdt = tp_curr.get_derivative("h", "T", "P")

    r_dag = - t_curr * dpdt ** 2 / dpdv

    dp = - rho_curr * g
    dv = (1 - r_dag / cp) / dpdv * dp
    dt = dtdp * dp + dtdv * dv

    drho = dp/dpdrho + dt/dtdrho

    points_overall["z"].append(z)

    points_overall["p"].append(p_curr)
    points_overall["h"].append(h_curr)
    points_overall["s"].append(tp_curr.get_variable("s"))
    points_overall["T"].append(t_curr)
    points_overall["rho"].append(rho_curr)

    points_overall["dp"].append(dp)
    points_overall["dT"].append(dt)
    points_overall["drho"].append(drho)

    points_overall["dhdp"].append(dhdp)
    points_overall["dhdt"].append(dhdt)

    points_overall["C0"].append(cp_0)

    return [dp, drho]


# %%-------------------------------------   RK INTEGRATION                      -------------------------------------> #
tp_in.set_variable("P", 7.4*1e6)
tp_in.set_variable("T", 300)

rho0 = tp_in.get_variable("rho")
p0 = tp_in.get_variable("P")
v0 = 1 / rho0
depth = -10000

points_simple = new_dict()
points_enthalpy = new_dict()
points_overall = new_dict()

integrator_simple = RK45(rk_funct_simple, 0, [p0], depth, rtol=1e-08, atol=1e-10)
integrator_enthalpy = RK45(rk_enthalpy_der, 0, [p0, rho0], depth, rtol=1e-5, atol=1e-6)
integrator_overall = RK45(rk_overall_der, 0, [p0, rho0], depth, rtol=1e-6, atol=1e-7)

while integrator_simple.status == 'running':
    integrator_simple.step()

while integrator_enthalpy.status == 'running':
    integrator_enthalpy.step()

while integrator_overall.status == 'running':
    integrator_overall.step()

output_simple = integrator_simple.y
output_enthalpy = integrator_enthalpy.y
output_overall = integrator_overall.y

err_abs_enthalpy = (output_enthalpy[0] - output_simple[0])
err_rel_enthalpy = err_abs_enthalpy / output_simple[0]
print("the overall error was {:0.2f}%".format(err_rel_enthalpy*100))

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
df_enthalpy = pd.DataFrame.from_dict(points_enthalpy)
df_overall = pd.DataFrame.from_dict(points_overall)

df_simple = df_simple.sort_values("z", ascending=True)
df_enthalpy = df_enthalpy.sort_values("z", ascending=True)
df_overall = df_overall.sort_values("z", ascending=True)

df_simple = df_simple.drop_duplicates(["z"])
base_x = np.linspace(min(df_simple["z"]), max(df_simple["z"]), 100)

axs = [ax, ax_twn]
keys = ["p", "rho"]
labels = ["Pressure [bar]", "Density [kg/$m^3$]"]
colors = ["tab:blue", "orange"]
scales = ["log", "std"]
corrs = [1e-5, 1]
lines = list()
enthalpy_scatters = list()
overall_scatters = list()

for i in range(len(axs)):

    ax = axs[i]
    key = keys[i]
    label = labels[i].split(" [")[0]

    cs = scipy.interpolate.CubicSpline(df_simple["z"], df_simple[key])

    new_line = ax.plot(

        np.abs(base_x), cs(base_x)*corrs[i],
        label=label, color=colors[i], zorder=0

    )

    new_enth_scatter = ax.scatter(

        np.abs(df_enthalpy["z"]), df_enthalpy[key] * corrs[i],
        label="dh based", color=colors[i], marker="x", zorder=100,
        s=100

    )

    new_ovr_scatter = ax.scatter(

        np.abs(df_overall["z"]), df_overall[key] * corrs[i],
        label="$R^{\dagger}$ based", color=colors[i], marker="+", zorder=200,
        s=150

    )

    ax.set_xlabel("depth [m]", fontsize=17)
    ax.set_ylabel(labels[i], fontsize=17, labelpad=12)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)

    lines.append(new_line[0])
    enthalpy_scatters.append(new_enth_scatter)
    overall_scatters.append(new_ovr_scatter)

    if i == 1:

        ax.ticklabel_format(axis="y", style="plain")

    else:

        ax.set_yscale("log")

lines.append(enthalpy_scatters[0])
lines.append(overall_scatters[0])
labs = [l.get_label() for l in lines]
ax.legend(lines, labs, fontsize=15)

# %%-------------------------------------   SHOW PLOT                           -------------------------------------> #
plt.tight_layout(pad=5)
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "2 - Integration Check", "output", "adiabatic check.png"))