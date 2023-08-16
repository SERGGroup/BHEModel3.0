# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.support.simple_integrator import SimpleIntegrator
from main_classes.eos import SRKEOS, PREOS, convert_cp_coeff
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
from scipy.integrate import RK45
import matplotlib.pyplot as plt
import scipy.constants
import pandas as pd
import numpy as np
import os

def new_dict():

    return {

        "p": list(),
        "T": list(),
        "s": list(),
        "rho": list(),
        "h": list(),

        "cv": list(),
        "R": list(),

    }

fluid = "Carbon Dioxide"

# %%-------------------------------------   RP FUNCTION DEFINITION              -------------------------------------> #
R = scipy.constants.R
g = scipy.constants.g

tp_in = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
t_ref = tp_in.RPHandler.TC*2
p_ref = tp_in.RPHandler.PC*0.01
tp_ref = tp_in.duplicate()
tp_ref.set_variable("P", t_ref)
tp_ref.set_variable("T", p_ref)

tp_curr = tp_in.duplicate()
points_overall = new_dict()

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
    drho = -rho_curr**2 * dv

    return [drho]

def rk_overall_der_T(z, y):

    p_curr = z
    t_curr = y[0]

    tp_curr.set_variable("p", p_curr)
    tp_curr.set_variable("T", t_curr)

    cv = tp_curr.get_variable("cv")
    rho_curr = tp_curr.get_variable("rho")

    dpdrho = tp_curr.get_derivative("P", "rho", "T")
    dtdrho = tp_curr.get_derivative("T", "rho", "P")
    dpdv = - rho_curr ** 2 * dpdrho
    dtdv = - rho_curr ** 2 * dtdrho

    dpdt = tp_curr.get_derivative("P", "T", "rho")
    dtdp = tp_curr.get_derivative("T", "P", "rho")

    r_dag = - t_curr * dpdt ** 2 / dpdv

    dv = 1 / (dpdv * (1 + r_dag/cv))
    dt = dtdp + dtdv * dv

    return [dt]


# %%-------------------------------------   SRK FUNCTION DEFINITION             -------------------------------------> #
tp_in.set_variable("T", tp_in.RPHandler.TC)
tp_in.set_variable("P", tp_in.RPHandler.PC)
cp0_crit = tp_in.get_variable("CP0")

# base_coeff = [0.79655508, 0.21270104]
# coeff = convert_cp_coeff(base_coeff, cp0_crit, tp_in.RPHandler.TC)
coeff = [4.49897751e+02, 1.66780277e+00, -1.27243808e-03, 3.90820268e-07]
fluid_eos = PREOS(

    t_crit=tp_in.RPHandler.TC, p_crit=tp_in.RPHandler.PC,
    cp_ideal=coeff, m_molar=0.04401, acntr=0.239

)

ref_state = fluid_eos.get_state(t=t_ref, p=p_ref)
curr_state = fluid_eos.get_state(t=300, p=1e5)

points_srk = new_dict()

def rk_srk_der(z, y):

    p_curr = z
    v_curr = y[0]

    curr_state.update_state(p=p_curr, v=v_curr)
    dv = 1 / (curr_state.dpdv * (1 + curr_state.r / curr_state.cv))

    return [dv]

def rk_srk_der_T(z, y):

    p_curr = z
    t_curr = y[0]

    curr_state.update_state(p=p_curr, t=t_curr)
    dt = curr_state.r / curr_state.cp / curr_state.dpdt

    return [dt]


# %%-------------------------------------   RK INTEGRATION                      -------------------------------------> #
tp_in.set_variable("T", t_ref)
tp_in.set_variable("P", p_ref)
integrate_temperature = True

p_ratio = 2
p1 = p_ref * p_ratio

points_srk = new_dict()
points_overall = new_dict()

if integrate_temperature:

    integrator_srk = SimpleIntegrator(rk_srk_der_T, p_ref, [tp_in.get_variable("T")], p1, n_steps=30)
    integrator_overall = RK45(rk_overall_der_T, p_ref, [tp_in.get_variable("T")], p1)

else:

    integrator_srk = SimpleIntegrator(rk_srk_der, p_ref, [1/tp_in.get_variable("rho")], p1, n_steps=30)
    integrator_overall = RK45(rk_overall_der, p_ref, [tp_in.get_variable("rho")], p1)

while integrator_srk.status == 'running':

    integrator_srk.step()

    points_srk["p"].append(curr_state.p)
    points_srk["rho"].append(1 / curr_state.v)
    points_srk["T"].append(curr_state.t)
    points_srk["s"].append(curr_state.s - ref_state.s)
    points_srk["h"].append(curr_state.h - ref_state.h)

    points_srk["cv"].append(curr_state.cv)
    points_srk["R"].append(curr_state.r)

while integrator_overall.status == 'running':

    integrator_overall.step()

    points_overall["p"].append(tp_curr.get_variable("P"))
    points_overall["rho"].append(tp_curr.get_variable("rho"))
    points_overall["T"].append(tp_curr.get_variable("T"))
    points_overall["s"].append(tp_curr.get_variable("s") - tp_ref.get_variable("s"))
    points_overall["h"].append(tp_curr.get_variable("h") - tp_ref.get_variable("h"))

    points_overall["cv"].append(tp_curr.get_variable("cv"))
    points_overall["R"].append(tp_curr.get_variable("cp") - tp_curr.get_variable("cv"))

output_overall = integrator_overall.y
output_srk = integrator_srk.y

err_abs_overall = (output_overall[0] - output_srk[0])
err_rel_overall = err_abs_overall / output_srk[0]
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

df_simple = pd.DataFrame.from_dict(points_srk)
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
        label="EoS based", color=colors[i], marker="+", zorder=200,
        s=150

    )

    ax.set_xlabel("Pressure [bar]", fontsize=17)
    ax.set_ylabel(labels[i], fontsize=17, labelpad=12)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    #ax.set_xscale("log")
    lines.append(new_line[0])
    overall_scatters.append(new_ovr_scatter)


lines.append(overall_scatters[0])
labs = [l.get_label() for l in lines]
ax.legend(lines, labs, fontsize=15)


# %%-------------------------------------   SHOW PLOT                           -------------------------------------> #
plt.tight_layout(pad=5)
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "2 - Integration Check", "output", "direct comparison.png"))
