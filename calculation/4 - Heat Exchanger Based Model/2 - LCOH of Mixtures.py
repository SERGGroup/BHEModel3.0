# %%-------------------------------------   IMPORT MODULES                          ---------------------------------> #
from main_classes import BaseBHE, economicEvaluator
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


# %%-------------------------------------   INIT CALCULATION                        ---------------------------------> #
m_dot = 1
t_in = 10           # [°C]
depth = 500         # [m]
gradient = 0.07     # [°C/m]
mix_conc = [0, 0.25, 0.5, 1]
bhe_mix_list = []
other_fluid = "Ethane"
p_list = list()

for conc in mix_conc:

    bhe_in_mix = ThermodynamicPoint(["CarbonDioxide", other_fluid], [1-conc, conc], unit_system="MASS BASE SI")
    bhe_in_mix.set_variable("T", t_in + 273.15)
    bhe_in_mix.set_variable("Q", 0)
    p_list.append(bhe_in_mix.get_variable("P"))

    mix_bhe = BaseBHE(bhe_in_mix)
    mix_bhe.geom.depth = depth
    mix_bhe.res_prop.grad = gradient

    mix_bhe.set_HX_condition()
    bhe_mix_list.append(mix_bhe)

# %%-------------------------------------   ECONOMIC EVALUATOR                      ---------------------------------> #
evaluators = []

for mix_bhe in bhe_mix_list:
    evaluators.append(economicEvaluator(mix_bhe))

l_m_dot_kms = np.logspace(-2, 2, 50)
LCOH = list()
LCOex = list()

DT_list = list()
beta_list = list()

for l_m_dot_km in l_m_dot_kms:

    l_km = l_m_dot_km * m_dot
    sub_LCOH = list()
    sub_LCOex = list()

    for evaluator in evaluators:

        LCOx = evaluator.evaluate_LCOx(["w_out", "ex_out"], len_hor=l_km*1e3, m_dot=m_dot)

        sub_LCOH.append(LCOx[0])
        sub_LCOex.append(LCOx[1])

    LCOH.append(sub_LCOH)
    LCOex.append(sub_LCOex)

LCOH = np.array(LCOH)*1e4
LCOex = np.array(LCOex)*1e4


# %%-------------------------------------   PLOT ECONOMIC RESULTS                   ---------------------------------> #
lines = ["--", "-.", ":"]

line = plt.plot(l_m_dot_kms, LCOH[:, 0], label="LCOH")
for i in range(len(evaluators)-1):

    plt.plot(

        l_m_dot_kms, LCOH[:, i+1], lines[i],
        color=line[0].get_color(),
        label="LCOH - {}% {}".format(
            int(mix_conc[i+1]*100),
            other_fluid
        )

    )

line = plt.plot(l_m_dot_kms, LCOex[:, 0], label="LCOex")
for i in range(len(evaluators)-1):

    plt.plot(

        l_m_dot_kms, LCOex[:, i+1], lines[i],
        color=line[0].get_color(),
        label="LCOex - {}% {}".format(

            int(mix_conc[i+1]*100),
            other_fluid

        )

    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$l_{rel}$ (km*s/kg)")
plt.ylabel("LCOx (c€/kW)")
plt.legend()
plt.show()


# %%-------------------------------------   EVALUATE THERMODYNAMICS                 ---------------------------------> #
DT_list = list()
beta_list = list()
dh_prec_list = list()

pbar = tqdm(desc="Calculating Points", total=len(l_m_dot_kms)*len(bhe_mix_list))

for l_m_dot_km in l_m_dot_kms:

    DT_sub_list = list()
    beta_sub_list = list()
    dh_prec_sub_list = list()

    l_hor = l_m_dot_km * m_dot

    i = 0
    for mix_bhe in bhe_mix_list:

        mix_bhe.geom.l_horiz = l_hor
        dh_percs, out_points = mix_bhe.evaluate_HXG([10 * 3.154e+7], m_dot=m_dot)

        DT_sub_list.append(out_points[0].get_variable("T") - (t_in+273.15))
        beta_sub_list.append(out_points[0].get_variable("P") / p_list[i])
        i += 1

        dh_prec_sub_list.append(dh_percs[0])
        pbar.update(1)

    DT_list.append(DT_sub_list)
    beta_list.append(beta_sub_list)
    dh_prec_list.append(dh_prec_sub_list)

pbar.close()
values = [np.array(DT_list), np.array(beta_list), np.array(dh_prec_list)]


# %%-------------------------------------   PLOT OUTPUT PRESSURE AND TEMPERATURE    ---------------------------------> #
lines = ["--", "-.", ":"]

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()

line = ax.plot(l_m_dot_kms, values[0][:, 0], label="$beta$")
for i in range(len(evaluators)-1):

    ax.plot(

        l_m_dot_kms, values[0][:, i+1], lines[i],
        color=line[0].get_color(),
        label="$beta$ - {}% {}".format(
            int(mix_conc[i+1]*100),
            other_fluid
        )

    )

line = ax2.plot(l_m_dot_kms, values[2][:, 0], color="orange", label="$\Delta T$")
for i in range(len(evaluators)-1):

    ax2.plot(

        l_m_dot_kms, values[2][:, i+1], lines[i],
        color=line[0].get_color(),
        label="$\Delta T$ - {}% {}".format(

            int(mix_conc[i+1]*100),
            other_fluid

        )

    )

plt.xscale("log")
# plt.yscale("log")
plt.xlabel("$l_{rel}$ (km*s/kg)")
# plt.ylabel("LCOx (c€/kW)")
# plt.legend()
plt.show()


# %%
plt.plot(mix_conc, p_list)
plt.show()