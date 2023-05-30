# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.rk_fluid_classes import RKFluid, evaluate_system
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
p_mult = 1
grad_rel = 600

depths = np.logspace(3, 4, num=10)
t_crits = [250, 450, 650, 850]
t_rels = 1 - np.logspace(0, -5, num=100) * 0.5
betas = np.zeros((len(t_rels), len(depths), len(t_crits)))

pbar = tqdm(desc="Calculating Points", total=len(t_rels) * len(depths) * len(t_crits))

for k in range(len(t_crits)):

    fluid = RKFluid(

        p_crit=1e7,
        t_crit=t_crits[k],
        cp_ideal=850,
        m_molar=0.044

    )

    for j in range(len(depths)):

        depth = depths[j]
        grad = grad_rel * fluid.t_crit / depth

        for i in range(len(t_rels)):

            t_inj = t_rels[i] * fluid.t_crit
            t_rock = t_inj + grad * depth / 1000

            liq_state = fluid.get_sat_state(t=t_inj, liquid=True)
            vap_state = fluid.get_sat_state(t=t_inj, liquid=False)

            if p_mult == 1:

                input_state = fluid.get_state(t=liq_state.t, v=liq_state.v * 0.99)

            elif p_mult > 1:

                input_state = fluid.get_state(p=liq_state.p*p_mult, t=liq_state.t)

            else:

                input_state = fluid.get_state(p=liq_state.p * p_mult, t=liq_state.t)

            sys_states = evaluate_system(fluid, input_state, depth, t_rock)
            betas[i, j, k] = sys_states[-1].p / sys_states[0].p
            pbar.update(1)

        plt.plot(t_rels, betas[:, j, k], label="depth={}".format(int(np.ceil(depth))))

    plt.xlabel("$T_{inj}$ / $T_{crit}$ [-]")
    plt.ylabel("$beta$ [-]")
    plt.title("T_crit = {}".format(t_crits[k] - 273.15))
    plt.legend()
    plt.show()

pbar.close()


# %%-------------------------------------   PLOT RATIOS                         -------------------------------------> #
corr_d = 0.25
corr_t = -0.25
ratios = np.zeros((len(t_rels), len(depths), len(t_crits)))

for k in range(len(t_crits)):

    for j in range(len(depths)):

        ratios[:, j, k] = betas[:, j, k] / betas[:, 0, 0]

    d_mod = (depths[1] / depths[0]) ** corr_d
    t_mod = (t_crits[k] / t_crits[0]) ** corr_t

    plt.plot(t_rels, ratios[:, 1, k])
    plt.plot(t_rels, np.ones(t_rels.shape) * d_mod * t_mod, "--")

plt.show()


# %%-------------------------------------   PLOT CORRECTED                      -------------------------------------> #
for k in range(len(t_crits)):

    t_crit = t_crits[k]

    for j in range(len(depths)):

        depth = depths[j]
        plt.plot(t_rels, betas[:, j, k]/(depth**corr_d * t_crit**corr_t), label="depth={}".format(int(np.ceil(depth))))

plt.xlabel("$T_{inj}$ / $T_{crit}$ [-]")
plt.ylabel("$beta_{corr}$ [-]")
#plt.legend()
plt.show()
