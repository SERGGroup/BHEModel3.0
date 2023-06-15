# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.geothermal_system import evaluate_system
from main_classes.subclasses.redlich_kwong import RKEOS
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
fluid = RKEOS(

    p_crit=1e7,
    t_crit=350,
    cp_ideal=850,
    m_molar=0.044

)

depth = 3000.
grad_rel = 400
grad = grad_rel * fluid.t_crit / depth

t_rels = 1 - np.logspace(0, -5, num=100) * 0.5
sat_points = np.zeros((2*len(t_rels), 2))
betas = np.zeros(t_rels.shape)
dps = np.zeros(t_rels.shape)
states_list = list()

p_mult = 1
label = "p = p_sat"
pbar = tqdm(desc="Calculating Points", total=len(t_rels))

for i in range(len(t_rels)):

    t_inj = t_rels[i] * fluid.t_crit
    t_rock = t_inj + grad * depth / 1000

    liq_state, vap_state = fluid.get_sat_state(t=t_inj, which="both")

    if p_mult == 1:

        input_state = fluid.get_state(t=liq_state.t, v=liq_state.v * 0.999)

    elif p_mult > 1:

        label = "p = p_sat * {}".format(p_mult)
        input_state = fluid.get_state(p=liq_state.p*p_mult, t=liq_state.t)

    else:

        label = "p = p_sat * {}".format(p_mult)
        input_state = fluid.get_state(p=liq_state.p * p_mult, t=liq_state.t)

    sat_points[i, :] = np.array([liq_state.p, liq_state.v])
    sat_points[-(i + 1), :] = np.array([vap_state.p, vap_state.v])

    sys_states = evaluate_system(fluid, input_state, depth, t_rock)
    dps[i] = (sys_states[-1].p - sys_states[0].p) / fluid.p_crit
    betas[i] = sys_states[-1].p / sys_states[0].p
    states_list.append(sys_states)
    pbar.update(1)

pbar.close()
plt.plot(t_rels, betas, label=label)


# %%-------------------------------------   PLOT PROFILES                       -------------------------------------> #
plt.xlabel("$T_{inj}$ / $T_{crit}$ [-]")
plt.ylabel("$dp_{sys}$ / $p_{crit}$ [-]")
plt.legend()
plt.show()


# %%-------------------------------------   PLOT PROFILES                       -------------------------------------> #
plt.plot(sat_points[:, 1], sat_points[:, 0], color="black")

counter = 0
plt_every = 20

for sys_states in states_list:

    p_list = list()
    v_list = list()

    for state in sys_states[0:]:

        p_list.append(state.p)
        v_list.append(state.v)

    if np.mod(counter, plt_every) == 0:
        plt.plot(v_list, p_list)

    counter += 1

plt.xscale("log")
plt.yscale("log")

plt.show()
