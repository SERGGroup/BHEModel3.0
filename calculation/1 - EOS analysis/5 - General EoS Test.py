# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.cubic_eos import CubicEOS
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   FUNCTION DEFINITION                 -------------------------------------> #
fluid = CubicEOS(

    p_crit=1e6,
    t_crit=1e5,
    cp_ideal=850,
    m_molar=0.044

)

p_sat, z_l, z_v = fluid.p_sat(0.6*fluid.t_crit)
t_sat, z_l_2, z_v_2 = fluid.t_sat(p_sat)
check = t_sat / fluid.t_crit


# %%-------------------------------------   CALCULATE EOS                       -------------------------------------> #
n_points = 500
ps = np.logspace(-15, 0, n_points)
ts = np.array([0.2, 0.5, 0.7, 0.8, 0.9, 0.99])
p_rels, t_rels = np.meshgrid(ps, ts, indexing='ij')

v_liqs = np.empty((n_points, len(ts)))
z_liqs = np.empty((n_points, len(ts)))
fug_liqs = np.empty((n_points, len(ts)))
ddp_liqs = np.empty((n_points, len(ts)))

v_vaps = np.empty((n_points, len(ts)))
z_vaps = np.empty((n_points, len(ts)))
fug_vaps = np.empty((n_points, len(ts)))
ddp_vaps = np.empty((n_points, len(ts)))


# %%-------------------------------------   CALCULATE REFPROP                   -------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=n_points * len(ts))

for k in range(len(ts)):

    for j in range(n_points):

        t_curr = t_rels[j, k]*fluid.t_crit
        p_curr = p_rels[j, k]*fluid.p_crit
        z_l, z_v = fluid.z(t=t_curr, p=p_curr)

        z_liqs[j, k] = z_l
        v_liqs[j, k] = z_l * fluid.r_spc * t_curr / p_curr
        fug_liqs[j, k] = fluid.fug(t_curr, p_curr, z_l)
        alpha = fluid.a(t_curr) / (v_liqs[j, k] - fluid.b) - fluid.a(fluid.t_crit) / (fluid.v_crit - fluid.b)
        ddp_liqs[j, k] = alpha

        z_vaps[j, k] = z_v
        v_vaps[j, k] = z_v * fluid.r_spc * t_curr / p_curr
        fug_vaps[j, k] = fluid.fug(t_curr, p_curr, z_v)
        alpha = fluid.a(t_curr) / (v_vaps[j, k] - fluid.b) - fluid.a(fluid.t_crit) / (fluid.v_crit - fluid.b)
        ddp_vaps[j, k] = alpha

        pbar.update(1)

pbar.close()


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
i = 0
plt.plot(p_rels[:, i], fug_liqs[:, i])
plt.plot(p_rels[:, i], fug_vaps[:, i])

plt.title("fug")
plt.xscale("log")
#plt.yscale("log")
#plt.ylim((-25,0))
plt.show()

plt.plot(p_rels[:, i], ddp_liqs[:, i])
plt.plot(p_rels[:, i], ddp_vaps[:, i])

plt.title("z")
plt.xscale("log")
#plt.yscale("log")
#plt.ylim((-2,2))
plt.show()
