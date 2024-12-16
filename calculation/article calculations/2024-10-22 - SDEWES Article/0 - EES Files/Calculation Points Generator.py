# %%------------          Importing Libraries            -----------------------------------------------------------> #
import numpy as np

# %%------------   Generating Calculation Points         -----------------------------------------------------------> #
n_points = 8

depth_list = np.round(np.linspace(1000, 5000, n_points), 1)
grad_list = np.round(np.linspace(30, 100, n_points), 1) / 1000
h_rel = np.logspace(-1, -0.01, n_points)
alpha_hp = np.linspace(1.05, 2, n_points)

depth_list, grad_list, dt_HE, alpha_hp = np.meshgrid(depth_list, grad_list, h_rel, alpha_hp, indexing='ij')

depth_list = np.ravel(depth_list)
grad_list = np.ravel(grad_list)
dt_HE = np.ravel(dt_HE)
alpha_hp = np.ravel(alpha_hp)

table = np.stack((depth_list, grad_list, dt_HE, alpha_hp)).T
