# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import numpy as np
import os

n_points = 5

depth_list = np.round(np.linspace(1000, 5000, n_points), 1)
grad_list = np.round(np.linspace(30, 100, n_points), 1) / 1000
h_rel = 1 - np.logspace(-3, -0.25, n_points)
t_sat = np.linspace(10, 25, n_points)

depth_list, grad_list, h_rel, t_sat = np.meshgrid(depth_list, grad_list, h_rel, t_sat, indexing='ij')

depth_list = np.ravel(depth_list)
grad_list = np.ravel(grad_list)
h_rel = np.ravel(h_rel)
t_sat = np.ravel(t_sat)

table = np.stack((depth_list, grad_list, h_rel, t_sat)).T
