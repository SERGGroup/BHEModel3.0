# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import numpy as np
import os

n_points = 15

depth_list = np.round(np.linspace(1000, 5000, n_points), 1)
grad_list = np.round(np.linspace(30, 100, n_points), 1) / 1000
dt_HE = np.linspace(2, 20, n_points)

depth_list, grad_list, dt_HE = np.meshgrid(depth_list, grad_list, dt_HE, indexing='ij')

depth_list = np.ravel(depth_list)
grad_list = np.ravel(grad_list)
dt_HE = np.ravel(dt_HE)

table = np.stack((depth_list, grad_list, dt_HE)).T