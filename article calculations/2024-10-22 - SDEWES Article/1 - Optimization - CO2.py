# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import numpy as np
import os

n_points = 18

depth_list = np.round(np.linspace(1000, 5000, n_points), 1)
grad_list = np.round(np.linspace(30, 100, n_points), 1) / 1000
h_rel_list = 1 - np.logspace(-3, -0.05, n_points + 3)


depth, grad, h_rel = np.meshgrid(depth_list, grad_list, h_rel_list, indexing='ij')

depth = np.ravel(depth)
grad = np.ravel(grad)
h_rel = np.ravel(h_rel)


table = np.stack((depth, grad, h_rel)).T
