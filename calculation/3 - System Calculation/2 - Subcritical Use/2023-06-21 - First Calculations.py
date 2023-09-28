# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import RKEOS, evaluate_system, evaluate_surface, calculate_flash
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
fluid = RKEOS()

depth = 3000.
t_res = 500

in_state = calculate_flash(fluid, "TQ", 280, 0)
states = evaluate_system(fluid, in_state, depth, t_res)
res = evaluate_surface(fluid, states)
