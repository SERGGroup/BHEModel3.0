# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import RKEOS, evaluate_surface, calculate_flash
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
fluid = RKEOS()

depth = 3000.
t_res = 460

in_state = calculate_flash(fluid, "TQ", 285, 0)
res = evaluate_surface(fluid, in_state, depth, t_res)
