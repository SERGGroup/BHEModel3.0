# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes.geothermal_system.base_bhe import isobaricIntegral, ThermodynamicPoint
from main_classes.constant import ARTICLE_CALCULATION_DIR
from openpyxl import load_workbook
import numpy as np
import os


# %%------------   GENERATE CALCULATION POINTS            -----------------------------------------------------------> #
n_points = 15

depth_list = np.round(np.linspace(1000, 5000, n_points), 1)
grad_list = np.round(np.linspace(30, 100, n_points), 1) / 1000
dt_HE = np.linspace(2, 20, n_points)

depth_list, grad_list, dt_HE = np.meshgrid(depth_list, grad_list, dt_HE, indexing='ij')

depth_list = np.ravel(depth_list)
grad_list = np.ravel(grad_list)
dt_HE = np.ravel(dt_HE)

table = np.stack((depth_list, grad_list, dt_HE)).T


# %%------------   IMPORT RESULTS                         -----------------------------------------------------------> #
CURRENT_DIR = os.path.join(ARTICLE_CALCULATION_DIR, "2024-10-22 - SDEWES Article")
file_path = os.path.join(CURRENT_DIR, "0 - Results", "Water (Base).xlsx")

workbook = load_workbook(filename=file_path)

sheet = workbook['Results']
n_rows = sheet.max_row - 1

headers = [cell.value for cell in sheet[1]]
data_dict = {header: np.empty(n_rows) for header in headers}

for row in sheet.iter_rows(min_row=2, values_only=True):
    for key, cell_value in zip(headers, row):
        data_dict[key].append(cell_value)

workbook.close()

for key in headers:
    data_dict[key] = np.array(data_dict[key])


# %%------------   EVALUATE OVERALL LCOH                  -----------------------------------------------------------> #
depth_list = np.unique(data_dict['depth'])
grad_list = np.unique(data_dict['grad'])
dt_HE_list = np.unique(data_dict['dT_HE'])

curr_depth = depth_list[5]
curr_grad = grad_list[5]
curr_indices = np.where(np.logical_and(data_dict['depth'] == curr_depth, data_dict['grad'] == curr_grad))
curr_data = {header: [] for header in headers}

for key in headers:
    curr_data[key] = data_dict[key][curr_indices]

curr_data.update({"LCOH": np.empty(curr_data[headers[0]].shape)})
