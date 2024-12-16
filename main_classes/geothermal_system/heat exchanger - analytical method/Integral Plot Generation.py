# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import os


# %%------------   IMPORT DATA                            -----------------------------------------------------------> #
proj_folder = "/Users/PietroUngar/PycharmProjects/BHEModel3.0"
curr_folder = "main_classes/geothermal_system/heat exchanger - analytical method"
wb = openpyxl.load_workbook(os.path.join(proj_folder, curr_folder, 'Integral Description.xlsx'), data_only=True)
sheet = wb['Integration']

# Read data from a specific range
x_data = []
y1_data = []
y2_data = []

for row in sheet['F5:J1005']:
    x_data.append(row[0].value)
    y1_data.append(row[1].value)
    y2_data.append(row[4].value)

x_data = np.array(x_data) / 100
y1_data = np.array(y1_data)
y2_data = np.array(y2_data)


# %%------------   PLOT INTEGRAL                          -----------------------------------------------------------> #
# Creating the plot with dual y-axes
fig, ax1 = plt.subplots(figsize=(6, 5))
plt.xlim([0, 1])

# Plotting the first line
color = 'tab:blue'
ax1.set_ylim([0, 50])
label = '$\Delta T$ [°C]'
ax1.set_xlabel('$\Delta h_{\%}$')
ax1.tick_params(axis='x', labelsize=12)
ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
line1, = ax1.plot(x_data, y1_data, color=color, label=label)

# Creating a second y-axis for the second line
ax2 = ax1.twinx()
color = 'tab:red'
label = '$\\int_{0}^{\\Delta h_{\%}}\\frac{dh\'}{\Delta T(h\')}$'

#ax2.set_yscale('log')
ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
line2, = ax2.plot(x_data, y2_data, color=color, label=label)
ax2.set_ylim([0, np.nanmax(y2_data[:-1])])

# Title and grid
lines = [line1, line2]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper right', fontsize=12)

plt.title('Numerical Integration of $\Delta T(h\')$', fontsize=14)
plt.tight_layout(pad=2)
plt.savefig(os.path.join(proj_folder, curr_folder, "numerical_integration.png"), dpi=300)

plt.show()
