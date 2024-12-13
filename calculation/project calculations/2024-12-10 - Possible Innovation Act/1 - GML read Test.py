# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.interpolate import interp1d
from shapely.geometry import Point
import matplotlib.pyplot as plt
import geopandas as gpd
from tqdm import tqdm
import numpy as np
import os

if os.name == "nt":
    from main_classes.constant import PROJECT_CALCULATION_DIR

else:
    PROJECT_CALCULATION_DIR = "/Users/PietroUngar/PycharmProjects/BHEModel3.0/calculation/project calculations"

CURRENT_DIR = os.path.join(PROJECT_CALCULATION_DIR, "2024-12-10 - Possible Innovation Act")
GML_TEMP_FOLDER = os.path.join(CURRENT_DIR, "0 - Data", "GML Temperature Layers")


# %%------------   IMPORT GML FILES                       -----------------------------------------------------------> #
gdf = {}
depths = [1000, 2000, 3000]
depths_with_zero = [0, 1000, 2000, 3000]

for depth in [1000, 2000, 3000]:

    file_name = f'area_temp_{depth}.gml'
    # Read the GML file
    gdf.update({depth: gpd.read_file(os.path.join(GML_TEMP_FOLDER, file_name))})


# %%------------   EVALUATE CONDITION IN 1 POINT          -----------------------------------------------------------> #
florence = Point(11.2558, 43.7696)  # longitude, latitude
depth = 2000
mask = gdf[depth].contains(florence)
data_at_florence = gdf[depth][mask]

print(np.mean(data_at_florence.TEMP))


# %%------------   EVALUATE CONDITION WITH IDW            -----------------------------------------------------------> #
depth = 1500
n = 5


def idw_interpolation(longitude, latitude, depth_well, power=2):

    point = Point(longitude, latitude)
    inside_check = gdf[depths[0]].contains(point)

    if inside_check.any():

        result_list = [10]
        for depth in depths:

            dist = gdf[depth].boundary.apply(lambda x: point.distance(x))
            values = gdf[depth].values[:, 2]
            min_dist = np.min(dist)

            if min_dist == 0:
                result_list.append(values[np.where[dist == min_dist]])

            else:
                w = (1 / dist.replace(0, np.nan)) ** power
                sums = gdf[depth].values[:, 2] * w
                result_list.append(sums.sum() / w.sum())

        linear_interp = interp1d(depths_with_zero, result_list, kind='linear', fill_value="extrapolate")
        return linear_interp(depth_well)

    else:

        return np.nan


print(idw_interpolation(florence.x, florence.y, depth_well=1500, power=n))


# %%------------   EVALUATE MAP OF ITALY                  -----------------------------------------------------------> #
dz_well = 1500
n = 5

lat_min, lat_max = 35.5, 47.1
lon_min, lon_max = 6.6, 18.5
n_points = 400j
grid_lon, grid_lat = np.mgrid[lon_min:lon_max:n_points, lat_min:lat_max:n_points]

pbar = tqdm(total=grid_lon.shape[0] * grid_lon.shape[1], desc="Calculating map")
grid_z = np.zeros_like(grid_lon)
for i in range(grid_lon.shape[0]):
    for j in range(grid_lon.shape[1]):
        grid_z[i, j] = idw_interpolation(grid_lon[i, j], grid_lat[i, j], depth_well=dz_well, power=n)
        pbar.update(1)

pbar.close()


# %%------------   PLOT MAP OF ITALY                      -----------------------------------------------------------> #
plt.figure(figsize=(10, 10))
plt.pcolormesh(grid_lon, grid_lat, grid_z, shading='auto', cmap='viridis')
plt.colorbar(label='Temperature [°C]')
plt.title(f'Temperature @ depth = {dz_well}m')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.savefig(os.path.join(CURRENT_DIR, "0 - Output", f"Temperature Result - depth = {dz_well}.png"), dpi=300)
plt.show()
