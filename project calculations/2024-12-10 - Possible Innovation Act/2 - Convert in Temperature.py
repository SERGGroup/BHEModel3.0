# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.spatial.distance import cdist
from owslib.wms import WebMapService
import matplotlib.pyplot as plt
from tqdm import tqdm
from PIL import Image
import numpy as np
import requests
import io

color_code = [

[[0, 169, 230], [20, 30]],
[[115, 178, 255], [30, 40]],
[[190, 232, 255], [40, 50]],
[[137, 205, 102], [50, 60]],
[[165, 245, 122], [60, 70]],
[[215, 215, 158], [70, 80]],
[[215, 176, 158], [80, 90]],
[[255, 255, 115], [90, 100]],

[[230, 230, 0], [50, 100]],
[[255, 170, 0], [100, 150]],

[[211, 255, 190], [50, 75]],
[[230, 152, 0], [75, 100]],

[[255, 190, 232], [150, 200]],
[[245, 122, 182], [200, 250]],
[[255, 0, 197], [250, 300]],
[[230, 0, 0], [300, np.nan]],

]


# %%------------   TRY IMPORT MAP                         -----------------------------------------------------------> #
wms_url = "https://repo2.igg.cnr.it/geoserver/IGG/wms?service=WMS&version=1.1.0&request=GetMap"
wms = WebMapService(wms_url)
print("Available Layers:", list(wms.contents.keys()))

depth = 1000
layer = "IGG:area_temp_{}".format(depth)  # Replace with the layer name you need
bbox = (6.6273, 35.4922, 18.7845, 47.0921)  # Bounding box for Italy
width, height = 2*1024, 2*1024  # Image resolution


# %%------------   TRY IMPORT MAP                         -----------------------------------------------------------> #
# GetMap request parameters
params = {
    'service': 'WMS',
    'version': '1.1.0',
    'request': 'GetMap',
    'layers': layer,
    'bbox': ','.join(map(str, bbox)),  # Convert bbox to string format
    'width': width,
    'height': height,
    'srs': 'EPSG:4326',  # Default CRS is EPSG:4326
    'format': 'image/png',  # Image format
    'styles': '',  # Leave blank if no style is required
}

# Get the map image
getmap_url = wms.getOperationByName('GetMap').methods[0][('url')]


# %%------------   TRY IMPORT MAP                         ----------------------------------------------------------->
response = requests.get(getmap_url, params=params)

if response.status_code == 200:
    image = Image.open(io.BytesIO(response.content))
    print(f"Image retrieved")

else:
    print(f"Error retrieving map: {response.status_code}, {response.text}")

# %%------------   TRY IMPORT MAP                           ----------------------------------------------------------->
image_array = np.array(image)
plt.figure(figsize=(10, 10), dpi=300)
plt.imshow(image_array)
plt.title("Temperature @ {}m".format(depth))
plt.axis('off')  # Turn off axis labels
plt.show()


# %%------------   TRY IMPORT MAP                           ----------------------------------------------------------->
pixels = image_array.reshape(-1, 3)

temperature_array = np.empty(len(pixels))
temperature_array[:] = np.nan

for color in color_code:

    matching_pixels = np.all(pixels == color[0], axis=-1)
    temperature_array[matching_pixels] = np.nanmean(color[1])

temperature_array = temperature_array.reshape(image_array.shape[0:2])

min_t = np.nanmin(temperature_array)
max_t = np.nanmax(temperature_array)
temperature_img = (temperature_array - min_t) / (max_t - min_t)

plt.figure(figsize=(10, 10), dpi=500)
plt.imshow(temperature_img)
plt.title("Temperature @ {}m".format(depth))
plt.axis('off')  # Turn off axis labels
plt.show()