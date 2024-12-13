# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
import copy

from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.spatial.distance import cdist
from owslib.wms import WebMapService
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import requests
import io

from tqdm import tqdm

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
plt.figure(figsize=(10, 10), dpi=500)
plt.imshow(image_array)
plt.title("Temperature @ {}m".format(depth))
plt.axis('off')  # Turn off axis labels
plt.show()


# %%------------   IDENTIFY UNIQUE COLORS AND SORT BY APPEARANCE        ----------------------------------------------->
rgb_triplets = image_array.reshape(-1, image_array.shape[-1])
unique_rgb = np.unique(rgb_triplets, axis=0)
pixels = image_array.reshape(-1, 3)

pixel_counts = []
pbar = tqdm(total=len(unique_rgb))
for color in unique_rgb:

    matching_pixels = np.all(image_array == color, axis=-1)
    pixel_counts.append(np.sum(matching_pixels))
    pbar.update(1)

pbar.close()

pixel_counts = np.array(pixel_counts)
counts_sorted = pixel_counts[np.argsort(pixel_counts)][::-1]
color_sorted = unique_rgb[np.argsort(pixel_counts)][::-1]


# %%------------   GROUP SIMILAR COLORS                     ----------------------------------------------------------->
get_first = 25
tolerance = 30
not_similar_colors = [unique_rgb[0]]
not_similar_counts = [0]

for i, color in enumerate(color_sorted[0:get_first+1]):

    distances = cdist([color], not_similar_colors)
    if np.min(distances) > tolerance:  # Add only if it exceeds the tolerance
        not_similar_colors.append(color)
        not_similar_counts.append(counts_sorted[i])

    else:
        similar_j = np.where(distances[0] == np.min(distances))[0][0]
        not_similar_counts[similar_j] = not_similar_counts[similar_j] + counts_sorted[i]

not_similar_counts = np.array(not_similar_counts)
not_similar_colors = np.array(not_similar_colors)
counts_final = not_similar_counts[np.argsort(not_similar_counts)][::-1]
final_colors = not_similar_colors[np.argsort(not_similar_counts)]


# %%------------   TRY IMPORT MAP                         ----------------------------------------------------------->
cmap_colors = copy.deepcopy(final_colors)[1:]
cmap_colors[-1, :] = cmap_colors[-2, :]
cmap = ListedColormap(cmap_colors / 255.0)

data_range = np.array(range(len(cmap_colors)))
norm = BoundaryNorm(data_range, cmap.N)
fig, ax = plt.subplots(figsize=(6, 1))
fig.subplots_adjust(bottom=0.5)

# Add the colorbar
cb = plt.colorbar(
    plt.cm.ScalarMappable(cmap=cmap, norm=norm),
    cax=ax,
    orientation='horizontal',
    ticks=data_range,
)
cb.set_label('Custom Colorbar')
plt.show()
