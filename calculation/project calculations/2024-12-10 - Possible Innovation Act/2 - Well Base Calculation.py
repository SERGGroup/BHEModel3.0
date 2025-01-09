# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes import BaseBHE, ReservoirProperties, BHEGeometry, baseEconomicEvaluator
from REFPROPConnector import ThermodynamicPoint
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d
from functools import total_ordering
from matplotlib import pyplot as plt
from shapely.geometry import Point
from pykml import parser
import geopandas as gpd
from tqdm import tqdm
import numpy as np
import lxml.html
import requests
import os


if os.name == "nt":
    from main_classes.constant import PROJECT_CALCULATION_DIR

else:
    PROJECT_CALCULATION_DIR = "/Users/PietroUngar/PycharmProjects/BHEModel3.0/calculation/project calculations"

CURRENT_DIR = os.path.join(PROJECT_CALCULATION_DIR, "2024-12-10 - Possible Innovation Act")
DATA_FOLDER = os.path.join(CURRENT_DIR, "0 - Data")
GML_TEMP_FOLDER = os.path.join(DATA_FOLDER, "GML Temperature Layers")
kml_file = os.path.join(DATA_FOLDER, "pozzi_consultabili.kml")


# %%------------   DEFINE CLASSES                         -----------------------------------------------------------> #
class GMLInterpolator:

    __depths = [1000, 2000, 3000]
    __depths_with_zero = [0, 1000, 2000, 3000]
    __t_surf = 10

    def __init__(self, gml_folder, power=5):

        self.__power = power
        self.__read_gml_files(gml_folder)

    def __read_gml_files(self, gml_folder):

        self.__gdf = {}

        for depth in self.__depths:
            file_name = f'area_temp_{depth}.gml'
            self.__gdf.update({depth: gpd.read_file(os.path.join(gml_folder, file_name))})

    def get(self, longitude, latitude, depth_well):

        interp = self.get_interpolator(longitude, latitude)
        if interp is not None:

            return interp(depth_well)

        else:

            return np.nan

    def get_interpolator(self, longitude, latitude):

        point = Point(longitude, latitude)
        inside_check = self.__gdf[self.__depths[0]].contains(point)

        if inside_check.any():

            result_list = [self.__t_surf]
            for depth in self.__depths:

                dist = self.__gdf[depth].boundary.apply(lambda x: point.distance(x))
                values = self.__gdf[depth].values[:, 2]
                min_dist = np.min(dist)

                if min_dist == 0:
                    result_list.append(values[np.where[dist == min_dist]])

                else:
                    w = (1 / dist.replace(0, np.nan)) ** self.__power
                    sums = self.__gdf[depth].values[:, 2] * w
                    result_list.append(sums.sum() / w.sum())

            linear_interp = interp1d(self.__depths_with_zero, result_list, kind='linear', fill_value="extrapolate")
            return linear_interp

        else:

            return None

    def is_inside(self, point):

        for depth in self.__depths:

            if not self.__gdf[depth].contains(point):
                return False

        return True

@total_ordering
class WellOptimizer:

    m_dot = 10    # [kg/s]
    time = 10 * (365 * 24 * 60 * 60)  # [years] to [s]
    economic_evaluator = baseEconomicEvaluator()
    economic_evaluator.Le = 20

    def __init__(self, well_placemark, temp_gdf, temp_iterp, bhe_in):

        self.__placemark = well_placemark
        self.__temp_iterp = temp_iterp
        self.__gdf = temp_gdf

        self.__init_attributes()
        self.__set_main_information()
        self.__analyze_extended_attr()
        self.__retrieve_information_online()

        if self.can_be_optimized:

            self.__init_thermo_well(bhe_in)

    @property
    def can_be_optimized(self):
        return self.t_down is not None

    @property
    def w_dot_opt(self):

        """return w_dot in W"""

        return self.thermo_well.integrator.dh_max * self.m_dot * self.h_rel_opt

    def __init_attributes(self):

        self.link = None
        self.response_content = None
        self.year = None
        self.operator = None
        self.goal = None
        self.depth = None
        self.t_down = None
        self.outcome = None
        self.placement = None

        self.h_rel_opt = None
        self.lcoh_min = None
        self.l_horiz_opt = None

    def __set_main_information(self):

        self.name = self.__placemark.name.text if hasattr(self.__placemark, 'name') else "Unknown"
        i_gdf = np.where(self.__gdf.values[:, 0] == self.name)
        self.position = self.__gdf.values[i_gdf, 2][0][0]
        self.interpolator = self.__temp_iterp.get_interpolator(self.position.x, self.position.y)

    def __analyze_extended_attr(self):

        if hasattr(self.__placemark, 'ExtendedData'):

            for curr_data in self.__placemark.ExtendedData.SchemaData.findall('{http://www.opengis.net/kml/2.2}SimpleData'):

                if curr_data.attrib['name'] == 'LINK':
                    tree = lxml.html.fromstring(curr_data.text)
                    self.link = tree.xpath('//a/@href')[0]

                if curr_data.attrib['name'] == 'ANNO':
                    self.year = curr_data.text

                if curr_data.attrib['name'] == 'OPERATORE':
                    self.operator = curr_data.text

    def __retrieve_information_online(self):

        if self.link is not None:

            response = requests.get(self.link)

            if response.status_code == 200:

                doc = lxml.html.fromstring(response.text)
                table = doc.xpath("//span[contains(text(), 'Dati generali del pozzo')]/following-sibling::table[1]")[0]

                # Extract table rows
                self.response_content = []
                for row in table.xpath(".//tr"):

                    cells = row.xpath(".//td/text()")
                    self.response_content.append(cells)

                    if cells[0] == 'Scopo':
                        self.goal = cells[1]

                    if cells[0] == 'Esito':
                        self.outcome = cells[1]

                    if cells[0] == 'Profondità':
                        try:
                            self.depth = float(self.response_content[4][1].strip("m").replace(".", ""))
                        except:
                            pass
                        else:
                            if self.interpolator is not None:
                                self.t_down = float(self.interpolator(self.depth))

                placement_table = doc.xpath("//span[contains(text(), 'Ubicazione')]/following-sibling::table[1]")[0]
                for row in placement_table.xpath(".//tr"):
                    cells = row.xpath(".//td/text()")
                    if cells[0] == 'Zona marina':
                        self.placement = "Sea"
                    elif cells[0] == 'Provincia':
                        self.placement = "Land"

    def __init_thermo_well(self, bhe_in):

        # Geometric Params
        bhe_geom = BHEGeometry()
        bhe_geom.depth = self.depth

        # Reservoir Params
        res_prop = ReservoirProperties()
        res_prop.grad = self.t_down / self.depth

        self.thermo_well = BaseBHE(bhe_in, reservoir_properties=res_prop, geometry=bhe_geom)
        self.thermo_well.set_HX_condition()

    def optimize_h_rel(self):

        result = minimize_scalar(

            (lambda x: self.evaluate_LCOx(x)[0]),
            bounds=(0.001, 0.99), method='bounded'

        )

        if result.success:

            self.h_rel_opt = result.x
            self.lcoh_min, self.l_horiz_opt = self.evaluate_LCOx(self.h_rel_opt)

    def evaluate_LCOx(self, h_rel):

        # <-- EVALUATE L_HORIZ ------------------------------->
        d_well = self.thermo_well.geom.d_well
        integrator = self.thermo_well.integrator

        UA = integrator.evaluate_integral(h_rel)
        UdAs = np.pi * d_well / self.thermo_well.res_prop.evaluate_rel_resistance(times=[self.time], d=d_well)[0]
        l_horiz = UA / UdAs * integrator.dh_max * self.m_dot
        q_out = integrator.dh_max * self.m_dot * h_rel / 1e3

        if l_horiz > 0:

            # <-- EVALUATE LCOH ---------------------------------->
            lcoh, c_well = self.economic_evaluator.LCOx(

                useful_effects=q_out,
                l_overall=l_horiz,
                d_well=d_well,
                other_costs=1e6,
                w_net_el=0.

            )

        else:

            lcoh = np.inf

        return lcoh, l_horiz

    def __lt__(self, other):

        if self.can_be_optimized and not other.can_be_optimized:
            return True

        elif not self.can_be_optimized and other.can_be_optimized:
            return False

        elif not self.can_be_optimized and not other.can_be_optimized:
            return self.name < other.name

        else:
            return self.lcoh_min < other.lcoh_min

    def __eq__(self, other):

        if self.can_be_optimized and not other.can_be_optimized:
            return False

        elif not self.can_be_optimized and other.can_be_optimized:
            return False

        elif not self.can_be_optimized and not other.can_be_optimized:
            return self.name == other.name

        else:
            return self.lcoh_min == other.lcoh_min


class WellOptimizerList:

    fluid = "Water"
    t_in = 10   # [°C]
    bhe_in = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    bhe_in.set_variable("T", t_in + 273.15)
    bhe_in.set_variable("P", 1e5)

    def __init__(self, kml_file, gml_fonder, stop_after=-1):

        self.kml_file = kml_file
        self.temp_interp = GMLInterpolator(gml_fonder)

        self.__init_wells(stop_after)

    def __init_wells(self, stop_after):

        self.wells = list()
        gdf = gpd.read_file(self.kml_file, driver='KML')

        with open(self.kml_file, 'r') as file:
            doc = parser.parse(file).getroot()

        placemarks = doc.findall('.//{http://www.opengis.net/kml/2.2}Placemark')

        pbar = tqdm(total=len(placemarks[:stop_after]))
        for placemark in placemarks[:stop_after]:

            self.wells.append(WellOptimizer(placemark, gdf, self.temp_interp, self.bhe_in))
            if self.wells[-1].can_be_optimized:
                self.wells[-1].optimize_h_rel()

            pbar.update(1)

        pbar.close()

    @property
    def land_based_wells(self):

        return_list = list()
        for well in self.wells:
            if well.placement == "Land":
                return_list.append(well)

        return return_list

    @property
    def optimizable_wells(self):

        return_list = list()
        for well in self.wells:
            if well.can_be_optimized:
                return_list.append(well)

        return return_list

    def get_best_wells(self, n_best=5):

        well_sorted = sorted(optimizer.wells)
        return well_sorted[:n_best]

    def get_summary(self):

        optimizable_wells = self.optimizable_wells
        summary = np.empty((len(optimizable_wells), 5))
        summary[:, :] = np.nan

        for i, well in enumerate(optimizable_wells):

            summary[i, 0] = well.depth
            summary[i, 1] = well.thermo_well.res_prop.grad
            summary[i, 2] = well.lcoh_min
            summary[i, 3] = well.l_horiz_opt
            summary[i, 4] = well.h_rel_opt

        return summary


# %%------------   INITIALIZE OPTIMIZER                   -----------------------------------------------------------> #
optimizer = WellOptimizerList(kml_file, GML_TEMP_FOLDER, stop_after=-1)
best_wells = optimizer.get_best_wells(n_best=5)


# %%------------   TEST WELL CLASS                        -----------------------------------------------------------> #
h_rels = np.linspace(0, 1, 350)[1:-1]
lcoh_res = np.zeros(h_rels.shape)
l_horiz_res = np.zeros(h_rels.shape)

fig, axs = plt.subplots(1, 2, figsize=(12, 4))

# FIRST AX (OPTIMIZATION PROCESS)
ax = axs[0]
ax_horiz = ax.twinx()

for m_dot in [5, 10, 15]:

    optimizer.wells[0].m_dot = m_dot
    for i, h_rel in enumerate(h_rels):

        lcoh_res[i], l_horiz_res[i] = optimizer.wells[0].evaluate_LCOx(h_rel)

    line, = ax.plot(h_rels, lcoh_res*100)
    ax_horiz.plot(h_rels, l_horiz_res/1e3, "--", color=line.get_color(), label='m_dot = {}kg/s'.format(m_dot))

ax.set_xlabel("$h_{rel}$")
ax.set_ylabel("LCOH [c€/kWh]")
ax_horiz.set_ylabel("$l_{horiz}$ [km]")
ax.set_yscale("log")
ax_horiz.set_yscale("log")
ax_horiz.legend()

# SECOND AX (OPTIMAL POINT DESCRIPTION)
ax = axs[1]
ax_horiz = ax.twinx()
m_dot_list = np.linspace(5, 15, 30)
opt_lcoh = np.zeros(m_dot_list.shape)
opt_l_horiz = np.zeros(m_dot_list.shape)

for j, m_dot in enumerate(m_dot_list):

    optimizer.wells[0].m_dot = m_dot
    result = minimize_scalar((lambda x: optimizer.wells[0].evaluate_LCOx(x)[0]), bounds=(0.001, 0.99), method='bounded')
    opt_lcoh[j], opt_l_horiz[j] = optimizer.wells[0].evaluate_LCOx(result.x)

ax.plot(m_dot_list, opt_lcoh*100)
ax_horiz.plot(m_dot_list, opt_l_horiz / 1e3, "--")
ax.set_xlabel("m_dot [kg/s]")
ax.set_ylabel("LCOH [c€/kWh]")
ax_horiz.set_ylabel("$l_{horiz}$ [km]")
plt.suptitle(f"Optimization of the \"{optimizer.wells[0].name}\" well")
plt.tight_layout()
plt.savefig(os.path.join(CURRENT_DIR, "0 - Output", f"optimization_process.png"), dpi=300)
plt.show()


# %%------------   EVALUATE MAP OF ITALY                  -----------------------------------------------------------> #
n = 5
optimizable_wells = sorted(optimizer.optimizable_wells)[1:]

def idw_interpolation(longitude, latitude, power=2):

    point = Point(longitude, latitude)

    if optimizer.temp_interp.get_interpolator(longitude, latitude) is not None:

        values = 0.
        weights = 0.

        for well in optimizable_wells:

            dist = well.position.distance(point)
            if dist == 0:
                values = well.lcoh_min
                weights = 1
                break

            else:

                w = (1 / dist) ** power
                values += well.lcoh_min * w
                weights += w


        return values/weights

    else:

        return np.nan

lat_min, lat_max = 35.5, 47.1
lon_min, lon_max = 6.6, 18.5
n_points = 400j
grid_lon, grid_lat = np.mgrid[lon_min:lon_max:n_points, lat_min:lat_max:n_points]

pbar = tqdm(total=grid_lon.shape[0] * grid_lon.shape[1], desc="Calculating map")
grid_z = np.zeros_like(grid_lon)
for i in range(grid_lon.shape[0]):
    for j in range(grid_lon.shape[1]):
        grid_z[i, j] = idw_interpolation(grid_lon[i, j], grid_lat[i, j], power=n)
        pbar.update(1)

pbar.close()


# %%------------   PLOT MAP OF ITALY                      -----------------------------------------------------------> #
best_wells = optimizer.get_best_wells(n_best=61)[1:]
plt.figure(figsize=(10, 10))

for well in best_wells[30:]:
    plt.plot(well.position.x, well.position.y, marker='*', markersize=15, color='#e38614', alpha=1)

for well in best_wells[10:30]:
    plt.plot(well.position.x, well.position.y, marker='*', markersize=15, color='#adadad', alpha=1)

for well in best_wells[:10]:
    plt.plot(well.position.x, well.position.y, marker='*', markersize=15, color='#FFD700')

plt.pcolormesh(grid_lon, grid_lat, np.log10(grid_z*100), shading='auto', cmap='viridis')
plt.colorbar(label='$log_{10}(LCOH)$ [c€/kWh]')

plt.title(f'LCOH map')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.savefig(os.path.join(CURRENT_DIR, "0 - Output", f"LCOH map.png"), dpi=300)
plt.show()


# %%------------   GENERATE BEST WELL LIST                -----------------------------------------------------------> #
output_txt = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(

    "Name", "Depth[m]", "T_down [C]",
    "l_horiz [m]", "LCOH [cEURO/kWh]", "W_dot [MW]",
    "link"

)
best_wells = optimizer.get_best_wells(n_best=500)[1:]

for well in best_wells:

    output_txt += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(

        well.name, well.depth, well.t_down,
        well.l_horiz_opt, well.lcoh_min, well.w_dot_opt / 1e6,
        well.link

    )

with open(os.path.join(CURRENT_DIR, "0 - Output", f"best_wells.csv"), "w") as f:
    f.write(output_txt)