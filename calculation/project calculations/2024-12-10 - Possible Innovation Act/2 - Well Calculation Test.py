# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes.geothermal_system.base_bhe import (

    isobaricIntegral, ThermodynamicPoint,
    ReservoirProperties, baseEconomicEvaluator

)
from scipy.interpolate import interp1d
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


class WellOptimizer:

    fluids = ["Water", "CO2"]

    def __init__(self, well_placemark, temp_gdf, temp_iterp):

        self.__placemark = well_placemark
        self.__temp_iterp = temp_iterp
        self.__gdf = temp_gdf

        self.__init_attributes()
        self.__set_main_information()
        self.__analyze_extended_attr()
        self.__retrieve_information_online()

        self.__init_thermo_well()

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

    def __init_thermo_well(self):

        pass

    def optimize_h_rel(self):

        pass


class WellOptimizerList:

    def __init__(self, kml_file, gml_fonder):

        self.kml_file = kml_file
        self.temp_interp = GMLInterpolator(gml_fonder)

        self.__init_wells()

    def __init_wells(self):

        self.wells = list()
        gdf = gpd.read_file(self.kml_file, driver='KML')

        with open(self.kml_file, 'r') as file:
            doc = parser.parse(file).getroot()

        placemarks = doc.findall('.//{http://www.opengis.net/kml/2.2}Placemark')

        pbar = tqdm(total=len(placemarks[:10]))
        for placemark in placemarks[:10]:

            self.wells.append(WellOptimizer(placemark, gdf, self.temp_interp))
            pbar.update(1)

        pbar.close()

    @property
    def land_based_wells(self):

        return_list = list()
        for well in self.wells:
            if well.placement == "Land":
                return_list.append(well)

        return return_list


# %%------------   TEST WELL CLASS                        -----------------------------------------------------------> #
optimizer = WellOptimizerList(kml_file, GML_TEMP_FOLDER)
