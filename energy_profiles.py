import atlite
from atlite.gis import ExclusionContainer
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.plot import show
import os
from pathlib import Path


### --------- Input data and configuration ------------- ###

class M(object):
    pass

# Store all parameters in model object
def InstantiateModel():

    # Make model object from class
    m = M()

    m.years = [str(y) for y in range(1990,2020)]

    m.tech = "OnshoreWind"

    m.solarpanel = atlite.solarpanels.CSi
    m.windturbine_onshore = Path("K:/CodeModules/LAVA_china/LW_farm.yaml")

    m.data_path = "E:/CETO 2025 Spatial Analysis data/GIS/ERA5/China_mainland_1990_2.nc"

    m.output_path = f"E:/CETO 2025 Spatial Analysis data/NeiMongol/{m.tech}_profile.csv"

    m.derate = 0.95    

    return m

# --------------------- Init --------------------- #
dirname = os.getcwd() 
data_path = os.path.join(dirname, 'data', region_folder_name)

EPSG_custom = config['EPSG_manual']
with open(os.path.join(data_path, region_name+'_EPSG.pkl'), 'rb') as file:
    EPSG = pickle.load(file)
if EPSG_custom:
    EPSG = int(EPSG_custom)
EPSG = "EPSG:32650"

atlite_year = 2010

# Load region geometry
#regionPath =os.path.join(data_path, f'{region_name}_{EPSG}.geojson')
regionPath = "E:/CETO 2025 Spatial Analysis data/GIS/Admin boundaries/China_provinces.gpkg"
region = gpd.read_file(regionPath)
region = region[region['name'] == "East-InnerMongolia"]

# Load potential
#potentialPath = os.path.join(data_path, f'available_land_filtered-min{min_pixels_connected}_{region_name}_EPSG{EPSG}.tif')
potentialPath = "E:/CETO 2025 Spatial Analysis data/NeiMongol/available_land_filtered-min5_gradeE_China_EPSG32650.tif"


# --------------------- Main script --------------------- #
m = InstantiateModel()

# Find bounds of the region
bounds = region.to_crs(EPSG).total_bounds

# Adjust
x1, y1, x2, y2 = bounds.values[0]
x=slice(x1 - 0.2, x2 + 0.2),
y=slice(y1 - 0.2, y2 + 0.2),

# Load weather data cutout
cutout = atlite.Cutout(path=m.data_path).sel(x=x, y=y)

# Loading potential
excluder = ExclusionContainer(crs=EPSG)
excluder.add_raster(potentialPath, codes=1, invert=True)

# Availability of the area
masked, transform = excluder.compute_shape_availability(region)
fig, ax = plt.subplots()
excluder.plot_shape_availability(region)
eligible_share = masked.sum() * excluder.res**2 / region.geometry.item().area
print(f"The eligibility share is: {eligible_share:.2%}")

# `A` is an DataArray with 3 dimensions (`shape`, `x`, `y`) and very sparse data. 
# It indicates the relative overlap of weather cell `(x, y)` with geometry `region` while excluding the area specified by the `excluder`. 
A = cutout.availabilitymatrix(region, excluder)

# Aggregation of potential to weather data cells
fig, ax = plt.subplots()
region.plot(ax=ax, edgecolor="k", color="None")
A.plot(ax=ax, cmap="Greens")
#cutout.grid.plot(ax=ax, color="None", edgecolor="grey")

capacity_matrix = A.stack(spatial=["y", "x"])


#---- Heat demand -----

I = cutout.indicatormatrix(region)

pop_layout = xr.open_dataarray(pop_layout)

stacked_pop = pop_layout.stack(spatial=("y", "x"))
M = I.T.dot(np.diag(I.dot(stacked_pop)))

pop_matrix = sp.sparse.csr_matrix(pop_map.T)




# Match case statement for different technologies
match m.tech:
    case "SolarPV":
        df_tech = cutout.pv(
            matrix=capacity_matrix,
            panel=m.solarpanel,
            orientation="latitude_optimal",
            index=region.index,
            per_unit=True
        )
        # Apply derate
        df_tech = df_tech * m.derate

    case "SolarPVTracking":
        df_tech = cutout.pv(
            matrix=capacity_matrix,
            panel=m.solarpanel,
            orientation="latitude_optimal",
            index=region.index,
            tracking="horizontal",
            per_unit=True
        )
        # Apply derate
        df_tech = df_tech * m.derate

    case "OnshoreWind":
        df_tech = cutout.wind(
            matrix=capacity_matrix,
            turbine=m.windturbine_onshore,
            index=region.index,
            per_unit=True
        )
        # Apply derate
        df_tech = df_tech * m.derate

    case "OffshoreWind":
        df_tech = cutout.wind(
            matrix=capacity_matrix,
            turbine=m.windturbine_offshore,
            index=region.index,
            per_unit=True
        )
        # Apply derate
        df_tech = df_tech * m.derate

    case "Hydro":
        plants, basins = hydro()
        df_tech = cutout.hydro(
            plants=plants,
            hydrobasins=basins
        )
        # Apply derate
        df_tech = df_tech * m.derate

    case "HeatDemand":
        df_tech = cutout.heat_demand(
            threshold=15,
            a=1,
            constant=0 # Corresponding to no temperature independent heat demand
        )

    case _:
        raise ValueError(f"Unknown technology: {tech}")



# Set demand to zero outside of heating season
start_day = config["heat_demand_start_day"]
end_day = config["heat_demand_end_day"]
regonal_daily_hd.loc[f"{atlite_year}-{start_day}":f"{atlite_year}-{end_day}"] = 0

# Convert to DataFrame for export 
df_tech = df_tech.to_pandas().to_csv(output_path)




















# --------------------- Misc --------------------- #
df_tech.to_pandas().plot(ylabel="Power [GW]", ls="--", figsize=(10, 4))
wnd.to_pandas().plot(ylabel="Power [GW]", ls="--", figsize=(10, 4))
# Apply derate
df_tech = df_tech * m.derate

df_tech.name = 'Nei-Mongol'

from atlite.windturbines import load_turbine, turbineconfig

def hydro():
    plants = gpd.read_file("E:/CETO 2025 Spatial Analysis data/GIS/Hydro/Hydro_plants.gpkg")
    basins = gpd.read_file("E:/CETO 2025 Spatial Analysis data/GIS/Hydro/Hydrobasin_level4.gpkg")

    return plants, basins

### ---------- Produce profiles ----------- ###



capacity_matrix(m, cutout, excluder)


    EPSG = "EPSG:32650"
    

    # Balmorel time sets and geography
    m.region_file = "E:/CETO 2025 Spatial Analysis data/GIS/Admin boundaries/China_provinces.gpkg"
    m.region_name = 'East-InnerMongolia'
    m.region = gpd.read_file(m.region_file)
    m.region = m.region[m.region['name'] == m.region_name].geometry.to_crs(EPSG)


