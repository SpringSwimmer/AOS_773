from MCEq.core import MCEqRun
import crflux.models as pm
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colorcet as cc

era5_dataset = xr.load_dataset("/home/leoseen/aos_773/final_project/era5_reanalysis_2025.grib", engine = 'cfgrib')
spole_data = era5_dataset.t.sel(latitude = '-90', longitude = '-180.0')

geopotential_dataset = xr.load_dataset("/home/leoseen/aos_773/final_project/geopotential.grib", engine = 'cfgrib')
geopotential_data = geopotential_dataset.sel(latitude = '-90', longitude = '-180.0')
g = 9.8
R = 6378388
geom_height = R*geopotential_data.z/(g*R-geopotential_data.z)

def calculate_density_from_temperature(temperature_data):
    density_data = temperature_data.isobaricInhPa.values*1e2*28.9647/(8.314*temperature_data.values*1e6)
    return density_data

def calculate_density_from_temperature(temperature_data):
    density_data = temperature_data.isobaricInhPa.values*1e2*28.9647/(8.314*temperature_data.values*1e6)
    return density_data

def get_one_solver(height_data, temperature_data):
    density_data = calculate_density_from_temperature(temperature_data)
    length = np.max(height_data.values)-np.min(height_data.values)
    traverse_length = (np.max(height_data.values) - height_data.values)*1e2

    mceq_run = MCEqRun(
            interaction_model='SIBYLL2.3e',
            primary_model = (pm.HillasGaisser2012, "H4a"),
            theta_deg=0.0,
            density_model = ('GeneralizedTarget', None)
            )
    mceq_run.density_model.set_length(length*1e2)
    for tl, rho in zip(np.flip(traverse_length)[:-1], np.flip(density_data)[:-1]):
        if tl == traverse_length[1]:
            print('hi')
        mceq_run.density_model.add_material(tl, rho, 'air')
    mceq_run.solve()
    mu_conv = mceq_run.get_solution('conv_mu+', 3)+ mceq_run.get_solution('conv_mu-', 3)
    return mu_conv

for h_data, temp_data in zip(geom_height, spole_data):
    if not temp_data.time.values == h_data.time.values:
        print('Not the same day')
    if not np.array_equal(temp_data.isobaricInhPa.values,h_data.isobaricInhPa.values):
        print('Not same pressure')
    get_one_solver(h_data, temp_data)