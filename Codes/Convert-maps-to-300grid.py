#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 16:55:47 2021

@author: christiellyborges

### thanks to https://gist.github.com/michellemho/43e817e8e85c2246d75c40277b7a4968
"""

import pandas as pd
import geopandas as gpd
import os
import numpy as np
from shapely import wkt



### import 300km2 grid ###

path_grid="/Data/LatinAmerica/"
grid300=gpd.read_file(os.path.join(path_grid, "Grid_style.shp"))
grid300.crs


######################
### rain data ###
######################

rain_grid="Data/Data-ready/"
rain=gpd.read_file(os.path.join(rain_grid, "LA-HexicRain.shp", ))
obs_shp = rain

grid_obs = grid300.copy(deep=True)
grid_obs['Precipitation'] = 0
empty = {key: np.inf for key in obs_shp._mean.unique()}
grid_obs['Precipitation'] = [[] for i in range(len(grid_obs))]
spatial_index = grid_obs.sindex
for idx, lang in obs_shp.iterrows():
    possible_matches_index = list(spatial_index.intersection(lang.geometry.bounds))
    possible_matches = grid_obs.iloc[possible_matches_index]    
    precise_matches = possible_matches[possible_matches.intersects(lang.geometry)]
    grid_obs.iloc[precise_matches.index].apply(lambda x: x['Precipitation'].append(lang._mean), axis=1)
grid_obs['Precipitation'] = grid_obs.apply(lambda x: len(set(x['Precipitation'])), axis=1)

# save
grid_obs.to_file(os.path.join(rain_grid, "Rain-300g.shp"))

######################
### observed data ###
######################

path_obs = '/Results/Layers/'
obs_shp = gpd.read_file(os.path.join(path_obs,"Lang-obs.shp"))
obs_shp.crs == grid300.crs

grid300.head(2)
obs_shp.head(2)

grid300.shape
obs_shp.shape

# SPATIAL JOIN
#join=gpd.sjoin(grid300, obs_shp, how="inner")

# INDEX ON GRID, LOOP THROUGH CELLS
grid_obs = grid300.copy(deep=True)
grid_obs['Richness'] = 0
empty = {key: np.inf for key in obs_shp.id_1.unique()}
grid_obs['Richness'] = [[] for i in range(len(grid_obs))]
spatial_index = grid_obs.sindex
for idx, lang in obs_shp.iterrows():
    possible_matches_index = list(spatial_index.intersection(lang.geometry.bounds))
    possible_matches = grid_obs.iloc[possible_matches_index]    
    precise_matches = possible_matches[possible_matches.intersects(lang.geometry)]
    grid_obs.iloc[precise_matches.index].apply(lambda x: x['Richness'].append(lang.id_1), axis=1)
grid_obs['Richness'] = grid_obs.apply(lambda x: len(set(x['Richness'])), axis=1)

# save
grid_obs.to_file(os.path.join(path_obs, "OBS-Rich.shp"))


######################
### predicted data ###
######################

path_mod = '/Results/'

#####  save the best models ##
total_lang = pd.read_csv(os.path.join(path_mod,"Total_Lang_2trial.csv"))
#best_langs=total_lang[(total_lang["Languages"] > 963) & (total_lang["Languages"] < 1005)]

# remove _ from names
#best_langs["Model"] = best_langs["Model"].replace({'_':''}, regex=True)

n_files = total_lang.shape[0]
filenames = total_lang["Model"]

filenames = filenames.reset_index()
filenames = filenames["Model"]

mod_list = []
for f in range(len(filenames)):
    print(f)
    mod_list.append(pd.read_csv(os.path.join(path_mod, "Models-hexic-map", filenames[f]+".csv")))


mod_list[0]
mod_list[1]
mod_list[100]
mod_list[199]

# save back only the best
#for save in range(len(filenames)):
#    mod_list[save].to_csv(os.path.join(path_mod, "Best-Models-300-grid", "Model"+str(save)+".csv"))
##### end save best model ####

    
## import best models

import glob

os.chdir(path_mod+"Best-Models-300-grid/")
mod_list = []
for file in glob.glob("*"):
    #mod_list.append(file)
    mod_list.append(pd.read_csv(file))


# SPATIAL JOIN

pred_rich=[]
for j in range(len(mod_list)):
    print("Running model", j)
    pred_shp = mod_list[j]
    pred_shp = pred_shp.drop(pred_shp.columns[[0, 2, 4, 5]], axis=1) # delete unwanted columns
    pred_shp['geometry'] = pred_shp['geometry'].apply(wkt.loads)
    pred_shp = gpd.GeoDataFrame(pred_shp, geometry='geometry')
    
    grid_pred = grid300.copy(deep=True)
    grid_pred['Richness'] = 0
    empty = {key: np.inf for key in pred_shp.Language.unique()}
    grid_pred['Richness'] = [[] for i in range(len(grid_pred))]
    spatial_index = grid_pred.sindex
    for idx, lang in pred_shp.iterrows():
        possible_matches_index = list(spatial_index.intersection(lang.geometry.bounds))
        possible_matches = grid_pred.iloc[possible_matches_index]    
        precise_matches = possible_matches[possible_matches.intersects(lang.geometry)]
        grid_pred.iloc[precise_matches.index].apply(lambda x: x['Richness'].append(lang.Language), axis=1)
    grid_pred['Richness'] = grid_pred.apply(lambda x: len(set(x['Richness'])), axis=1)
    
    pred_rich.append(grid_pred)
    # save files
    grid_pred.to_file(os.path.join(path_mod, "Models-Richness/", "Model"+str(j)+".shp"))


######################
#### language family ####
######################

path_obs = '/Data/Data-ready/'
obs_shp = gpd.read_file(os.path.join(path_obs,"Family-LA.shp"))
obs_shp.crs == grid300.crs

grid300.head(2)
obs_shp.head(2)

grid300.shape
obs_shp.shape

# SPATIAL JOIN
#join=gpd.sjoin(grid300, obs_shp, how="inner")

# INDEX ON GRID, LOOP THROUGH CELLS
grid_obs = grid300.copy(deep=True)
grid_obs['Richness'] = 0
empty = {key: np.inf for key in obs_shp.family.unique()}
grid_obs['Richness'] = [[] for i in range(len(grid_obs))]
spatial_index = grid_obs.sindex
for idx, lang in obs_shp.iterrows():
    possible_matches_index = list(spatial_index.intersection(lang.geometry.bounds))
    possible_matches = grid_obs.iloc[possible_matches_index]    
    precise_matches = possible_matches[possible_matches.intersects(lang.geometry)]
    grid_obs.iloc[precise_matches.index].apply(lambda x: x['Richness'].append(lang.family), axis=1)
grid_obs['Richness'] = grid_obs.apply(lambda x: len(set(x['Richness'])), axis=1)

# save
grid_obs.to_file(os.path.join(path_obs, "OBS-FamRich.shp"))
