#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:59:15 2021

@author: christiellyborges
"""

"""
    Supplementary information for
    
    # article title #
    
    This model was created as a process-based simulation modelling
    to examine the distribution of human groups and their languages
    in Neotropical realm (present-day Latin America). 
  
"""

########################
### Import libraries ###
########################

import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import random
import os
from sklearn.linear_model import LinearRegression



log = logging.getLogger(__name__)
pd.options.mode.chained_assignment = None  # default='warn'


######################################
### Carrying capacity (k) function ###
######################################


# function for K
class k:
    
    @property
    def data(self):
        if self.df is None:
            return None
        return (self.df)
    
    
    def __init__(self):
        self.df = gpd.read_file(os.path.join(path,"RAINFALL.shp")) # import rainfall data
        self.df["K"] = None # empty column for K
        
    def showData(self):
        print(self.df)
        
    def power(
            self,
            alpha: float=None,
            beta: float=None
        ):
        """"
        Calculates the carrying capacity as a power function.
        Formula: K = alpha * (Precipitation ^ beta)
        --------------
        alpha: unkown parameter, serves as simple scaling factor
        beta: unkown parameter, governs the rate of growth of k with the increase in precipitation
        env: mean annual precipitation
        Because alpha is a very small number, we estimate log10(alpha). 
        So, to calculate K with the estimated numbers, one has to compute K = (10^alpha) * (Precipitation ^ beta)
        Returns
        -------
        : pandas.Dataframe
            table with a k for each cell
        """
        df=self.df
        env=round(self.df._mean)
        k = pow(10, alpha) * pow(env, beta)
        if ((k==0) | (k<1)).any: # cells must be occupied by at least one individual
            pass
    
        df = df.assign(K=k)
        return df
    
    def exponential(
            self,
            alpha: float=None,
            beta: float=None
        ):
        """
        Calculates the carrying capacity as an exponential function.
        Formula: K = alpha + e ^ (beta * Precipitation)
        --------------
        alpha: unkown parameter, serves as simple scaling factor
        beta: unkown parameter, governs the rate of growth of k with the increase in precipitation. Is a 
        multiplier of precipitation, which is in the exponent of euler number (natural exponent)
        env: mean annual precipitation
        Returns
        -------
        : pandas.Dataframe
            table with a k for each cell
        """
        df=self.df
        env=self.df._mean
        k = alpha + np.exp(pow(10,beta) * env) # formula do MT
        if ((k==0) | (k<1)).any: # cells must be occupied by at least one individual
            pass
        df = df.assign(K=k)
        return df
    
    def logistic(
            self,
            x0: float=None,
            s: float=None,
            MaxY: float=None
        ):
        """
        Calculates the carrying capacity as a logistic function.
        Formula: K = x0 / (1 + EXP((-1*(10^s)))*(Precipitation-MaxY)))
        --------------
        env: mean annual precipitation
        x0: the value of precipitation in which the curve is in its midpoint (point of inflection of the curve)
        s: is the steepness of the curve ; because s is a very small number, we estimate log10(s)
        m: is the maximum carrying capacity (saturation level), in units of carrying capacity
        Returns
        -------
        : pandas.Dataframe
            table with a k for each cell
        """
        df=self.df
        env=self.df._mean
        k = x0 / (1 + np.exp((-1 * (pow(10,s))) * (env - MaxY)))
        if ((k==0) | (k<1)).any: # cells must be occupied by at least one individual
            pass
        df = df.assign(K=k)
        return df


################################
### Define a k for each cell ###
################################

# identify the neighbors and store the carrying capacity
# define a k for each cell
def id_neighbors_and_k(neigh, k):

    # cells with no neighbors must equal zero    
    for i in range(0, len(neigh)):
        if (neigh.NEIGHBORS[i] is None):
            neigh.NEIGHBORS[i] = "0"
    
    k = getattr(k, "tolist", lambda: k)()

    # splitting neighbors into separate rows 
    neighsplit = pd.DataFrame(neigh.NEIGHBORS.str.split(',').tolist(), index=neigh.Id).stack()
    neighsplit = neighsplit.reset_index([0, 'Id'])
    neighsplit.columns = ['Id', 'NEIGHBORS']

    neighsplit["K"] = ""
    for i in range(len(neighsplit)):
        cellid = int(neighsplit.NEIGHBORS[i])
        if cellid == 0:
            neighsplit.K[i] = 0
        else:
            kvalue = k[k['Id']==cellid]
            neighsplit.K[i] = kvalue.K.iloc[0].round(4)
    
    return(neighsplit)
   

###############################
### Maximum population size ###
###############################

def get_maxpop(population):
    pop = population.tlpop.sample(n=233)
    random.shuffle(pop)
    return(pop)


#########################
### Mechanistic model ###
#########################

# populate Central + South America!    
def populate(neighk, k, grid, population):
    # cleaning the neighbors dataframe
    neighk = getattr(neighk, "tolist", lambda: neighk)()
    neighk["NEIGHBORS"] = neighk["NEIGHBORS"].astype(int)
    neighk["Id"] = neighk["Id"].astype(int)
    neighk["Occupy"] = 0 # column to keep track of occupancy

    # starting grid ks
    k["Occupy"] = 0 # column to keep track of occupancy
    k["Language"] = None
    k["K"] = k.K.round(0) #rounding K values
    
    # virtual languages
    langs = 1000000
    
    # start anywhere
    start = int(grid.Id.sample(n=1))    
    
    counter = 0
 
    while counter < len(k):
        
        for l in range(langs):
            
            #print("Virtual language:", l)
            #print("Cells occupied:", counter)
            
            # get language population
            pop = get_maxpop(population)[0] # random maximum population size
            #print("Population is:", pop)
                         
            while pop > 0:
                # occupy cells with one language population
                if int(k.loc[k.Id.eq(start), "Occupy"]) == 0:
                    #occupy cell and subtract population
                    k.loc[k.Id.eq(start), "Occupy"] = 1
                    k.loc[k.Id.eq(start), "Language"] = "lang"+str(l)
                    
                    # also fill in cell in neighbors df
                    mask = neighk[neighk.NEIGHBORS.eq(start)].index 
                    neighk.loc[mask,'Occupy'] = 1
                    
                    remainpop = int(k[k.Id.eq(start)].K)
                    pop -= remainpop
                    #print("1 Remaining population is:", pop)
                    if pop < 0:
                        #print("Population is filled!")
                        break
                    
                    # occupy neighbors                    
                    if k["Language"].str.contains("lang"+str(l)).any() == True:
                        #print("Occupy neighbor cells.")
                        colonized = k.loc[k.Language.eq("lang"+str(l))]
                        ids = colonized.Id.tolist()
                        mask = neighk["Id"].isin(ids)
                        neighbors = neighk[mask]
                        
                        if (neighbors["Occupy"]==0).any():
                            neighbors = neighbors.loc[neighbors.Occupy.eq(0)]
                            whichcell = neighbors.NEIGHBORS
                            moves = neighbors.shape[0]
                                
                            for m in range(moves):
                                if (neighbors.loc[neighbors.index[m], "Occupy"] == 0):
                                    neighk.loc[neighk.NEIGHBORS.eq(whichcell.iloc[m]), 'Occupy'] = 1
                                    k.loc[k.Id.eq(whichcell.iloc[m]), "Occupy"] = 1
                                    k.loc[k.Id.eq(whichcell.iloc[m]), "Language"] = "lang"+str(l)
                                    neighbors.loc[neighbors.index[m], "Occupy"] = 1
                                    remainpop = neighbors.loc[neighbors.index[m], "K"]
                                    pop -= remainpop
                                    #print("2 Remaining population is:", pop)
                                    if pop < 0:
                                        #print("Population reached its maximum!")
                                        break
                                    
                                    if (neighbors.Occupy.sum() == len(neighbors)):
                                        #print("needs another go")
                                        start = int(neighbors.NEIGHBORS.sample(n=1))
                                        break
                                
                # if starting cell is already occupied move to other neighbors                                      
                else:
                    #print("Starting cell is occupied. Getting a new one.")
                    
                    if k["Language"].str.contains("lang"+str(l)).any() == True:
                        colonized = k.loc[k.Language.eq("lang"+str(l))]
                        ids = colonized.Id.tolist()
                        mask = neighk["Id"].isin(ids)
                        neighbors = neighk[mask]
                        
                        if (neighbors["Occupy"]==0).any():
                            neighbors = neighbors.loc[neighbors.Occupy.eq(0)]
                            whichcell = neighbors.NEIGHBORS
                            moves = neighbors.shape[0]
                                
                            for j in range(moves):
                                if (neighbors.loc[neighbors.index[j], "Occupy"] == 0):
                                    neighk.loc[neighk.NEIGHBORS.eq(whichcell.iloc[j]), 'Occupy'] = 1
                                    k.loc[k.Id.eq(whichcell.iloc[j]), "Occupy"] = 1
                                    k.loc[k.Id.eq(whichcell.iloc[j]), "Language"] = "lang"+str(l)
                                    neighbors.loc[neighbors.index[j], "Occupy"] = 1
                                    remainpop = neighbors.loc[neighbors.index[j], "K"]
                                    pop -= remainpop
                                    #print("3 Remaining population is:", pop)
                                    if pop < 0:
                                        #print("Population reached its maximum!")
                                        break
                                    
                                    if (neighbors.Occupy.sum() == len(neighbors)):
                                        #print("needs another go")
                                        start = int(neighbors.NEIGHBORS.sample(n=1))
                                        break
                            
                        # find an empty neighbor cell!
                        else:
                            #print("No neighbors are available.")
                            remainpop = pop+1
                            pop -= remainpop
                            break
                                                                   
                           
                    
            if pop <= 0:
                #print("Population is filled. Moving to next language.")
                # should be from ANY neighbor cell that is already occupied
                counter = k.Occupy.sum()
                occupied = True
                colonize = 0
                while occupied == True:
                    poscell = neighk.loc[neighk.Occupy.eq(1)]
                    ids = poscell.NEIGHBORS.tolist()
                    mask = neighk['Id'].isin(ids)
                    neighbors = neighk[mask]
                    colonize += 1
                    #start = int(poscell.NEIGHBORS.sample(n=1))
                    #neighbors = neighk[neighk.Id.eq(start)]
                    if (neighbors["Occupy"]==0).any():
                        neighbors = neighbors.loc[neighbors.Occupy.eq(0)]
                        start = int(neighbors.NEIGHBORS.sample(n=1))
                        
                        # check availability of neighbors
                        neighbors = neighk[neighk.Id.eq(start)]
                        if (neighbors.Occupy.sum() < len(neighbors)):
                            occupied = False
                        else:
                            occupied = True
                            
                    if colonize == 1000:
                        print("Populations cannot move.")
                        counter = counter + 10000
                        break
                        
            
                if counter >= len(k):
                    print("End simulation.")
                    break                           
        
    return k, l


########################
### Model validation ###
########################

# create richness maps
def create_richness(grid, model):
    # SPATIAL JOIN
    grid_pred = grid.copy(deep=True)
    grid_pred['Richness'] = 0
    #empty = {key: np.inf for key in model.Language.unique()}
    grid_pred['Richness'] = [[] for i in range(len(grid_pred))]
    spatial_index = grid_pred.sindex
    for idx, lang in model.iterrows():
        possible_matches_index = list(spatial_index.intersection(lang.geometry.bounds))
        possible_matches = grid_pred.iloc[possible_matches_index]    
        precise_matches = possible_matches[possible_matches.intersects(lang.geometry)]
        grid_pred.iloc[precise_matches.index].apply(lambda x: x['Richness'].append(lang.Language), axis=1)
    grid_pred['Richness'] = grid_pred.apply(lambda x: len(set(x['Richness'])), axis=1)
    
    return(grid_pred)
   

# ad hoc goodness of fit index (f)
# from f, later choose the best 200 models (parameter combination)
def goodness_of_fit(obs_map, pred_map, lang):
    obs=986 #obs languages
    s=1 - (abs(obs - lang)/lang)
    
    # coeff determination
    x=np.array(obs_map.Richness).reshape((-1, 1))
    y=np.array(pred_map.Richness)
    reg=LinearRegression().fit(x, y)
    r_sq=reg.score(x, y)
    
    #fit index (f)
    f=(r_sq+s)/2
    return s, r_sq, f


######
# implement f

def simulate_languages(parameters, hexic_grid, grid_300, obs_grid, population, start, save="no"):   
    
    replicates = parameters.shape[0]
    fit_index = pd.DataFrame(columns = ["Model","Alpha","Beta","Languages","s","r2","fit_index"])
    pred_mods = []
    pred_rich = []
    
    for i in range(start, replicates):
        print("Running model:", i)
        alpha, beta = parameters.iloc[i]
        k = carrycap.power(alpha=alpha, beta=beta)
        neighk = id_neighbors_and_k(hexic_grid, k)
        populated_df, lang = populate(neighk, k, neigh, population)
        model_rich = create_richness(grid300, populated_df)
        if lang <= 0:
            print ("Languages simulated were less than 0")
            lang = 1
            next
        s, r2, f = goodness_of_fit(obs_map, model_rich, lang)
        fit_index.loc[i] = ["Model"+str(i), alpha, beta, (lang), s, r2, f]
        #pred_mods.append(populated_df)
        #pred_rich.append(model_rich)
        pred_mods=populated_df
        pred_rich=model_rich
        
        # save total languages per model df
        fit_index.to_csv(os.path.join(path_results,"General-Results/","Goodness_of_fit.csv"))

    
        # SAVE BEST OUTPUTS
        if (save=="yes"):
            if( f > 0.3):
                # predicted models
                pred_mods.to_csv((os.path.join(path_results, "Models-hexic-map/", "Model"+str(i)+".csv")))
                # richness maps
                pred_rich.to_file(os.path.join(path_results, "Models-Richness/", "Model"+str(i)+".shp"))
            
    return fit_index, pred_rich


    


########################
### Running Functions ##
########################

path = "/path/to/somewhere"
path_grid="/path/to/layers"
path_results="/path/to/results/"

carrycap=k()
neigh = gpd.read_file(os.path.join(path, "file.shp"))
grid300=gpd.read_file(os.path.join(path_grid, "fileshp"))
obs_map=gpd.read_file(os.path.join(path_results, "file.shp"))
population = pd.read_csv(os.path.join(path, "population.txt"), sep=" ")
parameters = pd.read_csv(os.path.join(path, "parameters.csv"), sep=",")


# run model
model_fit, mods = simulate_languages(parameters = parameters,
                                     hexic_grid = neigh, 
                                     grid_300 = grid300, 
                                     obs_grid = obs_map,
                                     population = population,
                                     start = 1,
                                     save = "yes")



