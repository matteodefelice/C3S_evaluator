from .C3S_data_loader import *
import xarray as xr
import xskillscore as xs       
import numpy as np
import os.path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

class evaluator:
    """
    This class is used to evaluate and store the skill data of a seasonal forecast created by 
    `C3S_data_loader`. This class can compute deterministic and probabilistic scores (provided by the `xskillscore` package).

    The skill information is saved in the same folder of the processed seasonal data with the suffix `.solved`
    """
    def __init__(self, dl:C3S_data_loader, quiet = False):
        """
        Initialise using an object `C3S_data_loader`. If computed skill data is found then it is loaded.

        Arguments:
            dl (C3S_data_loader): the object with the retrieved forecast/reanalysis pair
            quite (boolean): define if the execution should print out information  
        """
        self._basename = dl.get_basename()
        # Get the forecast/observation data from `dl`
        self._data = dl.get_data()
        # Initialise the xarray Dataset with the skill information copying the `data` object
        # removing the two fields containing forecasts and observations
        self._skill_data = self._data.drop(['fct_var', 'obs_var']) 
        self.quiet = quiet
        # Check if exists "solved" data file
        if os.path.exists(dl.get_basename() + '.solved.nc'):
            if not self.quiet:
                print(f"Loading existing computed skill data {dl.get_basename() + '.solved.nc'}")
            self._skill_data = xr.open_dataset(dl.get_basename() + '.solved.nc')
           

    def compute_deterministic(self):
        """
        Calculate deterministic metrics comparing the ensemble mean with the reanalysis data.
        
        The following metrics are computed:
        - Spearman correlation coefficient (both with the p-value and a masked correlation using p < 0.001)
        - Climatologies 
        - Bias (observation minus forecast)
        - RMSE
        - MAE
        - Normalised MAE (NMAE): ratio between MAE and climatology


        """
        # DETERMINISTIC
        if not self.quiet:
            print('Computing ensemble mean')
        fct_ens_mean = self._data['fct_var'].mean('number')
        if not self.quiet:
            print('Computing correlation')
        self._skill_data = self._skill_data.assign(spearman = xs.spearman_r(fct_ens_mean, self._data['obs_var'], dim = 'year'))
        self._skill_data = self._skill_data.assign(spearman_pvalue = xs.spearman_r_p_value(fct_ens_mean, self._data['obs_var'], dim = 'year'))
        self._skill_data = self._skill_data.assign(spearman_p001 = self._skill_data.spearman.where(self._skill_data.spearman_pvalue <= 1e-2))
        if not self.quiet:
            print('Computing climatology and bias')
        self._skill_data = self._skill_data.assign(clim_obs = self._data.obs_var.mean(dim = 'year'))
        self._skill_data = self._skill_data.assign(clim_fct = fct_ens_mean.mean(dim = 'year'))
        self._skill_data = self._skill_data.assign(clim_bias = self._skill_data.clim_obs - self._skill_data.clim_fct)
        if not self.quiet:
            print('Computing RMSE and MAE/NMAE')
        self._skill_data = self._skill_data.assign(rmse = xs.rmse(fct_ens_mean, self._data['obs_var'], dim = 'year'))
        self._skill_data = self._skill_data.assign(mae  = xs.mae(fct_ens_mean, self._data['obs_var'], dim = 'year'))
        self._skill_data = self._skill_data.assign(nmae  = xs.mae(fct_ens_mean, self._data['obs_var'], dim = 'year') / self._skill_data.clim_obs)
        if not self.quiet:
            print(f"Save skill data to {self._basename + '.solved.nc'}")
        self._skill_data.to_netcdf(self._basename + '.solved.nc')
        
    def compute_probabilistic(self):
        """
        Calculate probabilistic metrics comparing the seasonal forecast with the reanalysis data.
        
        The following metrics are computed:
        - CRPS
        - upper/lower Brier Skill Scores (using terciles)
        - RPSS (using terciles)
        """

        # check if nan exists
        if (np.isnan(self._data['obs_var']).any() or np.isnan(self._data['fct_var']).any()):
            nan_are_present = True
        else:
            nan_are_present = False
    	# PROBABILISTIC scores -----------------------------------------------------------
        if not self.quiet:
            print('Computing CRPS')
        self._skill_data = self._skill_data.assign(crps = xs.crps_ensemble(self._data['obs_var'], self._data['fct_var'], member_dim='number', dim = 'year'))
        
        if not self.quiet:
            print('Computing upper BSS')
        obs_thres = self._data['obs_var'] >  self._data['obs_var'].quantile(2/3, dim = 'year', skipna = nan_are_present)
        ubs = xs.brier_score(obs_thres, (self._data['fct_var'] > self._data['fct_var'].quantile(2/3, dim = 'year', skipna = nan_are_present)).mean("number"), dim = 'year')
        ubaseline = xs.brier_score(obs_thres, xr.DataArray(1/3), dim = 'year')
        self._skill_data = self._skill_data.assign(ubss = 1 - (ubs/ubaseline))
        if not self.quiet:
            print('Computing lower BSS')
        obs_thres = self._data['obs_var'] <  self._data['obs_var'].quantile(1/3, dim = 'year', skipna = nan_are_present)
        lbs = xs.brier_score(obs_thres, (self._data['fct_var'] < self._data['fct_var'].quantile(1/3, dim = 'year', skipna = nan_are_present)).mean("number"), dim = 'year')
        lbaseline = xs.brier_score(obs_thres, xr.DataArray(1/3), dim = 'year')
        self._skill_data = self._skill_data.assign(lbss = 1 - (lbs/lbaseline))
        # Save also scores 
        self._skill_data = self._skill_data.assign(ubs = ubs)
        self._skill_data = self._skill_data.assign(lbs = lbs)

        if not self.quiet:
            print('Computing RPSS')
        obs_zeromean = self._data['obs_var'] - self._data['obs_var'].mean(dim = 'year')
        fct_zeromean = self._data['fct_var'] - self._data['fct_var'].mean(dim = 'year')
        cat_edges = obs_zeromean.quantile(q = [1/3, 2/3], dim = 'year', skipna = nan_are_present).rename({'quantile':'category_edge'})
        self._skill_data = self._skill_data.assign(rps = xs.rps(obs_zeromean, fct_zeromean, cat_edges, member_dim='number', dim = 'year'))
        rps_cl = xs.rps(obs_zeromean, obs_zeromean.mean(dim = 'year').expand_dims({'number':1}), cat_edges, dim = 'year', member_dim='number')
        self._skill_data = self._skill_data.assign(rpss = 1 - (self._skill_data.rps/rps_cl))
        
        if not self.quiet:
            print(f"Save skill data to {self._basename + '.solved.nc'}")
        self._skill_data.to_netcdf(self._basename + '.solved.nc')

    def compute_classification_scores(self):
        """
        Calculate classification metrics comparing the seasonal forecast with the reanalysis data. 
        The prediction from the seasonal forecasts' ensemble is considered True when more than 50% of the members agree. 
        
        The following metrics are computed:
        - Hit Rate (Upper and Lower)
        """

        # check if nan exists
        if (np.isnan(self._data['obs_var']).any() or np.isnan(self._data['fct_var']).any()):
            nan_are_present = True
        else:
            nan_are_present = False
    	# CLASSIFICATION scores -----------------------------------------------------------
        if not self.quiet:
            print('Computing Upper Hit Rate')
        
        obs_thres = self._data['obs_var'] > self._data['obs_var'].quantile(2/3, dim = 'year', skipna = nan_are_present)
        fct_thres = (self._data['fct_var'] > self._data['fct_var'].quantile(2/3, dim = 'year', skipna = False)).mean("number")

        o_category_edges = np.array([0, 1, 2]) 
        f_category_edges = np.array([0, 0.5, 1]) # It's True when >50% of the members agree

        # Calculate contingency
        cc = xs.Contingency(obs_thres, fct_thres, 
               o_category_edges, f_category_edges,
               dim='year') 

        self._skill_data = self._skill_data.assign(uhr = cc.hit_rate())
        if not self.quiet:
            print('Computing Lower Hit Rate')

        obs_thres = self._data['obs_var'] < self._data['obs_var'].quantile(1/3, dim = 'year', skipna = nan_are_present)
        fct_thres = (self._data['fct_var'] < self._data['fct_var'].quantile(1/3, dim = 'year', skipna = False)).mean("number")

        o_category_edges = np.array([0, 1, 2]) 
        f_category_edges = np.array([0, 0.5, 1]) # It's True when >50% of the members agree

        # Calculate contingency
        cc = xs.Contingency(obs_thres, fct_thres, 
               o_category_edges, f_category_edges,
               dim='year') 

        self._skill_data = self._skill_data.assign(lhr = cc.hit_rate())

        if not self.quiet:
            print(f"Save skill data to {self._basename + '.solved.nc'}")
        self._skill_data.to_netcdf(self._basename + '.solved.nc')

    def get_data(self) -> xr.Dataset:
        """
        Return the xarray Dataset containing seasonal data and ERA5
        """
        return(self._data)
    def set_data(self, data:xr.Dataset):
        """
        Set a new xarray Dataset containing seasonal data and ERA5. 
        """
        self._data = data
        if not self.quiet:
            print('Forecast/observation data has changed, please recompute metrics and scores for consistency')
    def get_skill_data(self) -> xr.Dataset:
        """
        Return the xarray Dataset containing only the skill and metrics 
        """
        return(self._skill_data)
    
    def save_maps(self):
        """
        Save to disk a set of plots for the computed metrics. 
        """
        
        variables   = ['spearman', 'rmse', 'mae', 'crps', 'ubss', 'lbss', 'clim_bias', 'nmae', 'spearman_p001', 'rps', 'rpss', 'uhr', 'lhr']
        labels      = ['spearman correlation', 'RMSE', 'MAE', 'CRPS', 'Upper BSS', 'Lower BSS', 'Clim. bias (obs-fct)', 'NMAE', 'spearman corr. (p-value < 1e-2)', 'RPS', 'RPSS', 'Hit Rate (Upper)', 'Hit Rate (Lower)']
        for this_fig in zip(variables, labels):
            if this_fig[0] in self._skill_data.data_vars.keys():
                plt.figure(figsize=(1200/300, 800/300), dpi=300)
                ax = plt.axes(projection=ccrs.Robinson())
                p = self._skill_data[this_fig[0]].plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar = False)
                ax.set_global(); ax.coastlines(linewidth=0.3);
                plt.title(os.path.basename(self._basename),fontsize=8)
                plt.suptitle(this_fig[1],fontsize=14, y=0.9)
                plt.colorbar(p,fraction=0.025, pad=0.04)
                plt.savefig(f'{self._basename}-{this_fig[0]}.png', dpi = 300)
                plt.close()
            else:
                if not self.quiet:
                    print(f"{this_fig[0]} not present")
                
