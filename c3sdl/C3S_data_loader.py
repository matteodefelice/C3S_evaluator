import cdsapi
import xarray as xr
import numpy as np
import xesmf as xe 
import os.path, sys, yaml, warnings, datetime

class C3S_data_loader:
    """
    This class represents a single (single-level) variable from a seasonal forecast downloaded from the C3S Climate Data Store (CDS)
    (cds.climate.copernicus.eu/). ERA5 reanalysis is associated to the seasonal data for validation purposes.
    Both the fields have the same grid.
    """
    def read_config(self):
        """
        Read and parse the config file (config.yml). 
        The following keys are read:
        - DATA_PATH: path where processed data will be stored. Default value is the execution dir. 
        """
        try:
            file = open(r'config.yml')
            documents = yaml.full_load(file)
            # DATA PATH
            read_path = documents.get('DATA_PATH')
            if read_path is None:
                warnings.warn('DATA_PATH not found in Config.yml, using local directory')
                self.DATA_DIR = ""
            else:
                self.DATA_DIR = read_path

        except FileNotFoundError:
            warnings.warn('Config.yml not found, using default settings')
            self.DATA_DIR = ""

    def _retrieve_cds_and_merge(self, centre:str, variable:str, start_month:int, lead_time:list):
        """
        Retrieving and procesing seasonal+reanalysis monthly data from the CDS. The function retrieves
        all the available hindcast years using the most recent system. Please, be aware that not all the 
        starting dates and variables are available for all the systems, check the ECMWF wiki for further information: 
        - Starting dates: https://confluence.ecmwf.int/display/CKB/Summary+of+available+data
        - Variables: https://confluence.ecmwf.int/display/CKB/Detailed+list+of+parameters

        This function download the GRIB files for the seasonal forecasts and ERA5, regrid the latter to the 
        seasonal's grid, calculate the seasonal averages (i.e. one sample per year) and then save them
        into the same data structure

        NOTE: this function does not yet support time-ranges crossing the 31th December, (for example DJF)

        Arguments:
            centre (str): Originating centre (one in ['ecmwf', 'meteo_france', 'dwd', 'cmcc'', 'ncep', 'jma'])
            variable (str): one of the variable (e.g. '2m_temperature')
            start_month (int): starting month 1-12
            lead_time (list): list of the lead times that will be averaged to calculate the seasonal average
        """
        # Define start/end year according to hindcast period
        start_year, end_year = self._MODEL_DIC.get(centre).get('range_years')
        # define system
        system = self._MODEL_DIC.get(centre).get('system')
        # check if the time range cross the year
        if any([(start_month+x)>12 for x in lead_time]):
            sys.exit('This function does not support yet forecast periods spanning two years')

        # SEASONAL FORECAST ------------------------------------------------------------------
        request_dict = {
                'originating_centre': centre,
                'variable':variable,
                'product_type':'monthly_mean',
                'year':[x for x in map(str, range(start_year, end_year))],
                'month': str(start_month),
                'leadtime_month':[str(x) for x in lead_time],
                'format':'grib',
                'system': system
            }
        if not self.quiet:
            print('Retrieving seasonal forecasts...')
        c = cdsapi.Client()
        try:
            r = c.retrieve(
                    'seasonal-monthly-single-levels',
                    request_dict) 
        except:
            print(request_dict)
            sys.exit('Exception from cdsapi, check the request or the status of CDS')
        # Update the request and get the reply    
        r.update()
        reply = r.reply

        if reply.get('state') != 'completed':
            print(request_dict)
            # TODO replace sys.exit with raise and check climetlab implementation
            sys.exit('Request not completed, check the request of the status of CDS')
        else:
            r.download('fct_temp_out.grib')

        # OBSERVATIONS ----------------------------------------------------------
        if not self.quiet:
            print('Retrieving ERA5 data...')
        request_dict = {
                'variable':variable,
                'product_type':'monthly_averaged_reanalysis',
                'year':[x for x in map(str, range(start_year, end_year))],
                'month': [x for x in map(str, range(start_month, start_month+lead_time[-1]))],
                'format':'grib',
                'time': '00:00'
            }
        try:
            r = c.retrieve(
                'reanalysis-era5-single-levels-monthly-means',
                request_dict)
        except:
            print(request_dict)
            sys.exit('Exception from cdsapi, check the request of the status of CDS')
            
        # Update the request and get the reply    
        r.update()
        reply = r.reply
        if reply.get('state') != 'completed':
            print(request_dict)
            sys.exit('Request not completed, check the request or the status of CDS')
        else:
            r.download('obs_temp_out.grib')
        
        # Reading GRIB files into xarray Datasets and calculate the seasonal averages
        if not self.quiet:
            print('Processing seasonal forecasts: renaming coordinates and computing annual average')
        fct = xr.open_dataset('fct_temp_out.grib', engine='cfgrib', backend_kwargs=dict(time_dims = ('verifying_time',)))
        fct_final = fct.rename({'latitude':'lat', 'longitude': 'lon'}).groupby('verifying_time.year').mean()
        # Reading ERA5 GRIB file and calculate seasonal average
        if not self.quiet:
            print('Processing ERA5: renaming coordinates and computing annual average')
        obs = xr.open_dataset('obs_temp_out.grib', engine = 'cfgrib')
        # This to remove the spurious months for accumulated variables
        obs = obs.sel(time=np.isin(obs['valid_time.month'], lead_time))
        obs_y = obs.groupby('time.year').mean('time').rename({'latitude':'lat', 'longitude': 'lon'})
        # Regrid ERA5 on forecast's grid
        if not self.quiet:
            print("Regridding ERA5 on forecasts' grid")
        
        regridder = xe.Regridder(obs_y, fct_final, 'bilinear')
        obs_final = regridder(obs_y)
        
        # Check if the variable names and numbers are consistent
        fct_name_var = [x for x in fct_final.data_vars]
        obs_name_var = [x for x in obs_final.data_vars]
        if len(fct_name_var) > 1 or len(obs_name_var) > 1:
        	sys.exit('Error in grib files: variables > 1')
        if fct_name_var[0] != obs_name_var[0]:
            # Manage various cases
            if fct_name_var[0] == 'tprate' and obs_name_var[0] == 'tp':
                print('Convert total precipitation rate to total precipitation')
                fct_final['tprate'] *= 86400
            else:
        	    raise Exception(f'Error in grib files: variables names do not match ({fct_name_var[0]=} != {obs_name_var[0]=}')
        
        # Put the data into the instance structure
        self._data = xr.merge([fct_final.rename({fct_name_var[0]: 'fct_var'}), obs_final.rename({obs_name_var[0]: 'obs_var'})])
        self._data.attrs = {'created': datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
            'hostname': os.uname().nodename,
            'centre': centre,
            'variable': variable
        }

    def get_data(self)-> xr.Dataset:
        """
        Return the xarray Dataset containing seasonal data and ERA5
        """
        if not hasattr(self, '_data'):
            return xr.Dataset()
        else:
            return(self._data)

    def get_basename(self) -> str:
        """
        Return the basename used to save the retrieved data
        """
        return(self.file_out)

    def __init__(self, centre:str, variable:str, start_month:int, lead_time:list, force_download=False, quiet = False):
        """
        Initialise a `C3S_data_loader` object, which represents seasonal+reanalysis monthly data from the CDS. 
        
        The function retrieves the hindcast data using the most recent system. Please, be aware that not all the 
        starting dates and variables are available for all the systems, check the ECMWF wiki for further information: 
        - Starting dates: https://confluence.ecmwf.int/display/CKB/Summary+of+available+data
        - Variables: https://confluence.ecmwf.int/display/CKB/Detailed+list+of+parameters

        This function download the GRIB files for the seasonal forecasts and ERA5, regrid the latter to the 
        seasonal's grid, calculate the seasonal averages (i.e. one sample per year) and then save them
        into the same data structure.

        A unique file name is chosen according to the passed parameters and it is used to save the processed 
        data (saved in the DATA_PATH folder defined in config.yml) in NetCDF format. If the file exists, the constructor will load it
        avoiding to retrieve and process CDS data.

        NOTE: time-ranges crossing the 31th December, (for example DJF), are not yet supported

        Arguments:
            centre (str): Originating centre (one in ['ecmwf', 'meteo_france', 'dwd', 'cmcc'', 'ncep', 'jma'])
            variable (str): one of the variable (e.g. '2m_temperature')
            start_month (int): starting month 1-12
            lead_time (list): list of the lead times that will be averaged to calculate the seasonal average
            force_download (boolean): force the retrieving of the data even if an existing processed file is found
            quiet (boolean): define if the execution should print out some information or not 
        """
        self._MODEL_DIC = {
            'ecmwf':{'system': '5', 'range_years':(1993, 2016), 'start_dates':range(1, 13)}, 
            'ukmo':{'system': 600, 'range_years':(1993, 2016), 'start_dates':[3,4]},
            'meteo_france':{'system': '7', 'range_years':(1993,2016), 'start_dates':range(1, 13)},
            'dwd':{'system': '21', 'range_years':(1993,2016), 'start_dates': [1, 2, 3, 4, 11, 12]},
            'cmcc':{'system': '35', 'range_years':(1993,2016), 'start_dates': [1, 2, 3, 4, 10, 11, 12]},
            'ncep':{'system': '2', 'range_years':(1993,2016), 'start_dates': range(1, 13)},
            'jma':{'system': '2', 'range_years':(1993,2016), 'start_dates': [1, 2, 3, 4, 10, 11, 12]}
            }
        self.quiet = quiet
        self.read_config()
        # File with SEASONAL + OBS
        self.file_out = os.path.join(self.DATA_DIR,  f'{centre}-{variable}-S{start_month}-L{lead_time[0]}-{lead_time[-1]}')
        if not self.quiet:
            print(f"Target filename {self.file_out}")

        if centre not in ['ecmwf', 'ukmo', 'meteo_france', 'dwd', 'cmcc', 'ncep', 'jma']:
            print(f'Centre {centre} not recognised')
            raise ValueError()
        elif start_month not in self._MODEL_DIC.get(centre).get('start_dates'):
            print(f'The starting month {start_month} is not available for {centre} system, please check https://confluence.ecmwf.int/display/CKB/Summary+of+available+data')
            raise ValueError()
        else:
            if force_download or not os.path.exists(self.file_out+'.nc'):
                try:
                    self._retrieve_cds_and_merge(centre, variable, start_month, lead_time)
                    self._data.to_netcdf(self.file_out+'.nc')
                    self._data = xr.open_dataset(self.file_out+'.nc')
                except:
                    print('An error occured while retrieving and processing data from the CDS')
            else:
                if not self.quiet: 
                    print(f"Loading existing file {self.file_out+'.nc'}")
                self._data = xr.open_dataset(self.file_out+'.nc')
            


