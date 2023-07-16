# -*- coding: utf-8 -*-
__author__ = "Zhongwang Wei / zhongwang007@gmail.com"
__version__ = "0.1"
__release__ = "0.1"
__date__ = "Mar 2023"
import shutil 
import sys
from Geo_information import get_general_info
from Makefiles_parallel import Makefiles_parallel
from Validation import Validation
def strtobool (val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))
def read_namelist(file_path):
    """
    Read a namelist from a text file.

    Args:
        file_path (str): Path to the text file.

    Returns:
        dict: A dictionary containing the keys and values from the namelist.
    """
    namelist = {}
    
    current_dict = None

    with open(file_path, 'r') as f:
        for line in f:
            # Ignore comments (lines starting with '#')
            if line.startswith('#'):
                continue
            elif not line:
                continue
            # Ignore if line is empty，or only contains whitespace
            elif not line.strip():
                continue
            else:
                # Check if line starts with '&', indicating a new dictionary
                if line.startswith('&'):
                    dict_name = line.strip()[1:]
                    current_dict = {}
                    namelist[dict_name] = current_dict
                # Check if line starts with '/', indicating the end of the current dictionary
                elif line.startswith('/'):
                    current_dict = None
                # Otherwise, add the key-value pair to the current dictionary
                elif current_dict is not None:
                    line = line.split("#")[0]
                    if  not line.strip():
                        continue
                    else:
                        key, value = line.split('=')
                        current_dict[key.strip()] = value.strip()
                        #if the key value in dict is True or False, set its type to bool
                        if value.strip() == 'True' or value.strip() == 'true' or value.strip() == 'False' or value.strip() == 'false':
                            current_dict[key.strip()] = bool(strtobool(value.strip()))
                        #if the key str value in dict is a positive or negative int (), set its type to int
                        elif value.strip().isdigit() or value.strip().replace('-','',1).isdigit():
                            current_dict[key.strip()] = int(value.strip())
                        #if the key str value in dict is a positive or negative float (), set its type to int or float
                        elif value.strip().replace('.','',1).isdigit() or value.strip().replace('.','',1).replace('-','',1).isdigit():
                            current_dict[key.strip()] = float(value.strip())
                        #else set its type to str
                        else:
                            current_dict[key.strip()] = str(value.strip())
    try:
        #if namelist['General']['Syear'] type is not int, report error and exit
        if not isinstance(namelist['General']['Syear'], int) or not isinstance(namelist['General']['Syear'], int) :
            print('Error: the Syear or Eyear type is not int!')
            sys.exit(1)

        #if namelist['General']['Min_year'] type is not float, report error and exit
        if (not isinstance(namelist['General']['Min_year'], float)):
            print('Error: the Min_year type is not float!')
            sys.exit(1)

        if (not isinstance(namelist['General']['Min_year'], float) or not isinstance(namelist['General']['Min_lat'], float) 
            or not isinstance(namelist['General']['Max_lat'], float) or not isinstance(namelist['General']['Min_lon'], float) or not isinstance(namelist['General']['Max_lon'], float)):
            print('Error: the Min_year or Max_lat or Min_lat or Max_lat or Min_lon or Max_lon type is not float!')
            sys.exit(1)

            

    except KeyError:
        print('Error: the namelist is not complete!')
        sys.exit(1)
    return namelist

if __name__=='__main__':
    print("Welcome to the Geo module of the validation system!")
    print("This module is used to validate the Geo information of the model output data")
    print("===============================================================================")
    print("Start running Geo module...")
    
    print("-------------------------------------Caution-----------------------------------")
    print("Please make sure the time axis of the simulation data is consistent with the time axis of the validation data!")
    #input("Press Enter to continue...")
    print("...............................................................................")
    argv                      = sys.argv
    nml                       = str(argv[1])
    namelist                  = read_namelist(f'{nml}')

    def run_validation(module_name, namelist):
        if namelist['General'][module_name]:
            print("Start running %s module..." % module_name)
            module=get_general_info(module_name,namelist)
            print(module.Obs_Dir)

            makefiles = Makefiles_parallel(module)
            makefiles.Makefiles_parallel()

            validation = Validation(module)
            validation.make_validation()
            validation.make_plot_index()

    modules = [
            'Biomass',
            'LAI',
            'BurnedArea',
            'Global_Net_Ecosystem_Carbon_Balance',
            'Gross_Primary_Productivity',
            'Ecosystem_Respiration',
            'Soil_Carbon',
            'Nitrogen_Fixation',

            'Evapotranspiration',
            'Transpiration',
            'Interception',
            'Soil_Evaporation',
            'Soil_Moisture',
            'Runoff',
            'Inundation',
            'Latent_Heat',
            'Sensible_Heat',
            'Terrestrial_Water_Storage_Anomaly',
            'Snow_Water_Equivalent',
            'Permafrost_LIST',

            'Albedo',
            'Surface_Upward_SW_Radiation',
            'Surface_Upward_LW_Radiation',
            'Surface_Net_SW_Radiation',
            'Surface_Net_LW_Radiation',
            'Surface_Net_Radiation',
            'Ground_Heat_Flux',

            'Diurnal_Temperature_Range',
            'Diurnal_Max_Temperature',
            'Diurnal_Min_Temperature',
            'Surface_Downward_SW_Radiation',
            'Surface_Downward_LW_Radiation',
            'Surface_Relative_Humidity',
            'Precipitation',
            'Surface_Air_Temperature'
    ]

    for module in modules:
        if module not in namelist['General']:
            pass
        else:
            run_validation(module, namelist)
