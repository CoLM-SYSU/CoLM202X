# -*- coding: utf-8 -*-
__author__ = "Zhongwang Wei / zhongwang007@gmail.com"
__version__ = "0.1"
__release__ = "0.1"
__date__ = "Mar 2023"
import shutil 
import sys
from Geo_information import Inundation,Evapotranspiration,Transpiration,Interception,SoilEvaporation,Runoff,SoilMoisture,LAI
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
    print("--------------------------------------------=-----------------------------------")

    input("Press Enter to continue...")
    print("...............................................................................")
    namelist = read_namelist('namelist_geo.txt')
    if(namelist['General']['Run_Inundation']):
        p1=Inundation(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()
    if (namelist['General']['Run_Evapotranspiration']):
        p1=Evapotranspiration(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
        p3.make_plot_index()
    if (namelist['General']['Run_Transpiration']):
        p1=Transpiration(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()
    if (namelist['General']['Run_Interception']):
        p1=Interception(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()
    if (namelist['General']['Run_SoilEvaporation']):
        p1=SoilEvaporation(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()
    if (namelist['General']['Run_SoilMoisture']):
        p1=SoilMoisture(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
	    p3.make_plot_index()
    if (namelist['General']['Run_Runoff']):
        p1=Runoff(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()
    if (namelist['General']['Run_LAI']):
        p1=LAI(namelist)
        p2=Makefiles_parallel(p1)
        p2.Makefiles_parallel()
        p3=Validation(p1)
        p3.make_validation()
		p3.make_plot_index()

