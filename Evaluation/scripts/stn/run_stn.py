# -*- coding: utf-8 -*-
__author__ = "Zhongwang Wei / zhongwang007@gmail.com"
__version__ = "0.1"
__release__ = "0.1"
__date__ = "Mar 2023"
import sys
import sys
from Station_information import FLUXNET, StreamFlow, SoilMoisture, Transpiration, Evapotranspiration, LAI, Altimetry, Dam
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
    print("Welcome to the stn module of the validation system!")
    print("This module is used to validate the station information of the model output data")
    print("===============================================================================")
    print("Start running stn module...")
    
    print("-------------------------------------Caution-----------------------------------")
    print("Please make sure the time axis of the simulation data is consistent with the time axis of the validation data!")
    print("Please make sure the unit of the simulation is consistent with the unit of the validation data!")

    #input("Press Enter to continue...")
    print("...............................................................................")
    argv                      = sys.argv
    nml                       = str(argv[1])
    namelist                  = read_namelist(f'{nml}')
    try:
        if  (namelist['General']['FLUXNET']):
            p1=FLUXNET(namelist)
            p1.makelist()
            ppp1=Makefiles_parallel(namelist,p1)
            ppp1.makefiles_parallel()
            k=Validation(p1.casedir, p1.variables,p1.metrics,p1.Pltstn,p1.Max_lat,p1.Min_lat,p1.Max_lon,p1.Min_lon)
            k.make_validation_P()
            k.make_plot_index()
    except:
        pass
    try:
        if  (namelist['General']['StreamFlow']):
            if namelist['General']['compare_res'] in ["Day", "Month"]:
                pp1=StreamFlow(namelist)
                pp1.makelist()
                ppp1=Makefiles_parallel(namelist,pp1)
                ppp1.makefiles_parallel()
                k1=Validation(pp1.casedir, pp1.variables,pp1.metrics,pp1.Pltstn,pp1.Max_lat,pp1.Min_lat,pp1.Max_lon,pp1.Min_lon)
                k1.make_validation_P()
                k1.make_plot_index()
            else:
                print("Caution: the compare_res is not Day or Month! SoilMoisture validation is not run!")
    except:
        pass
    try:
        if  (namelist['General']['SoilMoisture']):
            if namelist['General']['compare_res'] in ["Day", "Month"]:
                pp2=SoilMoisture(namelist)
                pp2.makelist()
                ppp1=Makefiles_parallel(namelist,pp2)
                ppp1.makefiles_parallel()
                k2=Validation(pp2.casedir, pp2.variables,pp2.metrics,pp2.Pltstn,pp2.Max_lat,pp2.Min_lat,pp2.Max_lon,pp2.Min_lon)
                k2.make_validation_P()
                k2.make_plot_index()
            else:
                print("Caution: the compare_res is not Day or Month! SoilMoisture validation is not run!")
    except:
        pass

    try:
        if  (namelist['General']['Transpiration']):
            if namelist['General']['compare_res']=="Day":
                pp3 =Transpiration(namelist)
                pp3.makelist()
                ppp1=Makefiles_parallel(namelist,pp3)
                ppp1.makefiles_parallel()
                k3  =Validation(pp3.casedir, pp3.variables,pp3.metrics,pp3.Pltstn,pp3.Max_lat,pp3.Min_lat,pp3.Max_lon,pp3.Min_lon)
                k3.make_validation_P()
                k3.make_plot_index()
            else:
                print("Caution: the compare_res is not Day! Transpiration validation is not run!")
    except:
        pass
        
    try:
        if  (namelist['General']['Evapotranspiration']):
            if namelist['General']['compare_res'] in ["Day", "Month"]:
                pp3 =Evapotranspiration(namelist)
                pp3.makelist()
                ppp1=Makefiles_parallel(namelist,pp3)
                ppp1.makefiles_parallel()
                k3  =Validation(pp3.casedir, pp3.variables,pp3.metrics,pp3.Pltstn,pp3.Max_lat,pp3.Min_lat,pp3.Max_lon,pp3.Min_lon)
                k3.make_validation_P()
                k3.make_plot_index()
            else:
                print("Caution: the compare_res is not Day! Evapotranspiration validation is not run!")
    except:
        pass

    try:    
        if  (namelist['General']['LAI']):
            if namelist['General']['compare_res'] in ["Day", "Month"]:
                pp3 =LAI(namelist)
                pp3.makelist()
                ppp1=Makefiles_parallel(namelist,pp3)
                ppp1.makefiles_parallel()
                k3  =Validation(pp3.casedir, pp3.variables,pp3.metrics,pp3.Pltstn,pp3.Max_lat,pp3.Min_lat,pp3.Max_lon,pp3.Min_lon)
                k3.make_validation_P()
                k3.make_plot_index()
            else:
                print("Caution: the compare_res is not Day or Month! LAI validation is not run!")
    except:
        pass
    
    
    try:    
        if  (namelist['General']['Altimetry']):
            if namelist['General']['compare_res']=="Day":
                pp3 =Altimetry(namelist)
                pp3.makelist()
                ppp1=Makefiles_parallel(namelist,pp3)
                ppp1.makefiles_parallel()
                k3  =Validation(pp3.casedir, pp3.variables,pp3.metrics,pp3.Pltstn,pp3.Max_lat,pp3.Min_lat,pp3.Max_lon,pp3.Min_lon)
                k3.make_validation_P()
                k3.make_plot_index()
            else:
                print("Caution: the compare_res is not Day ! Altimetry validation is not run!")
    except:
        pass
    
    try:
        if  (namelist['General']['Dam']):
            if namelist['General']['compare_res'] in ["Day", "Month"]:
                pp3 =Dam(namelist)
                pp3.makelist()
                ppp1=Makefiles_parallel(namelist,pp3)
                ppp1.makefiles_parallel()
                k3  =Validation(pp3.casedir, pp3.variables,pp3.metrics,pp3.Pltstn,pp3.Max_lat,pp3.Min_lat,pp3.Max_lon,pp3.Min_lon)
                k3.make_validation_P()
                k3.make_plot_index()
            else:
                print("Caution: the compare_res is not Day or Month! Dam validation is not run!")
    except:
        pass
    



