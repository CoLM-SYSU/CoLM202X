import shutil
import sys
import multiprocessing
from p01_damflw_validation import DamFlow_validation
from multiprocessing import Process, Queue

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
        #if namelist['General']['Start_Year'] type is not int, report error and exit
        if not isinstance(namelist['General']['Start_Year'], int) or not isinstance(namelist['General']['Start_Year'], int) :
            print('Error: the Start_Year or End_Year type is not int!')
            sys.exit(1)
        if not isinstance(namelist['General']['Start_Month'], int) or not isinstance(namelist['General']['Start_Month'], int) :
            print('Error: the Start_Month or End_Month type is not int!')
            sys.exit(1)
        if not isinstance(namelist['General']['Start_Day'], int) or not isinstance(namelist['General']['End_Day'], int) :
            print('Error: the Start_Day or End_Day type is not int!')
            sys.exit(1)

        #if namelist['General']['Min_year'] type is not float, report error and exit


    except KeyError:
        print('Error: the namelist is not complete!')
        sys.exit(1)
    return namelist

if __name__=='__main__':
    namelist = read_namelist('namelist.txt')
    # print(namelist)
    mainf = DamFlow_validation(namelist)
    p = mainf.main_function()
    # p = multiprocessing.Process(target=mainf.main_function)
    # p.start()
    # p.join()