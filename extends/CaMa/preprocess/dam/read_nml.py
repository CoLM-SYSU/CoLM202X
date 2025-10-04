import sys

def strtobool (val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1', 'ture','ok'):
        return True
    elif val in ('n', 'no', 'f', 'false', 'off', '0','flase','no'):
        return False
    else:
        raise ValueError("invalid truth value %r" % (val,))

def echo_error():
    '''
    Print message to stdout
    '''
    print("\n-"*40)
    print("\nError: If you updated the namelist, then please also update the parameter\
           attributes in the nml_info dictionary in read_nml")


def convert_type(dict_name, key, value, file_path):
    '''
    Convert the type of value to the type of key in dict_name
    '''
    nml_info = {
        'General': {
            'Para_Tag'        : 'bool',
            'Debug_Tag'       : 'bool',
            'Num_Cores'       : 'int',
            'Map_Tag'         : 'str',
            'Map_Dir'         : 'str',
            'Save_Dir'        : 'str',
        },

        'dam_basicInfo': {
            'Min_Error'       : 'float',
            'Min_Uparea'      : 'float',
            'GRanD_If'        : 'str',
        },

        'dam_discharge': {
            'Start_Year'      : 'int',
            'End_Year'        : 'int',
            'Period_Year'     : 'int',
            'Max_Days'        : 'int',
            'Qf1'             : 'float',
            'Qf2'             : 'float',
            'Sim_Dir'         : 'str',
        },

        'dam_storage': {
            'Pc_Fld'          : 'float',
            'Pc_Nor'          : 'float',
            'Pc_Con'          : 'float',
            'GRSADdir'        : 'str',
            'ReGeomdir'       : 'str',
            'ReGeom_ErrorFile': 'str',
        },
    }
    if dict_name not in nml_info:
        print("-"*40)
        print("\tError: the dictionary name is not supported!")
        print("\tThe supported dictionary names are:", nml_info.keys())
        print("\tThe input dictionary name is:", dict_name)
        print("\tThe input namelist file :", file_path)
        echo_error()
        sys.exit(1)
    if key not in nml_info[dict_name]:
        print("-"*40)
        print("\tError: the parameter name is not supported!")
        print("\tThe input dictionary name is:", dict_name)
        print("\tThe supported parameter names are:", nml_info[dict_name].keys())
        print("\tThe input parameter name is:", dict_name)
        print("\tThe input namelist file :", file_path)
        echo_error()
        sys.exit(1)


    if nml_info[dict_name][key] == 'bool':
        return strtobool(value)
    elif nml_info[dict_name][key] == 'int':
        return int(value)
    elif nml_info[dict_name][key] == 'float':
        return float(value)
    elif nml_info[dict_name][key] == 'str':
        return value
    else:
        print("Error: the key value type is not supported!")
        sys.exit(1)



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
                        current_dict[key.strip()] = convert_type(dict_name, key.strip(), value.strip(), file_path)

    try:
        print("\nAll the namelist parameters are read successfully!")


    except KeyError:
        print('\nError: the namelist is not complete!\n')
        sys.exit(1)
    return namelist

