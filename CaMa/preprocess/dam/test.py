
import os
import read_nml as nml
import sys
import time

#sys.path.append(r"./dam_wuse")
#sys.path.append(r"./dam_basicInfo")
#sys.path.append(r"./dam_discharge")
#sys.path.append(r"./dam_storage")

from dam_basicInfo_Class import dam_basicInfo_Class
from dam_discharge_Class import dam_discharge_Class
from dam_storage_Class import dam_storage_Class
from dam_wuse_Class import dam_wuse_Class


# read namelist
namelist = nml.read_namelist('./dam.nml')



# dam_basicInfo
print("\n")
print("Start dam_basicInfo...")
start_time = time.time()
mainf = dam_basicInfo_Class(namelist)
mainf.main_func()
print("End dam_basicInfo...")
print("--- %s seconds ---" % (time.time() - start_time))



# dam_discharge
print("\n\n\n")
print("Start dam_discharge...")
start_time = time.time()
mainf = dam_discharge_Class(namelist)
mainf.main_func()
print("End dam_discharge...")
print("--- %s seconds ---" % (time.time() - start_time))



# dam_storage
print("\n\n\n")
print("Start dam_storage...")
start_time = time.time()
mainf = dam_storage_Class(namelist)
mainf.main_func()
print("End dam_storage...")
print("--- %s seconds ---" % (time.time() - start_time))
# print start_time




# dam_wuse
print("\n\n\n")
print("Start dam_wuse...")
start_time = time.time()
mainf = dam_wuse_Class(namelist)
p = mainf.main_function()
print("End dam_wuse...")
print("--- %s seconds ---" % (time.time() - start_time))



