# -*- coding: utf-8 -*-
import os
import pandas as pd
import xarray as xr
import numpy as np
import sys
import shutil 

class StreamFlow:
    def __init__(self,namelist):
        self.name = 'StreamFlow'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+'/StreamFlow/'

        self.Sim_Dir                 =  namelist['StreamFlowLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['StreamFlowLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['StreamFlowLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['StreamFlowLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['StreamFlowLIST']['Sim_Prefix']
        self.Sim_SpRes               =  namelist['StreamFlowLIST']['Sim_SpRes']
        self.sim_Syear               =  namelist['StreamFlowLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['StreamFlowLIST']['Sim_Eyear']

        #self.FULLLIST                =  namelist['StreamFlowLIST']['FULLLIST'] 
        self.obs_source              =  namelist['StreamFlowLIST']['Obs_source']
        self.OBSDIR                  =  namelist['StreamFlowLIST']['OBSDIR'] 
        self.Pltstn                  =  namelist['StreamFlowLIST']['Pltstn']
        self.Max_UpArea              =  namelist['StreamFlowLIST']['Max_UpArea']
        self.Min_UpArea              =  namelist['StreamFlowLIST']['Min_UpArea']
        if (self.obs_source=='GRDC'):
            if ((namelist['General']['compare_res']=="Hour")):
                print('compare_res="Hour", the compare_res should be "Day","Month" or longer ')
                sys.exit(1) 
        elif (self.obs_source=='GSIM'): 
          #  self.OBSDIR             =  f"./StreamFlow/GSIM"
            if (namelist['General']['compare_res']=="Hour")|(namelist['General']['compare_res']=="Day"):
                print(' the compare_res should be "Month" or longer ')
                sys.exit(1) 
        elif (self.obs_source=='Dai_Trenberth'):
         #   self.OBSDIR             =  f"./StreamFlow/Dai_Trenberth"
            if (namelist['General']['compare_res']=="Hour")|(namelist['General']['compare_res']=="Day"):
                print(' the compare_res should be "Month" or longer ')
                sys.exit(1) 
        elif (self.obs_source=='Hylke'):
        #    self.OBSDIR             =  f"./StreamFlow/Hylke"
            if (namelist['General']['compare_res']=="Hour")|(namelist['General']['compare_res']=="Day"):
                print(' the compare_res should be "Month" or longer ')
                sys.exit(1) 
        else:
            print('The source is not exist.')
            sys.exit(1)

        if (namelist['StreamFlowLIST']['Sim_TimRes'] =="Month"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Month" but compare_res="Hour", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res'] =="Day"):
                print('Sim_TimRes="Month" but compare_res="Day", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['StreamFlowLIST']['Sim_TimRes']=="Year"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Year" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Day"):
                print('Sim_TimRes="Year" but compare_res="Day", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Month"):
                print('Sim_TimRes="Year" but compare_res="Month", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['StreamFlowLIST']['Sim_TimRes']=="Day"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Day" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['StreamFlowLIST']['OBSDIR']    , namelist['StreamFlowLIST']['Pltstn'],namelist['StreamFlowLIST']['Sim_Syear'] , namelist['StreamFlowLIST']['Sim_Eyear']
        del  namelist['StreamFlowLIST']['Sim_Dir']   , namelist['StreamFlowLIST']['Sim_TimRes']  , namelist['StreamFlowLIST']['Sim_DataGroupby']
        del  namelist['StreamFlowLIST']['Sim_Suffix'], namelist['StreamFlowLIST']['Sim_Prefix']  , namelist['StreamFlowLIST']['Sim_SpRes']
        del  namelist['StreamFlowLIST']['Min_UpArea'], namelist['StreamFlowLIST']['Max_UpArea']  , namelist['StreamFlowLIST']['Obs_source']
                 
        self.variables               =  namelist['StreamFlowLIST']
        print(self.variables)
        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}

        if len(self.variables) == 0:
            print("Error: StreamFlowLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)
        #renew case directory
        if (self.obs_source=='GRDC'):
            self.casedir=self.casedir+'/GRDC/'

    def makelist(self):
        #creat list file for streamflow valiation 
        shutil.rmtree(self.casedir,ignore_errors=True)
        #creat tmp directory
        os.makedirs(self.casedir, exist_ok=True)

        if self.Pltstn==True:
            shutil.rmtree (f'{self.casedir}/tmp/plt',ignore_errors=True)
            os.makedirs   (f'{self.casedir}/tmp/plt', exist_ok=True)
        output  =  f"{self.casedir}/selected_list.txt"
        
        if (self.obs_source=='GRDC'):
            self.FULLLIST                =  f"{self.OBSDIR}/list/GRDC_alloc_{self.Sim_SpRes}.txt"
            print(self.FULLLIST)
            station_list                 =  pd.read_csv(f"{self.FULLLIST}",delimiter=r"\s+",header=0)

            station_list['Flag']     =[False] * len(station_list['lon']) #[False for i in range(len(station_list['lon']))] #[False] * len(station_list['lon'])  #(station_list['lon']*0 -9999)*False
            station_list['use_Syear']=[-9999] * len(station_list['lon'])  #int(station_list['lon']*0 -9999)
            station_list['use_Eyear']=[-9999] * len(station_list['lon'])
            station_list['obs_Syear']=[-9999] * len(station_list['lon']) #must be integer
            station_list['obs_Eyear']=[-9999] * len(station_list['lon']) #must be integer

            if (self.compare_res == 'Month'):
                for i in range(len(station_list['ID'])):
                    if(os.path.exists('%s/GRDC_Month/%s_Q_Month.nc'%(self.OBSDIR,station_list['ID'][i]))):
                        with xr.open_dataset('%s/GRDC_Month/%s_Q_Month.nc'%(self.OBSDIR,station_list['ID'][i])) as df:
                                station_list['obs_Syear'].values[i]=df["time.year"].values[0]
                                station_list['obs_Eyear'].values[i]=df["time.year"].values[-1]
                                station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
                                station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
                                if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                    (station_list['lon'].values[i]>=self.Min_lon) &\
                                    (station_list['lon'].values[i]<=self.Max_lon) &\
                                    (station_list['lat'].values[i]>=self.Min_lat) &\
                                    (station_list['lat'].values[i]<=self.Max_lat) &\
                                    (station_list['area1'].values[i]>=self.Min_UpArea) &\
                                    (station_list['area1'].values[i]<=self.Max_UpArea) &\
                                    (station_list['ix2'].values[i] == -9999) 
                                    ): 
                                    station_list['Flag'].values[i]=True
                                    print(f"Station ID : {station_list['ID'].values[i]} is selected")
            elif (self.compare_res == 'Day'):
                for i in range(len(station_list['ID'])):
                    if(os.path.exists('%s/GRDC_Day/%s_Q_Day.Cmd.nc'%(self.OBSDIR,station_list['ID'][i]))):
                        with xr.open_dataset('%s/GRDC_Day/%s_Q_Day.Cmd.nc'%(self.OBSDIR,station_list['ID'][i])) as df:
                                station_list['obs_Syear'].values[i]=df["time.year"].values[0]
                                station_list['obs_Eyear'].values[i]=df["time.year"].values[-1]
                                station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
                                station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
                                if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                    (station_list['lon'].values[i]>=self.Min_lon) &\
                                    (station_list['lon'].values[i]<=self.Max_lon) &\
                                    (station_list['lat'].values[i]>=self.Min_lat) &\
                                    (station_list['lat'].values[i]<=self.Max_lat) &\
                                    (station_list['area1'].values[i]>=self.Min_UpArea) &\
                                    (station_list['area1'].values[i]<=self.Max_UpArea) &\
                                    (station_list['ix2'].values[i] == -9999) 
                                    ): 
                                    station_list['Flag'].values[i]=True
                                    print(f"Station ID : {station_list['ID'].values[i]} is selected")

            ind = station_list[station_list['Flag']==True].index
            data_select = station_list.loc[ind]
            print(data_select)
            if self.Sim_SpRes=="15min":
                lat0=np.arange(89.875,-90,-0.25)
                lon0=np.arange(-179.875,180,0.25)
            elif self.Sim_SpRes=="01min":
                lat0=np.arange(89.9916666666666600,-90,-0.0166666666666667)
                lon0=np.arange(-179.9916666666666742,180,0.0166666666666667)
            elif self.Sim_SpRes=="05min":
                lat0=np.arange(89.9583333333333286,-90,-0.0833333333333333)
                lon0=np.arange(-179.9583333333333428,180,0.0833333333333333)
            elif self.Sim_SpRes=="06min":
                lat0=np.arange(89.95,-90,-0.1)
                lon0=np.arange(-179.95,180,0.1)
            elif self.Sim_SpRes=="03min":
                lat0=np.arange(89.975,-90,-0.05)
                lon0=np.arange(-179.975,180,0.05)
   
            data_select['lon_cama']=[-9999.] * len(data_select['lon']) 
            data_select['lat_cama']=[-9999.] * len(data_select['lat']) 
            for iii in range(len(data_select['ID'])):
                print(iii,len(data_select['ID']))
                data_select['lon_cama'].values[iii]=float(lon0[int(data_select['ix1'].values[iii])-1])
                data_select['lat_cama'].values[iii]=float(lat0[int(data_select['iy1'].values[iii])-1])
                if abs(data_select['lat_cama'].values[iii]-data_select['lat'].values[iii])>1:
                    print(f"Warning: ID {data_select['ID'][iii]} lat is not match")
                if abs(data_select['lon_cama'].values[iii]-data_select['lon'].values[iii])>1:
                    print(f"Warning: ID {data_select['ID'].values[iii]} lon is not match")

            # print(data_select)
            print(f"In total: {len(data_select['ID'])} stations are selected")
            data_select.to_csv(output,index=False)
        elif (self.obs_source=='GSIM'):
            print("Not ready yet");exit()
        elif (self.obs_source=='GSIM'):
            print("Not ready yet");exit()
        elif (self.obs_source=='***8'):
            print("Not ready yet");exit()
        else:
            print("please check the obs_source in namelist");exit()
        return

class SoilMoisture:
    def __init__(self,namelist):
        self.name = 'SoilMoisture'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+'/SoilMoisture/'

        self.Sim_Dir                 =  namelist['SoilMoistureLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['SoilMoistureLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['SoilMoistureLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['SoilMoistureLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['SoilMoistureLIST']['Sim_Prefix']
        self.Sim_SpRes               =  namelist['SoilMoistureLIST']['Sim_SpRes']
        self.sim_Syear               =  namelist['SoilMoistureLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['SoilMoistureLIST']['Sim_Eyear']

        self.obs_source              =  namelist['SoilMoistureLIST']['Obs_source']
        self.OBSDIR                  =  namelist['SoilMoistureLIST']['OBSDIR'] 
        self.Pltstn                  =  namelist['SoilMoistureLIST']['Pltstn']

        if (namelist['SoilMoistureLIST']['Sim_TimRes'] =="Month"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Month" but compare_res="Hour", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res'] =="Day"):
                print('Sim_TimRes="Month" but compare_res="Day", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['SoilMoistureLIST']['Sim_TimRes']=="Year"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Year" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Day"):
                print('Sim_TimRes="Year" but compare_res="Day", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Month"):
                print('Sim_TimRes="Year" but compare_res="Month", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['SoilMoistureLIST']['Sim_TimRes']=="Day"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Day" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['SoilMoistureLIST']['OBSDIR']    , namelist['SoilMoistureLIST']['Pltstn'] 
        del  namelist['SoilMoistureLIST']['Sim_Dir']   , namelist['SoilMoistureLIST']['Sim_TimRes'], namelist['SoilMoistureLIST']['Sim_DataGroupby']
        del  namelist['SoilMoistureLIST']['Sim_Suffix'], namelist['SoilMoistureLIST']['Sim_Prefix'], namelist['SoilMoistureLIST']['Sim_SpRes']
        del  namelist['SoilMoistureLIST']['Sim_Syear'] , namelist['SoilMoistureLIST']['Sim_Eyear'],  namelist['SoilMoistureLIST']['Obs_source']  
        self.variables               =  namelist['SoilMoistureLIST']
        
        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: SoilMoistureLIST is empty!")
            sys.exit()
        else:
            print(self.variables)
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}
        print(self.metrics)

        if (self.obs_source=='ismn'):
            #creat list file for SoilMoisture valiation 
            #creat directory if not exist
            #remove tmp directory if exist
            self.casedir=self.casedir+'/ismn/'
            shutil.rmtree(self.casedir,ignore_errors=True)

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')
        if (self.obs_source=='ismn'):
            #creat tmp directory
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            
            #/tera04/zhwei/validation/SoilMoisture/ismn/
            self.FULLLIST =  f"{self.casedir}/scratch/ismn_{self.Sim_SpRes}_surf.txt"

            # Get a list of all NetCDF files in the folder
            nc_files = [f for f in os.listdir(self.OBSDIR) if f.endswith('.nc')]

            # Create an empty list to store the lon-lat and time information
            data = []

            # Loop through each NetCDF file in the folder
            for file_name in nc_files:
                # Open the NetCDF file using xarray
                ds = xr.open_dataset(os.path.join(self.OBSDIR, file_name))
                for varname in ds.variables:
                    # Check if the variable name contains "sm_0cm_6cm" or "sm_0cm_5cm"
                    if "sm_0cm_6cm" in varname or "sm_0cm_5cm" in varname or "sm_0cm_7cm" in varname or "sm_0cm_10cm" in varname:
                        # If the variable name contains "sm_0cm_6cm" or "sm_0cm_5cm", add the file name and variable to the list    
                        # Get the lon-lat information from the NetCDF file
                        lon = ds['longitude'].values[:]
                        lat = ds['latitude'].values[:]
                        ds=ds.rename({f'{varname}': 'sm_surf'})
                        ds=ds[["sm_surf", "longitude", "latitude"]]
                        # Get the time information from the NetCDF file
                        time = ds['time']
                        #ds=ds.resample(time='1D').pad()
                        ds.to_netcdf(f'{self.casedir}/scratch/{file_name}')
                        SYear = pd.to_datetime(time.values[0]).year
                        SMon = pd.to_datetime(time.values[0]).month
                        EYear = pd.to_datetime(time.values[-1]).year
                        EMon = pd.to_datetime(time.values[-1]).month
                        sitename=file_name[:-3]
                        # Append the data to the list
                        data.append([sitename, float(lon), float(lat), SYear,SMon, EYear, EMon])
                        break
                # Close the NetCDF file
                ds.close()

            # Convert the list of data into a pandas DataFrame
            df = pd.DataFrame(data, columns=['sitename', 'lon', 'lat', 'SYear', 'SMon', 'EYear', 'EMon'])
            # Save the pandas DataFrame to a CSV file
            df.to_csv(self.FULLLIST,index=False)
            print('save list done')
        else:
            print("Error: obs_source is not exist!")
            sys.exit()
        print(" ")
        print(" ")
        print("-------------------------------------------------------------------------------")

        station_list                 = pd.read_csv(f"{self.FULLLIST}",header=0)
        print(station_list)
        station_list['Flag']     =[False] * len(station_list['lon'])  
        station_list['use_Syear']=[-9999] * len(station_list['lon'])   
        station_list['use_Eyear']=[-9999] * len(station_list['lon'])
        station_list['obs_Syear']=station_list['SYear']
        station_list['obs_Eyear']=station_list['EYear']
        station_list['ID']       =station_list['sitename']
        for i in range(len(station_list['lon'])):
            station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
            station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
            if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                station_list['Flag'].values[i]=True

        
        ind = station_list[station_list['Flag']==True].index
        data_select = station_list.loc[ind]
        # print(data_select)
        data_select.to_csv(f"{self.casedir}/selected_list.txt",index=False)

class Transpiration:
    def __init__(self,namelist):
        self.name = 'Transpiration_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+'/Transpiration/'

        self.Sim_Dir                 =  namelist['TranspirationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['TranspirationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['TranspirationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['TranspirationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['TranspirationLIST']['Sim_Prefix']
        self.Sim_SpRes               =  namelist['TranspirationLIST']['Sim_SpRes']
        self.sim_Syear               =  namelist['TranspirationLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['TranspirationLIST']['Sim_Eyear']

        self.obs_source              =  namelist['TranspirationLIST']['Obs_source']
        self.OBSDIR                  =  namelist['TranspirationLIST']['OBSDIR'] 
        self.Pltstn                  =  namelist['TranspirationLIST']['Pltstn']

        if (namelist['General']['compare_res'] !="Day"):
            print('the compare_res should be "Day"!')
            sys.exit(1)
        if (namelist['TranspirationLIST']['Sim_TimRes'] =="Month" or namelist['TranspirationLIST']['Sim_TimRes']=="Year"):
            print('Sim_TimRes="Month" but compare_res="Day", the Sim_TimRes should be "Day" or "Hour"!')
            sys.exit(1) 
      
        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['TranspirationLIST']['OBSDIR']    , namelist['TranspirationLIST']['Pltstn'] 
        del  namelist['TranspirationLIST']['Sim_Dir']   , namelist['TranspirationLIST']['Sim_TimRes'], namelist['TranspirationLIST']['Sim_DataGroupby']
        del  namelist['TranspirationLIST']['Sim_Suffix'], namelist['TranspirationLIST']['Sim_Prefix'], namelist['TranspirationLIST']['Sim_SpRes']
        del  namelist['TranspirationLIST']['Sim_Syear'] , namelist['TranspirationLIST']['Sim_Eyear'] , namelist['TranspirationLIST']['Obs_source']   
        self.variables               =  namelist['TranspirationLIST']
        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: TranspirationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if (self.obs_source=='GLEAM_hybird'):
            #creat list file for Transpiration valiation 
            #creat directory if not exist
            #remove tmp directory if exist
            self.casedir=self.casedir+'/GLEAM_hybird/'
            shutil.rmtree(self.casedir,ignore_errors=True)
            #creat tmp directory
            os.makedirs(self.casedir, exist_ok=True)
            
            #/tera04/zhwei/validation/Transpiration/ismn/
            self.FULLLIST =  f"{self.casedir}/GLEAM_hybird_{self.Sim_SpRes}.txt"
            self.nl       =   f"{self.casedir}/selected_list.txt"
            # Get a list of all NetCDF files in the folder
            nc_files = [f for f in os.listdir(self.OBSDIR) if f.endswith('.nc')]

            # Create an empty list to store the lon-lat and time information
            data = []

            # Loop through each NetCDF file in the folder
            for file_name in nc_files:
                # Open the NetCDF file using xarray
                ds = xr.open_dataset(os.path.join(self.OBSDIR, file_name))
                # Get the lon-lat information from the NetCDF file
                lon = ds['lon'].values
                lat = ds['lat'].values
                # Get the time information from the NetCDF file
                time = ds['time']
                SYear = pd.to_datetime(time.values[0]).year
                SMon = pd.to_datetime(time.values[0]).month
                EYear = pd.to_datetime(time.values[-1]).year
                EMon = pd.to_datetime(time.values[-1]).month
                sitename=file_name[:-3]

                # Append the data to the list
                data.append([sitename, float(lon), float(lat), SYear,SMon, EYear, EMon])
                # Close the NetCDF file
                ds.close()
            #print(data)
            # Convert the list of data into a pandas DataFrame
            df = pd.DataFrame(data, columns=['sitename', 'lon', 'lat', 'SYear', 'SMon', 'EYear', 'EMon'])
            # Save the pandas DataFrame to a CSV file
            df.to_csv(self.FULLLIST,index=False)
            print('save list done')
        else:
            print("Error: obs_source is not exist!")
            sys.exit()
        del df,ds
        print(" ")
        print(" ")
        print("-------------------------------------------------------------------------------")

        station_list                 = pd.read_csv(f"{self.FULLLIST}",header=0)
        station_list['Flag']     =[False] * len(station_list['lon']) #[False for i in range(len(station_list['lon']))] #[False] * len(station_list['lon'])  #(station_list['lon']*0 -9999)*False
        station_list['use_Syear']=[-9999] * len(station_list['lon'])  #int(station_list['lon']*0 -9999)
        station_list['use_Eyear']=[-9999] * len(station_list['lon'])
        station_list['obs_Syear']=station_list['SYear']
        station_list['obs_Eyear']=station_list['EYear']
        station_list['ID']       =station_list['sitename']
        for i in range(len(station_list['lon'])):
            station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
            station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
            if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                station_list['Flag'].values[i]=True

        
        ind = station_list[station_list['Flag']==True].index
        data_select = station_list.loc[ind]
        # print(data_select)
        data_select.to_csv(self.nl,index=False)

class FLUXNET:
    def __init__(self,namelist):
        self.name = 'FLUXNET'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'Mar 2023'
        self.author  = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['FLUXNETLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['FLUXNETLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['FLUXNETLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['FLUXNETLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['FLUXNETLIST']['Sim_Prefix']
        self.sim_Syear               =  namelist['FLUXNETLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['FLUXNETLIST']['Sim_Eyear']

        self.FULLLIST                =  namelist['FLUXNETLIST']['FULLLIST'] 
        self.OBSDIR                  =  namelist['FLUXNETLIST']['OBSDIR'] 
        self.obs_source              =  namelist['FLUXNETLIST']['Obs_source']
        self.Pltstn                  =  namelist['FLUXNETLIST']['Pltstn']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/{self.obs_source}/'

        if (namelist['FLUXNETLIST']['Sim_TimRes'] =="Month"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Month" but compare_res="Hour", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res'] =="Day"):
                print('Sim_TimRes="Month" but compare_res="Day", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['FLUXNETLIST']['Sim_TimRes']=="Year"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Year" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Day"):
                print('Sim_TimRes="Year" but compare_res="Day", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Month"):
                print('Sim_TimRes="Year" but compare_res="Month", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['FLUXNETLIST']['Sim_TimRes']=="Day"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Day" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['FLUXNETLIST']['FULLLIST']  , namelist['FLUXNETLIST']['OBSDIR']      , namelist['FLUXNETLIST']['Pltstn']
        del  namelist['FLUXNETLIST']['Sim_Dir']   , namelist['FLUXNETLIST']['Sim_TimRes']  , namelist['FLUXNETLIST']['Sim_DataGroupby']
        del  namelist['FLUXNETLIST']['Sim_Suffix'], namelist['FLUXNETLIST']['Sim_Prefix']  , namelist['FLUXNETLIST']['Obs_source']

        self.variables               =  namelist['FLUXNETLIST']

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: FLUXNETLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')

        #creat list file for FLUXNET valiation 
        output  =  f"{self.casedir}/selected_list.txt"

        print(" ")
        print(" ")
        print("-------------------------------------------------------------------------------")

        station_list                 = pd.read_csv(self.FULLLIST,header=0)
        print(station_list)
                
        print("-------------------------------------------------------------------------------")
        print(" ")
        print(" ")
    
        station_list['use_Syear']=[-9999] * len(station_list['obs_Syear'])  #int(station_list['lon']*0 -9999)
        station_list['use_Eyear']=[-9999] * len(station_list['obs_Syear'])

        for i in range(len(station_list['ID'])):
            station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
            station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
            if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=(self.Minimum_lenghth-1)) &\
                        (station_list['lon'].values[i]>=self.Min_lon) &\
                        (station_list['lon'].values[i]<=self.Max_lon) &\
                        (station_list['lat'].values[i]>=self.Min_lat) &\
                        (station_list['lat'].values[i]<=self.Max_lat)
                ): 
                station_list['Flag'].values[i]=True
            else:
                station_list['Flag'].values[i]=False
                print("exclude: "+station_list['ID'].values[i],station_list['lon'].values[i],station_list['lat'].values[i],station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i])
        ind = station_list[station_list['Flag']==True].index
        data_select = station_list.loc[ind]
        print(data_select)

        print(len(data_select['ID']))
        data_select.to_csv(output,index=False)

class Evapotranspiration:
    def __init__(self,namelist):
        self.name = 'FLUXNET'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'Mar 2023'
        self.author  = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['EvapotranspirationLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['EvapotranspirationLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['EvapotranspirationLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['EvapotranspirationLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['EvapotranspirationLIST']['Sim_Prefix']
        self.sim_Syear               =  namelist['EvapotranspirationLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['EvapotranspirationLIST']['Sim_Eyear']
        self.Sim_SpRes               =  namelist['EvapotranspirationLIST']['Sim_SpRes']

        self.OBSDIR                  =  namelist['EvapotranspirationLIST']['OBSDIR'] 
        self.obs_source              =  namelist['EvapotranspirationLIST']['Obs_source']
        self.Pltstn                  =  namelist['EvapotranspirationLIST']['Pltstn']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Evapotranspiration/{self.obs_source}/'

        if (namelist['EvapotranspirationLIST']['Sim_TimRes'] =="Month"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Month" but compare_res="Hour", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res'] =="Day"):
                print('Sim_TimRes="Month" but compare_res="Day", the compare_res should be "Month" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['EvapotranspirationLIST']['Sim_TimRes']=="Year"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Year" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Day"):
                print('Sim_TimRes="Year" but compare_res="Day", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
            elif (namelist['General']['compare_res']=="Month"):
                print('Sim_TimRes="Year" but compare_res="Month", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 
        elif (namelist['EvapotranspirationLIST']['Sim_TimRes']=="Day"):
            if (namelist['General']['compare_res']=="Hour"):
                print('Sim_TimRes="Day" but compare_res="Hour", the compare_res should be "Year" or longer than Sim_TimRes!')
                sys.exit(1) 

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['EvapotranspirationLIST']['OBSDIR']      , namelist['EvapotranspirationLIST']['Pltstn']    , namelist['EvapotranspirationLIST']['Sim_SpRes']
        del  namelist['EvapotranspirationLIST']['Sim_Dir']   , namelist['EvapotranspirationLIST']['Sim_TimRes']  , namelist['EvapotranspirationLIST']['Sim_DataGroupby']
        del  namelist['EvapotranspirationLIST']['Sim_Suffix'], namelist['EvapotranspirationLIST']['Sim_Prefix']  , namelist['EvapotranspirationLIST']['Obs_source']

        self.variables               =  namelist['EvapotranspirationLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: EvapotranspirationLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')

        if (self.obs_source=='GLEAM_hybird_PLUMBER2'):
            #creat tmp directory
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            
            self.FULLLIST =  f"{self.casedir}/scratch/{self.obs_source}_{self.Sim_SpRes}.txt"

            # Get a list of all NetCDF files in the folder
            nc_files = [f for f in os.listdir(self.OBSDIR) if f.endswith('.nc')]

            # Create an empty list to store the lon-lat and time information
            data = []
            # Loop through each NetCDF file in the folder
            for file_name in nc_files:
                # Open the NetCDF file using xarray
                ds = xr.open_dataset(os.path.join(self.OBSDIR, file_name))
                # Get the lon-lat information from the NetCDF file
                lon = ds['lon']#.values 
                lat = ds['lat']#.values 
                # Get the time information from the NetCDF file
                time = ds['time']
                ds=ds.resample(time='1D').pad()
                ds.to_netcdf(f'{self.casedir}/scratch/{file_name}')
                SYear = pd.to_datetime(time.values[0]).year
                SMon = pd.to_datetime(time.values[0]).month
                EYear = pd.to_datetime(time.values[-1]).year
                EMon = pd.to_datetime(time.values[-1]).month
                sitename=file_name[:-3]
                # Append the data to the list
                data.append([sitename, float(lon), float(lat), SYear,SMon, EYear, EMon])
                # Close the NetCDF file
                ds.close()
            # Convert the list of data into a pandas DataFrame
            df = pd.DataFrame(data, columns=['sitename', 'lon', 'lat', 'SYear', 'SMon', 'EYear', 'EMon'])
            # Save the pandas DataFrame to a CSV file
            df.to_csv(self.FULLLIST,index=False)
            print('save list done')
        else:
            print("Error: obs_source is not exist!")
            sys.exit()
        print(" ")
        print(" ")
        print("-------------------------------------------------------------------------------")

        station_list                 = pd.read_csv(f"{self.FULLLIST}",header=0)
        print(station_list)
                
        print("-------------------------------------------------------------------------------")
        print(" ")
        print(" ")
        
        station_list['Flag']     =[False] * len(station_list['lon'])  
        station_list['use_Syear']=[-9999] * len(station_list['lon'])   
        station_list['use_Eyear']=[-9999] * len(station_list['lon'])
        station_list['obs_Syear']=station_list['SYear']
        station_list['obs_Eyear']=station_list['EYear']
        station_list['ID']       =station_list['sitename']
        for i in range(len(station_list['lon'])):
            station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
            station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
            if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                station_list['Flag'].values[i]=True

        
        ind = station_list[station_list['Flag']==True].index
        data_select = station_list.loc[ind]
        # print(data_select)
        data_select.to_csv(f"{self.casedir}/selected_list.txt",index=False)

class LAI:
    def __init__(self,namelist):
        self.name = 'LAI'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'Mar 2023'
        self.author  = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['LAILIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['LAILIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['LAILIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['LAILIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['LAILIST']['Sim_Prefix']
        self.sim_Syear               =  namelist['LAILIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['LAILIST']['Sim_Eyear']
        self.Sim_SpRes               =  namelist['LAILIST']['Sim_SpRes']

        self.OBSDIR                  =  namelist['LAILIST']['OBSDIR'] 
        self.obs_source              =  namelist['LAILIST']['Obs_source']
        self.Pltstn                  =  namelist['LAILIST']['Pltstn']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/LAI/{self.obs_source}/'
        
        if (self.obs_source=='Yuan2022'):
            if ((namelist['General']['compare_res']=="Hour") or (namelist['General']['compare_res']=="Year") ):
                print('compare_res="Hour", the compare_res should be "Day","Month" ')
                sys.exit(1) 
            
        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['LAILIST']['OBSDIR']      , namelist['LAILIST']['Pltstn']    , namelist['LAILIST']['Sim_SpRes']
        del  namelist['LAILIST']['Sim_Dir']     , namelist['LAILIST']['Sim_TimRes']  , namelist['LAILIST']['Sim_DataGroupby']
        del  namelist['LAILIST']['Sim_Suffix']  , namelist['LAILIST']['Sim_Prefix']  , namelist['LAILIST']['Obs_source']
        del  namelist['LAILIST']['Sim_Syear']   , namelist['LAILIST']['Sim_Eyear']  

        self.variables               =  namelist['LAILIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: LAILIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')

        if (self.obs_source=='Yuan2022'):
            #creat tmp directory
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            self.FULLLIST =  f"{self.casedir}/scratch/{self.obs_source}_{self.compare_res}.txt"
            # Create an empty list to store the lon-lat and time information
            data = []
            if (self.compare_res=="Day"):
                file_name = 'lai_8-day_500.nc'
                with xr.open_dataset(os.path.join(self.OBSDIR, file_name)) as ds:
                    for ii in range(len(ds['site'])):
                        print(ds['lai'][:,ii])
                        k=ds['lai'][:,ii].squeeze()
                        k.to_netcdf(f"{self.casedir}/scratch/{ds['site'][ii].values}.nc")
                        print(ds['site'][ii].values)
                        lon = ds['lon'][ii].values#.values 
                        lat = ds['lat'][ii].values#.values 
                        # Get the time information from the NetCDF file
                        time = ds['time']
                        SYear = pd.to_datetime(time.values[0]).year
                        SMon = pd.to_datetime(time.values[0]).month
                        EYear = pd.to_datetime(time.values[-1]).year
                        EMon = pd.to_datetime(time.values[-1]).month
                        sitename=ds['site'][ii].values
                        # Append the data to the list
                        data.append([sitename, float(lon), float(lat), SYear,SMon, EYear, EMon])
                        # Close the NetCDF file
                    ds.close()
                # Convert the list of data into a pandas DataFrame
                df = pd.DataFrame(data, columns=['sitename', 'lon', 'lat', 'SYear', 'SMon', 'EYear', 'EMon'])
                # Save the pandas DataFrame to a CSV file
                df.to_csv(self.FULLLIST,index=False)
                print('save list done')
            elif (self.compare_res=="Month"):
                file_name = 'lai_monthly_500.nc'
                with xr.open_dataset(os.path.join(self.OBSDIR, file_name)) as ds:
                    for ii in range(len(ds['site'])):
                        print(ds['lai'][:,ii])
                        k=ds['lai'][:,ii].squeeze()
                        print(ds['site'][ii].values)
                        k.to_netcdf(f"{self.casedir}/scratch/{ds['site'][ii].values}.nc")
                        lon = ds['lon'][ii].values#.values 
                        lat = ds['lat'][ii].values#.values 
                        # Get the time information from the NetCDF file
                        time = ds['time']
                        SYear = pd.to_datetime(time.values[0]).year
                        SMon = pd.to_datetime(time.values[0]).month
                        EYear = pd.to_datetime(time.values[-1]).year
                        EMon = pd.to_datetime(time.values[-1]).month
                        sitename=ds['site'][ii].values
                        # Append the data to the list
                        data.append([sitename, float(lon), float(lat), SYear,SMon, EYear, EMon])
                        # Close the NetCDF file
                    ds.close()
                # Convert the list of data into a pandas DataFrame
                df = pd.DataFrame(data, columns=['sitename', 'lon', 'lat', 'SYear', 'SMon', 'EYear', 'EMon'])
                # Save the pandas DataFrame to a CSV file
                df.to_csv(self.FULLLIST,index=False)
                print('save list done')
        
        station_list                 = pd.read_csv(f"{self.FULLLIST}",header=0)
        print(station_list)
        station_list['Flag']     =[False] * len(station_list['lon'])  
        station_list['use_Syear']=[-9999] * len(station_list['lon'])   
        station_list['use_Eyear']=[-9999] * len(station_list['lon'])
        station_list['obs_Syear']=station_list['SYear']
        station_list['obs_Eyear']=station_list['EYear']
        station_list['ID']       =station_list['sitename']
        for i in range(len(station_list['lon'])):
            station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
            station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
            if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                station_list['Flag'].values[i]=True

        
        ind = station_list[station_list['Flag']==True].index
        data_select = station_list.loc[ind]
        # print(data_select)
        data_select.to_csv(f"{self.casedir}/selected_list.txt",index=False)

class Altimetry:
    def __init__(self,namelist):
        self.name = 'Altimetry'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'Mar 2023'
        self.author  = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['AltimetryLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['AltimetryLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['AltimetryLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['AltimetryLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['AltimetryLIST']['Sim_Prefix']
        self.sim_Syear               =  namelist['AltimetryLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['AltimetryLIST']['Sim_Eyear']
        self.Sim_SpRes               =  namelist['AltimetryLIST']['Sim_SpRes']

        self.OBSDIR                  =  namelist['AltimetryLIST']['OBSDIR'] 
        self.obs_source              =  namelist['AltimetryLIST']['Obs_source']
        self.Pltstn                  =  namelist['AltimetryLIST']['Pltstn']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Altimetry/{self.obs_source}/'

        #remove namelist['FLUXNETLIST']['FULLLIST'],namelist['FLUXNETLIST']['OBSDIR'] , then self.variables is variables list
        del  namelist['AltimetryLIST']['OBSDIR']      , namelist['AltimetryLIST']['Pltstn']    , namelist['AltimetryLIST']['Sim_SpRes']
        del  namelist['AltimetryLIST']['Sim_Dir']   , namelist['AltimetryLIST']['Sim_TimRes']  , namelist['AltimetryLIST']['Sim_DataGroupby']
        del  namelist['AltimetryLIST']['Sim_Suffix'], namelist['AltimetryLIST']['Sim_Prefix']  , namelist['AltimetryLIST']['Obs_source']

        self.variables               =  namelist['AltimetryLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: AltimetryLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')

        if (self.obs_source=='HydroWeb_2.0'):
            #creat tmp directory
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            self.FULLLIST                =  f"{self.OBSDIR}/HydroWeb_VS"
            print(self.FULLLIST)
            station_list                 =  pd.read_csv(f"{self.FULLLIST}",delimiter=',',header=0)
            station_list['ID']=station_list['identifier']
            station_list['Flag']     =[False] * len(station_list['latitude']) #[False for i in range(len(station_list['lon']))] #[False] * len(station_list['lon'])  #(station_list['lon']*0 -9999)*False
            station_list['use_Syear']=[-9999] * len(station_list['latitude'])  #int(station_list['lon']*0 -9999)
            station_list['use_Eyear']=[-9999] * len(station_list['latitude'])
            station_list['obs_Syear']= pd.to_datetime(station_list['start_date']).dt.year #must be integer
            station_list['obs_Eyear']= pd.to_datetime(station_list['end_date']).dt.year #must be integer
            station_list['lon']=station_list['longitude'].astype(float)
            station_list['lat']=station_list['latitude'].astype(float)

            print(station_list)
            # Create an empty list to store the lon-lat and time information
            data = []
            for i in range(len(station_list['ID'])):
                if(os.path.exists(f"{self.OBSDIR}/hydroprd_{station_list['ID'][i]}.txt")):
                    station_list['ID'].values[i]=f"hydroprd_{station_list['ID'][i]}"
                     # read the dataset and specify column names, ignore lines starting with #
                    df = pd.read_csv((f"{self.OBSDIR}/{station_list['ID'][i]}.txt"), header=None, names=['date', 'time', 'height', 'uncertainty'], delim_whitespace=True, comment='#')
                    # combine date and time columns into a single datetime column
                    df['datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'])
                    # set datetime column as index and drop date and time columns
                    height=df['height'].values.reshape(len(df['height']))
                    print(station_list['ID'][i])
                    dss = xr.Dataset({
                    'height': (('time'), height)
                    },
                    coords={
                        'time':  df['datetime']})
                    comp = dict(zlib=True, complevel=6, _FillValue= -9999)
                    encoding = {var: comp for var in dss.data_vars}
                    # creat height attrs
                    dss.height.attrs['Long name']                            = "water level"
                    dss.height.attrs['units']                                = "m"
                    dss.to_netcdf(f"{self.casedir}/scratch/{station_list['ID'][i]}.nc",engine='netcdf4',encoding=encoding)
                    station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
                    station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
                    if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                        station_list['Flag'].values[i]=True

            ind = station_list[station_list['Flag']==True].index
            data_select = station_list.loc[ind]

            
            # print(data_select)
            data_select.to_csv(f"{self.casedir}/selected_list.txt",index=False)

class Dam:
    def __init__(self,namelist):
        self.name = 'Dam (data from Shulei Zhang)'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'Mar 2023'
        self.author  = "Zhongwang Wei / zhongwang007@gmail.com"
        self.Minimum_lenghth         =  namelist['General']['Min_year']
        self.Max_lat                 =  namelist['General']['Max_lat']
        self.Min_lat                 =  namelist['General']['Min_lat']
        self.Max_lon                 =  namelist['General']['Max_lon']
        self.Min_lon                 =  namelist['General']['Min_lon']
        self.Syear                   =  namelist['General']['Syear']
        self.Eyear                   =  namelist['General']['Eyear']
        self.compare_res             =  namelist['General']['compare_res']
        self.casename                =  namelist['General']['casename']

        self.Sim_Dir                 =  namelist['DamLIST']['Sim_Dir']
        self.Sim_TimRes              =  namelist['DamLIST']['Sim_TimRes']
        self.Sim_DataGroupby         =  namelist['DamLIST']['Sim_DataGroupby']
        self.Sim_Suffix              =  namelist['DamLIST']['Sim_Suffix']
        self.Sim_Prefix              =  namelist['DamLIST']['Sim_Prefix']
        self.sim_Syear               =  namelist['DamLIST']['Sim_Syear']
        self.sim_Eyear               =  namelist['DamLIST']['Sim_Eyear']
        self.Sim_SpRes               =  namelist['DamLIST']['Sim_SpRes']

        self.OBSDIR                  =  namelist['DamLIST']['OBSDIR'] 
        self.obs_source              =  namelist['DamLIST']['Obs_source']
        self.Pltstn                  =  namelist['DamLIST']['Pltstn']
        self.casedir                 =  os.path.join(namelist['General']['casedir'], namelist['General']['casename'])+f'/Dam/{self.obs_source}/'

        del  namelist['DamLIST']['OBSDIR'], namelist['DamLIST']['Pltstn']    , namelist['DamLIST']['Sim_SpRes']
        del  namelist['DamLIST']['Sim_Dir']   , namelist['DamLIST']['Sim_TimRes']  , namelist['DamLIST']['Sim_DataGroupby']
        del  namelist['DamLIST']['Sim_Suffix'], namelist['DamLIST']['Sim_Prefix']  , namelist['DamLIST']['Obs_source']
        del  namelist['DamLIST']['Sim_Syear'], namelist['DamLIST']['Sim_Eyear']   

        self.variables               =  namelist['DamLIST']
        print(self.variables)

        #select the key values are string in self.variables, and save them in self.variables
        self.variables = {k: v for k, v in self.variables.items() if isinstance(v, str)}
        if len(self.variables) == 0:
            print("Error: DamLIST is empty!")
            sys.exit()
        #select the key values are True in namelist['metrics'], and save them in self.metrics
        self.metrics = {k: v for k, v in namelist['metrics'].items() if v == True}

    def makelist(self):
        #creat directory if not exist
        if not os.path.isdir(self.casedir):
            #print('The directory is not present. Creating a new one..')
            os.makedirs(self.casedir)
        else:
            print('Be caution: The directory is exist.')
        
        if (self.obs_source=='ResOpsUS'):
            os.makedirs(self.casedir+'/scratch', exist_ok=True)
            self.FULLLIST                =  f"{self.OBSDIR}/../../dam_params_glb_{self.Sim_SpRes}.csv"
            print(self.FULLLIST)
            station_list = pd.read_csv(self.FULLLIST, skiprows=1, header=0,delimiter=",",)
            print(station_list)
            station_list['ID']=str(station_list['GRAND_ID'])
            station_list['Flag']     =[False]  * len(station_list['DamLon']) #[False for i in range(len(station_list['lon']))] #[False] * len(station_list['lon'])  #(station_list['lon']*0 -9999)*False
            station_list['use_Syear']=[-9999]  * len(station_list['DamLon'])  #int(station_list['lon']*0 -9999)
            station_list['use_Eyear']=[-9999]  * len(station_list['DamLon'])
            station_list['obs_Syear']= [-9999] * len(station_list['DamLon'])
            station_list['obs_Eyear']= [-9999] * len(station_list['DamLon'])
            station_list['lon']=station_list['DamLon'].astype(float)
            station_list['lat']=station_list['DamLat'].astype(float)
            data = []
            for i in range(len(station_list['GRAND_ID'])):
                if(os.path.exists(f"{self.OBSDIR}/ResOpsUS_{station_list['GRAND_ID'][i]}.csv")):
                    station_list['ID'].values[i]=f"ResOpsUS_{station_list['GRAND_ID'][i]}"
                     # read the dataset and specify column names, ignore lines starting with #
                    df = pd.read_csv((f"{self.OBSDIR}/{station_list['ID'][i]}.csv"), header=0)
                    df['datetime'] = pd.to_datetime(df['date'])

                    storage=df['storage'].values.reshape(len(df['storage']))
                    inflow=df['inflow'].values.reshape(len(df['inflow']))
                    outflow=df['outflow'].values.reshape(len(df['outflow']))
                    evaporation=df['evaporation'].values.reshape(len(df['evaporation']))
                    elevation=df['elevation'].values.reshape(len(df['elevation']))
                    print(f"deal with {station_list['ID'][i]}")
                    dss = xr.Dataset({
                    'storage': (('time'), storage),
                    'inflow': (('time'), inflow),
                    'outflow': (('time'), outflow),
                    'evaporation': (('time'), evaporation),
                    'elevation': (('time'), elevation)
                    },
                    coords={
                        'time':  df['datetime']})
                    comp = dict(zlib=True, complevel=6, _FillValue= -9999)
                    encoding = {var: comp for var in dss.data_vars}
                    # creat height attrs
                    dss.storage.attrs['Long name']                            = "water level"
                    dss.storage.attrs['units']                                = "m3"
                    dss.inflow.attrs['Long name']                             = "inflow"
                    dss.inflow.attrs['units']                                 = "m3/s"
                    dss.outflow.attrs['Long name']                            = "outflow"
                    dss.outflow.attrs['units']                                = "m3/s"
                    dss.evaporation.attrs['Long name']                        = "evaporation"
                    dss.evaporation.attrs['units']                            = "m3/s"
                    dss.elevation.attrs['Long name']                          = "elevation"
                    dss.elevation.attrs['units']                              = "m"
                    dss.to_netcdf(f"{self.casedir}/scratch/{station_list['ID'][i]}.nc",engine='netcdf4',encoding=encoding)
                    station_list['obs_Syear'].values[i]  = pd.to_datetime(df['datetime'].values[0]).year
                    obs_SMon = pd.to_datetime(df['datetime'].values[0]).month
                    station_list['obs_Eyear'].values[i]  = pd.to_datetime(df['datetime'].values[-1]).year
                    obs_EMon = pd.to_datetime(df['datetime'].values[-1]).month

                    station_list['use_Syear'].values[i]=max(station_list['obs_Syear'].values[i],self.sim_Syear,self.Syear)
                    station_list['use_Eyear'].values[i]=min(station_list['obs_Eyear'].values[i],self.sim_Eyear,self.Eyear)
                    if ((station_list['use_Eyear'].values[i]-station_list['use_Syear'].values[i]>=self.Minimum_lenghth) &\
                                (station_list['lon'].values[i]>=self.Min_lon) &\
                                (station_list['lon'].values[i]<=self.Max_lon) &\
                                (station_list['lat'].values[i]>=self.Min_lat) &\
                                (station_list['lat'].values[i]<=self.Max_lat) 
                                ): 
                        station_list['Flag'].values[i]=True
                        

            ind = station_list[station_list['Flag']==True].index
            data_select = station_list.loc[ind]
            if self.Sim_SpRes=="15min":
                lat0=np.arange(89.875,-90,-0.25)
                lon0=np.arange(-179.875,180,0.25)
            elif self.Sim_SpRes=="01min":
                exit()
            elif self.Sim_SpRes=="05min":
                exit()
            elif self.Sim_SpRes=="06min":
                lat0=np.arange(89.95,-90,-0.1)
                lon0=np.arange(-179.95,180,0.1)
            elif self.Sim_SpRes=="03min":
                lat0=np.arange(89.975,-90,-0.05)
                lon0=np.arange(-179.975,180,0.05)
   
            data_select['lon_cama']=[-9999.] * len(data_select['lon']) 
            data_select['lat_cama']=[-9999.] * len(data_select['lat']) 
            for iii in range(len(data_select['ID'])):
                data_select['lon_cama'].values[iii]=float(lon0[int(data_select['DamIX'].values[iii])-1])
                data_select['lat_cama'].values[iii]=float(lat0[int(data_select['DamIY'].values[iii])-1])
                if abs(data_select['lat_cama'].values[iii]-data_select['lat'].values[iii])>1:
                    print(f"Warning: ID {data_select['ID'][iii]} lat is not match")
                if abs(data_select['lon_cama'].values[iii]-data_select['lon'].values[iii])>1:
                    print(f"Warning: ID {data_select['ID'].values[iii]} lon is not match")

            # print(data_select)
            print(f"In total: {len(data_select['ID'])} stations are selected")
            # print(data_select)
            data_select.to_csv(f"{self.casedir}/selected_list.txt",index=False)




