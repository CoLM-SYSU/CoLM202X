# -*- coding: utf-8 -*-
import os
import pandas as pd
import xarray as xr
import numpy as np
import sys
import shutil 
from dask.diagnostics import ProgressBar
from joblib import Parallel, delayed

class Makefiles_parallel:
    def __init__(self,namelist,stn_info):
        self.name = 'Makefile'
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
        self.num_cores               =  namelist['General']['num_cores']

        self.casedir                 =  stn_info.casedir
        self.Sim_Dir                 =  stn_info.Sim_Dir
        self.obs_source              =  stn_info.obs_source 
        self.variables               =  stn_info.variables
        self.metrics                 =  stn_info.metrics
        self.Obs_Dir                 =  stn_info.Obs_Dir
        self.Pltstn                  =  stn_info.Pltstn
        self.Sim_Suffix              =  stn_info.Sim_Suffix
        self.Sim_Prefix              =  stn_info.Sim_Prefix
        self.Sim_TimRes              =  stn_info.Sim_TimRes
        self.Sim_SpRes               =  stn_info.Sim_SpRes
  


    def make_obs_parallel(self,station_list,i):
        if (self.obs_source=='GRDC'):
            if (self.compare_res == 'Month'):
                stn='%s/GRDC_Month/%s_Q_Month.nc'%(self.Obs_Dir,station_list['ID'][i])
            elif (self.compare_res == 'Day'):
                stn='%s/GRDC_Day/%s_Q_Day.Cmd.nc'%(self.Obs_Dir,station_list['ID'][i])
        elif (self.obs_source=='FLUXNET'):
            stn1=xr.open_dataset(f"{self.Obs_Dir}"+"/flux/"+f"{station_list['filename'][i]}")
            stn2=xr.open_dataset(f"{self.Obs_Dir}"+"/met/"+f"{station_list['filename'][i]}"[0:-7]+f"Met.nc")
            df = xr.merge([stn2,stn1], compat='override')
            if 'Qle_cor' not in df.data_vars and 'Qle' in df.data_vars:
                df = df.rename({'Qle': 'Qle_cor'})
            df.to_netcdf(f'{self.casedir}/tmp/obs/'+f"obs_{station_list['ID'][i]}.nc")
            stn=f'{self.casedir}/tmp/obs/'+f"obs_{station_list['ID'][i]}.nc"
        elif (self.obs_source=='ismn'):
            stn=f'{self.casedir}/scratch/'+station_list['ID'][i]+'.nc'
        elif (self.obs_source=='GLEAM_hybird'):
            stn=os.path.join(self.Obs_Dir,station_list['ID'][i])+'.nc'
        elif (self.obs_source=='GLEAM_hybird_PLUMBER2'):
            stn=f'{self.casedir}/scratch/'+station_list['ID'][i]+'.nc'
        elif (self.obs_source=='Yuan2022'):
            stn=f'{self.casedir}/scratch/'+station_list['ID'][i]+'.nc'
        elif (self.obs_source=='HydroWeb_2.0'):
            stn=f'{self.casedir}/scratch/'+station_list['ID'][i]+'.nc'
        elif (self.obs_source=='ResOpsUS'):
            stn=f'{self.casedir}/scratch/'+station_list['ID'][i]+'.nc'
        print(stn)
        with xr.open_dataset(stn) as df:
            startx=int(station_list['use_Syear'].values[i])
            endx  =int(station_list['use_Eyear'].values[i])
            dfx = df[self.variables.keys()] 

            dfx1=dfx.sel(time=slice(f'{startx}-01-01',f'{endx}-12-31')) 
            if (self.compare_res == 'Month'):
                dfx2=dfx1.resample(time='1M').mean()
                time_index = pd.date_range(start=f'{startx}-01-01', end=f'{endx}-12-31', freq='M')
            elif (self.compare_res == 'Day'):
                dfx2=dfx1.resample(time='1D').mean()
                time_index = pd.date_range(start=f'{startx}-01-01', end=f'{endx}-12-31', freq='D')
            elif (self.compare_res == 'Hour'):
                dfx2=dfx1.resample(time='1H').mean()
                time_index = pd.date_range(start=f'{startx}-01-01', end=f'{endx}-12-31', freq='H')       
            elif (self.compare_res == 'Year'):
                dfx2=dfx1.resample(time='1Y').mean()
                time_index = pd.date_range(start=f'{startx}-01-01', end=f'{endx}-12-31', freq='Y')   
            else:
                sys.exit(1)
            # Create empty xarray dataset with time index
            ds = xr.Dataset({'data': (['time'], np.nan*np.ones(len(time_index)))},coords={'time': time_index})
            # Reindex original dataset to match new time index
            orig_ds_reindexed = dfx2.reindex(time=ds.time)
            # Merge original and new datasets
            merged_ds = xr.merge([ds, orig_ds_reindexed]).drop_vars('data')
            merged_ds.to_netcdf(f'{self.casedir}/tmp/obs/'+f"obs_{station_list['ID'][i]}"+f"_{station_list['use_Syear'][i]}"+f"_{station_list['use_Eyear'][i]}.nc",engine='netcdf4')
            del startx,endx,dfx,dfx1,dfx2,ds,orig_ds_reindexed,merged_ds,time_index


    def make_sim_parallel(self,simx,station_list,ik):
        startx=int(station_list['use_Syear'].values[ik])
        endx  =int(station_list['use_Eyear'].values[ik])   
        simx['lon'] = xr.where(simx.lon < 0, simx.lon + 360, simx.lon)
        if (self.obs_source=='GRDC'):
            simx1=simx.sel(lat=[station_list['lat_cama'].values[ik]], lon=[station_list['lon_cama'].values[ik]], method="nearest")
        elif (self.obs_source=='ResOpsUS'):
            simx1=simx.sel(lat=[station_list['lat_cama'].values[ik]], lon=[station_list['lon_cama'].values[ik]], method="nearest")
        else:
            simx1=simx.sel(lat=[station_list['lat'].values[ik]], lon=[station_list['lon'].values[ik]], method="nearest")
        simx2=simx1.sel(time=slice(f'{startx}-01-01T00:00:00',f'{endx}-12-31T23:59:59'))
        simx2.to_netcdf(f"{self.casedir}/tmp/sim/sim_{station_list['ID'][ik]}"+f"_{station_list['use_Syear'][ik]}"+f"_{station_list['use_Eyear'][ik]}.nc",engine='netcdf4')
        del simx1,simx2,startx,endx,ik,simx,station_list

    def make_sim_combine_parallel(self,ii):
        VarFiles=(f'{self.Sim_Dir}/{self.Sim_Suffix}{ii}*{self.Sim_Prefix}.nc')
        print(VarFiles)
        with xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time",decode_times=False,preprocess=lambda dfx: dfx[list(self.variables.values())].astype('float32')) as dfx:         
            num=len(dfx['time'])

            if (self.Sim_TimRes=="Hour"):
                dfx['time'] = pd.date_range(f"{ii}-01-01", freq="H", periods=num)
            elif (self.Sim_TimRes=="Day"):
                dfx['time'] = pd.date_range(f"{ii}-01-01", freq="D", periods=num)
            elif (self.Sim_TimRes=="Month"):
                dfx['time'] = pd.date_range(f"{ii}-01-01", freq="M", periods=num)

            if (self.compare_res =="Month"):
                dfx=dfx.resample(time='1M').mean() 
            elif (self.compare_res =="Day"):
                dfx=dfx.resample(time='1D').mean() 
            elif (self.compare_res =="Hour"):
                dfx=dfx.resample(time='1H').mean() 
            else:
                sys.exit(1)     

            dfx=dfx.sel(time=slice(f'{ii}-01-01',f'{ii}-12-31'))
            mask_lon = (dfx.lon >= self.Min_lon) & (dfx.lon <= self.Max_lon)
            mask_lat = (dfx.lat >= self.Min_lat) & (dfx.lat <= self.Max_lat)
            cropped_ds = dfx.where(mask_lon & mask_lat, drop=True)
            delayed_obj=cropped_ds.to_netcdf(f'{self.casedir}/tmp/sim/'+f'sim_{ii}.nc', compute=False)
            with ProgressBar():
                delayed_obj.compute()
            print(f'Year {ii}: Files Combined')
            del  cropped_ds,mask_lon,mask_lat,delayed_obj,num,dfx
        del VarFiles

    def makefiles_parallel(self):
        print("=======================================")
        print("Create directory!")
        print(" ")
        print(" ")
        #remove tmp directory if exist
        shutil.rmtree(f'{self.casedir}/tmp/sim',ignore_errors=True)
        shutil.rmtree(f'{self.casedir}/tmp/obs',ignore_errors=True)
        #creat tmp directory
        os.makedirs(f'{self.casedir}/tmp/sim', exist_ok=True)
        os.makedirs(f'{self.casedir}/tmp/obs', exist_ok=True)
        if self.Pltstn==True:
            shutil.rmtree(f'{self.casedir}/tmp/plt',ignore_errors=True)
            os.makedirs(f'{self.casedir}/tmp/plt', exist_ok=True)

        print(f"tmp directory: {self.casedir}/tmp has been created!")
        print("=======================================")
        print(" ")
        print(" ")

        #read station information
        stnlist  =f"{self.casedir}/selected_list.txt"
        station_list = pd.read_csv(stnlist,header=0)

        
        minyear=min(station_list['use_Syear'].values[:])
        maxyear=max(station_list['use_Eyear'].values[:])
        
        num_cores = self.num_cores #os.cpu_count()  ##用来计算现在可以获得多少cpu核心。 也可以用multipocessing.cpu_count(),或者随意设定<=cpu核心数的数值

        #deal with observation data
        print("=======================================")
        print("deal with observation data")
        print(" ")
        print(" ")
        
        Parallel(n_jobs=num_cores)(delayed(self.make_obs_parallel)(station_list,i) for i in range(len(station_list['ID'])))

        print ('observation data prepared!')
        print("=======================================")
        print(" ")
        print(" ")   


        #deal with simulation data   
        print("=======================================")
        print("deal with simulation data")
        print(" ")
        print(" ") 
        if self.Sim_SpRes=='01min':
            num_cores=1
        # Increase timeout (tune this number to suit your use case).
        timeout=9999999
        Parallel(n_jobs=num_cores,timeout=timeout)(delayed(self.make_sim_combine_parallel)(ii) for ii in range((minyear),(maxyear)+1))
        VarFiles=(f'{self.casedir}/tmp/sim/sim_*.nc')
        print(VarFiles)
        with xr.open_mfdataset(VarFiles, combine='nested',concat_dim="time") as ds1: #,parallel=True,autoclose=True
            delayed_obj=ds1.to_netcdf(f'{self.casedir}/tmp/sim/sim.nc', compute=False)
            with ProgressBar():
                delayed_obj.compute()
        del ds1
        #delete VarFiles if exist
        #shutil.rmtree(f'{self.casedir}/tmp/sim/sim_*.nc',ignore_errors=True)

        with xr.open_dataset(f'{self.casedir}/tmp/sim/sim.nc') as simx:
            Parallel(n_jobs=num_cores,timeout=timeout)(delayed(self.make_sim_parallel)(simx,station_list,ik) for ik in range(len(station_list['ID'])))

        os.remove(f'{self.casedir}/tmp/sim/sim.nc')
        for ii in range((minyear),(maxyear)+1):
            os.remove(f'{self.casedir}/tmp/sim/sim_{ii}.nc')
        del simx
        print ('simulation data prepared!')
        print("=======================================")
        print(" ")
        print(" ")   
        return
    

