import numpy as np
import sys
import os, sys
import numpy as np
import xarray as xr
from metrics import metrics
import shutil 
# Check the platform

class Validation:
    def __init__(self,info):
        self.name = 'Validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.casedir                 =  info.casedir
        self.variables               =  info.variables
        self.metrics                 =  info.metrics
        self.Max_lat                 =  info.Max_lat
        self.Min_lat                 =  info.Min_lat
        self.Max_lon                 =  info.Max_lon
        self.Min_lon                 =  info.Min_lon
        self.compare_Gres            =  info.compare_Gres

        shutil.rmtree(self.casedir+'/output/',ignore_errors=True)
        os.makedirs(self.casedir+'/output/', exist_ok=True)

        print ('Validation processes starting!')
        print("=======================================")
        print(" ")
        print(" ")

    def make_validation(self):
        # loop the keys in self.variables
        for key in self.variables.keys():
            # loop the keys in self.variables to get the metric output
            for metric in self.metrics.keys():
                variable_value = self.variables[key]                
                o=xr.open_dataset(f'{self.casedir}/tmp/obs/'+f'obs_{key}.nc')[f'{key}']
                s=xr.open_dataset(f'{self.casedir}/tmp/sim/'+f'sim_{variable_value}.nc')[f'{variable_value}']
                if variable_value=='f_fevpa':
                    s=s*86400
                if variable_value=='QVEGT':
                    s=s*86400
                if variable_value=='QVEGE':
                    s=s*86400
                if variable_value=='QSOIL':
                    s=s*86400
                if (variable_value == 'H2OSOI') and (key == 'SMsurf'):
                    s=s[:,0,:,:]

                s['time']=o['time']

                mask1 = np.isnan(s) | np.isnan(o)
                s.values[mask1] = np.nan
                o.values[mask1] = np.nan

                pp2=metrics()
                for metric in self.metrics.keys():
                    if metric == 'pc_bias':
                        pb=pp2.pc_bias(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='pc_bias')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'apb':
                        pb=pp2.apb(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='apb')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'RMSE':
                        pb=pp2.rmse(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='RMSE')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'ubRMSE':
                        pb=pp2.ubRMSE(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='ubRMSE')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'mae':
                        pb=pp2.mae(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='mae')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'bias':
                        pb=pp2.bias(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='bias')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'L':
                        pb=pp2.L(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='L')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                        
                    elif metric == 'correlation':
                        pb=pp2.correlation(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='correlation')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'corrlation_R2':
                        pb=pp2.corrlation_R2(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='corrlation_R2')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                      
                    elif metric == 'NSE':
                        pb=pp2.NSE(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='NSE')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

                    elif metric == 'KGE':
                        pb=pp2.KGE(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='KGE')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
            
                    elif metric == 'KGESS':
                        pb=pp2.KGESS(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='KGESS')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
  
                    elif metric == 'index_agreement':
                        pb=pp2.KGESS(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='index_agreement')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                    elif metric == 'kappa_coeff':
                        pb=pp2.kappa_coeff(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='kappa_coeff')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                    elif metric == 'nBiasScore':
                        pb=pp2.nBiasScore(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='nBiasScore')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                    elif metric == 'nRMSEScore':
                        pb=pp2.nRMSEScore(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='nRMSEScore')
                        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')
                    else:
                        print('No such metric')
                        sys.exit(1)

            print("=======================================")
            print(" ")
            print(" ") 

        return

