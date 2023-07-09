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
        self.figplot                 =  info.figplot
        #shutil.rmtree(self.casedir+'/output/',ignore_errors=True)
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
                if key=='E':
                    o=o/86400
                if key=='Et':
                    o=o/86400
                if key=='Ei':
                    o=o/86400
                if key=='Es':
                    o=o/86400
                if (key == 'SMsurf'):
                    s=s[:,0,:,:].squeeze()
                if (key == 'SMroot'):
                    s=s[:,1,:,:].squeeze()

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
                        pb=pp2.RMSE(s,o)
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

                    elif metric == 'correlation_R2':
                        pb=pp2.correlation_R2(s,o)
                        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name='correlation_R2')
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
    def make_plot_index(self):
        from matplotlib import colors
        import matplotlib as mpl

        # read the data
        # loop the keys in self.variables
        for key in self.variables.keys():
            # loop the keys in self.variables to get the metric output
            for metric in self.metrics.keys():
                varfile=(f'{self.casedir}/output/{key}_{metric}.nc')
                if metric == 'pc_bias':
                    vmin=-100.0
                    vmax= 100.0
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'KGE':
                    vmin=-1
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'KGESS':
                    vmin=-1
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'NSE':
                    vmin=-1
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'correlation':
                    vmin=-1
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'correlation_R2':
                    vmin= 0
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    cmap4 = mpl.cm.get_cmap('Spectral_r',)  # coolwarm/bwr/twilight_shifted
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                elif metric == 'index_agreement':
                    vmin=-1
                    vmax= 1
                    bnd = np.linspace(vmin, vmax, 11)
                    cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                    cmap = colors.ListedColormap(cpool)
                    norm = colors.BoundaryNorm(bnd, cmap.N)
                self.plot_map(cmap, norm, key,bnd,metric)
    def plot_map_bk(self,cmap,norm,key,metric):
        from pylab import rcParams
        import cartopy.crs as ccrs
        import matplotlib
        import matplotlib.pyplot as plt
        import cartopy.feature as cfeature

        ### Plot settings
        font = {'family' : 'DejaVu Sans'}
        #font = {'family' : 'Myriad Pro'}
        matplotlib.rc('font', **font)

        params = {'backend': 'ps',
          'axes.labelsize': 12,
          'grid.linewidth': 0.2,
          'font.size': 15,
          'legend.fontsize': 12,
          'legend.frameon': False,
          'xtick.labelsize': 12,
          'xtick.direction': 'out',
          'ytick.labelsize': 12,
          'ytick.direction': 'out',
          'savefig.bbox': 'tight',
          'axes.unicode_minus': False,
          'text.usetex': False}
        rcParams.update(params)

        #set the region of the map based on self.Max_lat, self.Min_lat, self.Max_lon, self.Min_lon
        ds=xr.open_dataset(f'{self.casedir}/output/{key}_{metric}.nc')
        # Extract variables
        lat = ds.lat.values
        lon = ds.lon.values
        var = ds[metric].values
        ''''''
        # Define map projection and set up plot
        proj = ccrs.PlateCarree(central_longitude=0)
        fig, ax = plt.subplots(subplot_kw={'projection': proj})
        ax.set_global()
        # Plot river discharge data on map using pcolormesh function
        im = ax.contourf(lon, lat, var, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

        # Add colorbar to plot and customize it
        cbar = fig.colorbar(im, ax=ax, extend='max')
        cbar.set_label(f'{metric}', rotation=270, labelpad=15)

        # Add other map elements, such as coastlines, rivers, and land
        ax.coastlines()
        ax.add_feature(cfeature.LAND, edgecolor='black')
        ax.add_feature(cfeature.RIVERS)

        # Customize plot as needed
        #ax.set_title('Global River Discharge Map')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        plt.savefig(f'{self.casedir}/output/{key}_{metric}.png',  format='png',dpi=300)
        # Show plot
        plt.show()

        return
    
    def plot_map(self,cmap,norm,key,bnd,metric):
        from pylab import rcParams
        from mpl_toolkits.basemap import Basemap
        import matplotlib
        import matplotlib.pyplot as plt
        ### Plot settings
        font = {'family' : 'DejaVu Sans'}
        #font = {'family' : 'Myriad Pro'}
        matplotlib.rc('font', **font)

        params = {'backend': 'ps',
          'axes.labelsize': 12,
          'grid.linewidth': 0.2,
          'font.size': 15,
          'legend.fontsize': 12,
          'legend.frameon': False,
          'xtick.labelsize': 12,
          'xtick.direction': 'out',
          'ytick.labelsize': 12,
          'ytick.direction': 'out',
          'savefig.bbox': 'tight',
          'axes.unicode_minus': False,
          'text.usetex': False}
        rcParams.update(params)

        #set the region of the map based on self.Max_lat, self.Min_lat, self.Max_lon, self.Min_lon
        ds=xr.open_dataset(f'{self.casedir}/output/{key}_{metric}.nc')
        # Extract variables
        lat = ds.lat.values
        lon = ds.lon.values
        # Reshape lon and lat arrays to be 1-dimensional
        #lon = lon.flatten()
        #lat = lat.flatten()
        lat, lon = np.meshgrid(lat[::-1], lon)
        print(lat.shape)
        print(lon.shape)

        var = ds[metric].transpose("lon", "lat")[:, ::-1].values #.flatten()

        fig = plt.figure()
        #set the region of the map based on self.Max_lat, self.Min_lat, self.Max_lon, self.Min_lon

        M = Basemap(projection='cyl',llcrnrlat=self.Min_lat,urcrnrlat=self.Max_lat,\
                    llcrnrlon=self.Min_lon,urcrnrlon=self.Max_lon,resolution='l')

        
        #fig.set_tight_layout(True)
        #M = Basemap(projection='robin', resolution='l', lat_0=15, lon_0=0)
        M.drawmapboundary(fill_color='white', zorder=-1)
        M.fillcontinents(color='0.8', lake_color='white', zorder=0)
        M.drawcoastlines(color='0.6', linewidth=0.1)
        #M.drawcountries(color='0.6', linewidth=0.1)
       # M.drawparallels(np.arange(-60.,60.,30.), dashes=[1,1], linewidth=0.25, color='0.5')
        #M.drawmeridians(np.arange(0., 360., 60.), dashes=[1,1], linewidth=0.25, color='0.5')
        loc_lon, loc_lat = M(lon, lat)
        #print(loc_lon.shape)
        #var = var.reshape(lon.shape)

        cs = M.contourf(loc_lon, loc_lat,var, cmap=cmap, norm=norm,levels=bnd, extend='both')
        cbaxes = fig.add_axes([0.26, 0.31, 0.5, 0.015])
        cb = fig.colorbar(cs, cax=cbaxes, ticks=bnd, orientation='horizontal', spacing='uniform')
        cb.solids.set_edgecolor("face")
        cb.set_label('%s'%(metric), position=(0.5, 1.5), labelpad=-35)
        plt.savefig(f'{self.casedir}/output/{key}_{metric}.png',  format='png',dpi=300)
        plt.close()
