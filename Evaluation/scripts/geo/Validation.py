import numpy as np
import sys
import os, sys
import numpy as np
import xarray as xr
from metrics import metrics
import shutil 
# Check the platform

import os
import sys
import shutil
import numpy as np
import xarray as xr
os.environ['PYTHONWARNINGS']='ignore::FutureWarning'
os.environ['PYTHONWARNINGS']='ignore::RuntimeWarning'

class Validation(metrics):
    def __init__(self, info):
        self.name = 'Validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.casedir = info.casedir
        self.variables = info.variables
        self.metrics = info.metrics
        self.Max_lat = info.Max_lat
        self.Min_lat = info.Min_lat
        self.Max_lon = info.Max_lon
        self.Min_lon = info.Min_lon
        self.compare_Gres = info.compare_Gres
        self.figplot = info.figplot
        os.makedirs(self.casedir+'/output/', exist_ok=True)

        print('Validation processes starting!')
        print("=======================================")
        print(" ")
        print(" ")

    def process_metric(self, key, metric, s, o):
        pb = getattr(self, metric)(s, o)
        pb_da = xr.DataArray(pb, coords=[o.lat, o.lon], dims=['lat', 'lon'], name=metric)
        pb_da.to_netcdf(f'{self.casedir}/output/{key}_{metric}.nc')

    def make_validation(self, **kwargs):
        for key, variable_value in self.variables.items():
            o = xr.open_dataset(f'{self.casedir}/tmp/obs/' + f'obs_{key}.nc')[f'{key}']
            s = xr.open_dataset(f'{self.casedir}/tmp/sim/' + f'sim_{variable_value}.nc')[f'{variable_value}']

            if key in ['E', 'Et', 'Ei', 'Es']:
                o = o / 86400

            if key == 'SMsurf':
                s = s[:, 0, :, :].squeeze()

            if key == 'SMroot':
                s = s[:, 1, :, :].squeeze()

            s['time'] = o['time']

            mask1 = np.isnan(s) | np.isnan(o)
            s.values[mask1] = np.nan
            o.values[mask1] = np.nan

            for metric in self.metrics.keys():
                print(metric)
                if hasattr(self, metric):
                    self.process_metric(key, metric, s, o)
                else:
                    print(metric)
                    print('No such metric')
                    sys.exit(1)

            print("=======================================")
            print(" ")
            print(" ")

        return

    
    def make_plot_index(self):
        from matplotlib import colors

        for key in self.variables.keys():
            for metric in self.metrics.keys():
                varfile = f'{self.casedir}/output/{key}_{metric}.nc'
                if metric in ['bias', 'mae', 'ubRMSE', 'apb', 'RMSE', 'L','pc_bias','apb']:
                    vmin = -100.0
                    vmax = 100.0
                elif metric in ['KGE', 'KGESS', 'NSE', 'correlation']:
                    vmin = -1
                    vmax = 1
                elif metric in ['correlation_R2', 'index_agreement','nBiasScore','nRMSEScore']:
                    vmin = 0
                    vmax = 1
                else:
                    vmin = -1
                    vmax = 1

   

                bnd = np.linspace(vmin, vmax, 11)
                cpool = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
                cmap = colors.ListedColormap(cpool)
                norm = colors.BoundaryNorm(bnd, cmap.N)
                self.plot_map(cmap, norm, key,bnd,metric)

    def plot_map(self, colormap, normalize, key, levels, metric, **kwargs):
        # Plot settings
        import matplotlib.pyplot as plt
        import matplotlib
        import numpy as np
        import xarray as xr
        from mpl_toolkits.basemap import Basemap
        from matplotlib import rcParams
        font = {'family': 'DejaVu Sans'}
        matplotlib.rc('font', **font)

        params = {'backend': 'ps',
                  'axes.labelsize': 10,
                  'grid.linewidth': 0.2,
                  'font.size': 12,
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

        # Set the region of the map based on self.Max_lat, self.Min_lat, self.Max_lon, self.Min_lon
        ds=xr.open_dataset(f'{self.casedir}/output/{key}_{metric}.nc')
        # Extract variables
        lat = ds.lat.values
        lon = ds.lon.values
        lat, lon = np.meshgrid(lat[::-1], lon)

        var = ds[metric].transpose("lon", "lat")[:, ::-1].values

        fig = plt.figure()
        M = Basemap(projection='cyl', llcrnrlat=self.Min_lat, urcrnrlat=self.Max_lat,
                    llcrnrlon=self.Min_lon, urcrnrlon=self.Max_lon, resolution='l')

        M.drawmapboundary(fill_color='white', zorder=-1)
        M.fillcontinents(color='0.8', lake_color='white', zorder=0)
        M.drawcoastlines(color='0.6', linewidth=0.1)

        loc_lon, loc_lat = M(lon, lat)

        cs = M.contourf(loc_lon, loc_lat, var, cmap=colormap, norm=normalize, levels=levels, extend='both')
        cbaxes = fig.add_axes([0.26, 0.31, 0.5, 0.015])
        cb = fig.colorbar(cs, cax=cbaxes, ticks=levels, orientation='horizontal', spacing='uniform')
        cb.solids.set_edgecolor("face")
        cb.set_label('%s' % (metric), position=(0.5, 1.5), labelpad=-35)
        plt.savefig(f'{self.casedir}/output/{key}_{metric}.png', format='png', dpi=300)
        plt.close()


