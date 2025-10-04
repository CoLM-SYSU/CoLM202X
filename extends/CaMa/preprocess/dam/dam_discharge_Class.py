import calendar
from datetime import datetime
from multiprocessing import Pool, sharedctypes
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import pandas as pd
import matplotlib.dates as mdates
from matplotlib import colors
from scipy.signal import argrelmax
import sys
import netCDF4 as nc

class dam_discharge_Class:
    def __init__(self, namelist):

        self.name     = 'dam_discharge_Class'
        self.version  = '0.1'
        self.release  = '0.1'
        self.date     = 'May 2023'
        self.author   = 'Shulei Zhang, zhangshlei@mail.sysu.edu.cn'   

        #-- read namelist
        self.para_flag = namelist['General'      ]['Para_Tag'   ]
        self.para_flag = False
        self.num_cores = namelist['General'  ]['Num_Cores' ]
        self.outdir    = namelist['General'      ]['Save_Dir'   ]
        self.syear     = namelist['dam_discharge']['Start_Year' ]
        self.eyear     = namelist['dam_discharge']['End_Year'   ]
        self.pyear     = namelist['dam_discharge']['Period_Year']
        self.maxdays   = namelist['dam_discharge']['Max_Days'   ]
        self.simdir    = namelist['dam_discharge']['Sim_Dir'    ]
        self.frc_Qf1   = float(namelist['dam_discharge']['Qf1'        ])
        self.frc_Qf2   = float(namelist['dam_discharge']['Qf2'        ])        


        self.years     = self.eyear - self.syear + 1
        self.year_list = range(self.syear,self.eyear+1)

        self.dam_file    = self.outdir + '/damloc.csv'
        self.output_file = self.outdir + '/damflow.csv'
        self.max_outf    = self.outdir + '/tmp_p01_AnnualMax.bin'
        self.mean_outf   = self.outdir + '/tmp_p01_AnnualMean.bin'
        self.Qmean_file  = self.outdir + '/tmp_p01_AnnualMean.bin'
        self.Q100_file   = self.outdir + '/tmp_p02_100year.bin'
        self.fc_outf     = self.outdir + '/damfcperiod.csv'

        self.damcsv = pd.read_csv(self.dam_file,header=1)
        self.ndams  = len(self.damcsv)
        self.x_arr  = self.damcsv['ix'].values - 1
        self.y_arr  = self.damcsv['iy'].values - 1
        self.x_arr  = self.x_arr.astype(int)
        self.y_arr  = self.y_arr.astype(int)
        print('x_arr of dams: ', self.x_arr)
        print('y_arr of dams: ', self.y_arr)
        print('shape of x_arr: ', self.x_arr.shape)

        self.mean_yeararray  = np.ctypeslib.as_ctypes(np.zeros((self.years, self.ndams), dtype=np.float32))
        self.max_finarray    = np.ctypeslib.as_ctypes(np.zeros((self.years*self.maxdays, self.ndams), dtype=np.float32))
        self.mean_montharray = np.ctypeslib.as_ctypes(np.zeros((self.years, 12, self.ndams), dtype=np.float32))


        self.shared_array_mean_yeararray = sharedctypes.RawArray(self.mean_yeararray._type_ , self.mean_yeararray  )
        self.shared_array_max_finarray   = sharedctypes.RawArray(self.max_finarray._type_   , self.max_finarray    )
        self.shared_array_mean_montharray= sharedctypes.RawArray(self.mean_montharray._type_, self.mean_montharray )

        # self.max_data = np.fromfile(str(self.max_outf), 'float32').reshape(self.years, self.ndams)
        self.finarray = np.zeros((self.ndams))
        self.pps      = self.PlottingPosition(self.years)

        self.inputlist = np.arange(len(self.year_list))

    def read_outflw_p01(self, inp):
        year = self.year_list[inp]

        print(' ')
        print('read natsim outflw: year=', year)

        self.mean_yeararray   = np.ctypeslib.as_array(self.shared_array_mean_yeararray)
        self.max_finarray     = np.ctypeslib.as_array(self.shared_array_max_finarray  )

        ## read NAT outflw-bin
        # outflw_file = simdir + '/outflw' + str(year) + '.bin'
        # outflw_all = np.fromfile(outflw_file, 'float32').reshape(-1,ny,nx)

        ## read NAT outflw-nc
        outflw_file = self.simdir + '/o_outflw' + str(year) + '.nc'
        with nc.Dataset(outflw_file, "r") as cdf:
            outflw_all = cdf.variables["outflw"][:]
        # print(outflw_file)

        ## dam outflw: daily
        outflw_dam = outflw_all[:,self.y_arr,self.x_arr]

        ## dam outflw: annual average
        self.mean_yeararray[inp,:] = np.mean(outflw_dam, axis=0)

        ## dam outflw: annual maximum
        for j, row in self.damcsv.iterrows():
            outflw = outflw_dam[:,j]
            maxindex = argrelmax(outflw, order=8*7)
            maxarray = outflw[maxindex]
            maxarray_sorted = np.sort(maxarray)[::-1]
            if len(maxarray_sorted) > 0:
                self.max_finarray[inp * self.maxdays : (inp+1)*self.maxdays, j] = maxarray_sorted[0:self.maxdays]
            else:
                outflw_sorted = np.sort(outflw)[::-1]
                self.max_finarray[inp * self.maxdays : (inp+1)*self.maxdays, j] = outflw_sorted[0:self.maxdays]

    def read_outflw_opt(self, inp):
        year = self.year_list[inp]

        print(' ')
        print('read natsim outflw: year=', year)

        self.mean_montharray  = np.ctypeslib.as_array(self.shared_array_mean_montharray )

        ## read NAT outflw-bin
        # outflw_file = outdir + '/outflw' + str(year) + '.bin'
        # outflw_all = np.fromfile(outflw_file, 'float32').reshape(-1,ny,nx)

        ## read NAT outflw-nc
        outflw_file = self.simdir + '/o_outflw' + str(year) + '.nc'
        with nc.Dataset(outflw_file, "r") as cdf:
            outflw_all = cdf.variables["outflw"][:]
        # print(outflw_file)

        ## dam outflw: daily
        outflw_dam = outflw_all[:,self.y_arr,self.x_arr]

        ## dam outflw: monthly average
        start_date = f'{year}0101'
        end_date = f'{year}1231'
        dates = pd.date_range(start_date, end_date)

        df = pd.DataFrame(outflw_dam, index=dates)
        df_monthly = df.resample('M').mean()
        self.mean_montharray[inp,:,:] = df_monthly.values

    def PlottingPosition(slef, n):
        alpha = 0.0  #weibull
        ii = np.arange(n)+1
        pp = (ii-alpha)/(n+1-2*alpha)

        return pp

    def func_gum(self, xx):
        n=len(xx)
        b0=np.sum(xx)/n
        j=np.arange(0,n)
        b1=np.sum(j*xx)/n/(n-1)
        lam1=b0
        lam2=2*b1-b0
        aa=lam2/np.log(2)
        cc=lam1-0.5772*aa
        return aa,cc

    def est_gum(self, aa, cc, pp):
        return cc-aa*np.log(-np.log(pp))

    def gum(self, xx, pp, pyear):
        aa,cc=self.func_gum(xx)
        # estimation
        ye=self.est_gum(aa,cc,pp)
        rr=np.corrcoef(xx, ye)[0][1]
        # probable value
        prob=1.0-1.0/pyear
        yp=self.est_gum(aa,cc,prob)
        res=np.array([aa,cc,rr])
        ss='GUM\na={0:8.3f}\nc={1:8.3f}\nr={2:8.3f}\nn={3:4.0f}'.format(res[0],res[1],res[2],len(xx))

        #print(res,ye,yp,ss)
        return res,ye,yp,ss

    def p01_get_annual_discharge(self):

        if self.para_flag:
            p = Pool(self.num_cores)
            res = list(p.map(self.read_outflw_p01, self.inputlist))
            self.mean_yeararray = np.ctypeslib.as_array(self.shared_array_mean_yeararray)
            self.max_finarray = np.ctypeslib.as_array(self.shared_array_max_finarray)
            p.close()
        for inpi in self.inputlist:
            self.read_outflw_p01(inpi)
            self.mean_yeararray = np.ctypeslib.as_array(self.shared_array_mean_yeararray)
            self.max_finarray = np.ctypeslib.as_array(self.shared_array_max_finarray)

        ##------ save data[1]: annual average discharge ------------
        mean_finarray = np.nanmean(self.mean_yeararray, axis=0)
        print(' ')
        print('-- save mean discharge: ', self.mean_outf)
        # print(mean_finarray.shape)
        mean_finarray.astype('float32').tofile(self.mean_outf)

        ##------ save data[2]:flood discharge ------------------
        print(' ')
        print('-- save flood discharge: ', self.max_outf)
        # print(max_finarray.shape)
        self.max_finarray.astype('float32').tofile(self.max_outf)


    def p02_est_100yr_discharge(self):
        #### initial setting -------------------------------------
        outputpath = self.outdir + '/tmp_p02_'+str(self.pyear)+'year.bin'
        self.max_data = np.fromfile(str(self.max_outf), 'float32').reshape(self.years, self.ndams)
        for dam in range(self.ndams):
            site_arr = self.max_data[:, dam]
            if np.max(site_arr) >= 1e+20:
                self.finarray[dam] = np.nan
                continue
            if np.max(site_arr) == np.min(site_arr):
                self.finarray[dam] = np.nan
                continue
            site_arr = np.where(site_arr<0, 0, site_arr)
            site_arr = np.sort(site_arr)
            res, ye, yp, ss = self.gum(site_arr, self.pps, self.pyear)
            #print(str(pyear)+'yr discharge=', yp)

            if yp > 0:
                self.finarray[dam] = yp
            else:
                self.finarray[dam] = np.nan

            # print('damID:', dam+1, ", 100yr discharge:", "{:.1f}".format(yp))

        self.finarray.astype('float32').tofile(outputpath)
        print('file outputted:', outputpath)
        print('###########################################')
        print(' ')

    def p03_complete_discharge(self):
        Qn_all   = np.fromfile(self.Qmean_file, 'float32')
        Q100_all = np.fromfile(self.Q100_file, 'float32')

        damout = pd.DataFrame()
        damout['GRAND_ID'] = self.damcsv['GRAND_ID']
        damout['Qn_CMS'] = Qn_all
        damout['Qf_CMS'] = np.maximum(Q100_all * self.frc_Qf1, Qn_all * self.frc_Qf2)

        # check Qn & Qf
        damout.fillna(-999,inplace=True)
        print(damout)

        ## save output
        damout.to_csv(self.output_file, index=None)
        print(' ')
        print('###############################')
        print('dam parameters:', self.output_file)
        print(' ')

    def opt_dam_fcperiod(self):
        if self.para_flag:
            p = Pool(self.num_cores)
            res = list(p.map(self.read_outflw_opt, self.inputlist))
            self.mean_montharray = np.ctypeslib.as_array(self.shared_array_mean_montharray)
            p.close()
        else:
            for inpi in self.inputlist:
                self.read_outflw_opt(inpi)
                self.mean_montharray = np.ctypeslib.as_array(self.shared_array_mean_montharray)

        ##------ FC period -----------------------
        # calculate FC period: STFC NDFC STOP ------------
        ## mean annual monthly flow
        mean_outflw_dam_mon = np.mean(self.mean_montharray,axis=0)
        mean_value = mean_outflw_dam_mon.mean(axis=0)

        ## calculate fc period
        STFC = mean_outflw_dam_mon.argmin(axis=0)
        NDFC = np.zeros(mean_outflw_dam_mon.shape[1], dtype=int)
        STOP = np.zeros(mean_outflw_dam_mon.shape[1], dtype=int)

        ## calculate fc NDFC & STOP
        mean_outflw_dam_mon_T = mean_outflw_dam_mon.T
        for i, data in enumerate(mean_outflw_dam_mon_T):
            # NDFC
            idx = np.where(data >= mean_value[i])[0]
            if len(idx) > 0:
                for ci in range(0, len(idx)):
                    mi = idx[ci]
                    mi_1 = (idx[ci] - 1) % 11   # previous month
                    if data[mi_1] < mean_value[i]:  # before NDFC < mean_value
                        NDFC[i] = mi + 1
                        break
            else:
                NDFC[i] = -1

            # STOP
            idx = np.where(data <= mean_value[i])[0]
            if len(idx) > 0:
                for ci in range(0, len(idx)):
                    mi = idx[ci]
                    mi_1 = (idx[ci] - 1) % 11   # previous month
                    if data[mi_1] > mean_value[i]:  # before STOP > mean_value
                        STOP[i] = mi + 1
                        break
            else:
                STOP[i] = -1

        ## save FC period
        fc_data = pd.DataFrame()
        fc_data['GRAND_ID'] = self.damcsv['GRAND_ID']
        fc_data['STFC'] = STFC + 1
        fc_data['NDFC'] = NDFC + 1
        fc_data['STOP'] = STOP + 1
        fc_data.fillna(0,inplace=True)
        print(fc_data)
        print('-- save FC period: ', self.fc_outf)
        fc_data.to_csv(self.fc_outf, index=False)
    



    def main_func(self):
        self.p01_get_annual_discharge()
        self.p02_est_100yr_discharge()
        self.p03_complete_discharge()
        self.opt_dam_fcperiod()
