# -*- coding: utf-8 -*-
import calendar
import datetime
import multiprocessing
from multiprocessing import Pool, sharedctypes
import os
import sys
import pandas as pd
import numpy as np
from numpy import ma
import time

# =================================================
# ====  functions for validation  =================
# =================================================
class DamFlow_validation:
   
    def __init__(self, namelist):
        # self.manager = multiprocessing.Manager()
        # super(DamFlow_validation, self).__init__()
        # 
        self.name    = 'DamFlow_validation'
        self.version = '0.1'
        self.release = '0.1'
        self.date    = 'May 2023'
        self.author  = 'Shulei Zhang zhangshlei@mail.sysu.edu.cn'

        self.syear   = namelist['General']['Start_Year' ]
        self.smonth  = namelist['General']['Start_Month']
        self.sday    = namelist['General']['Start_Day'  ]
        self.eyear   = namelist['General']['End_Year'   ]
        self.emonth  = namelist['General']['End_Month'  ]
        self.eday    = namelist['General']['End_Day'    ]
        self.min_len = namelist['General']['Min_len'    ]
        self.otag    = namelist['General']['Out_Tag'    ]
        self.mtag    = namelist['General']['Map_Tag'    ]
        self.ptag    = namelist['General']['Para_Tag'   ]
        self.ncores  = namelist['General']['Num_Cores'  ]
        self.obsdir  = namelist['General']['Obs_Dir'    ]
        self.simdir  = namelist['General']['Sim_Dir'    ]
        self.damlist = namelist['General']['Dam_list'   ]
        self.output  = namelist['General']['Out_Dir'    ]

        self.start_dt = datetime.date(self.syear, self.smonth, self.sday)
        self.last_dt  = datetime.date(self.eyear, self.emonth, self.eday)
        self.last     = (self.last_dt - self.start_dt).days + 1

        dam_file = self.damlist

        #====================================  setting  ====================================
        self.dam_data = pd.read_csv(dam_file,skiprows=[0])
        self.dam_data = self.dam_data.sort_values('GRAND_ID')

        self.dam_id  = np.array(self.dam_data['GRAND_ID'])
        self.dam_name= np.array(self.dam_data['DamName' ])
        self.dam_lon = np.array(self.dam_data['DamLon'  ])
        self.dam_lat = np.array(self.dam_data['DamLat'  ])

        self.ndam    = len(self.dam_id)

        ## save results
        self.res_file = f'./validation_{self.syear}-{self.eyear}.csv'

        self.sim_sto   = np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        self.sim_inflw = np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        self.sim_outflw= np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        self.obs_sto   = np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        self.obs_inflw = np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        self.obs_outflw= np.ctypeslib.as_ctypes(np.full((self.last, self.ndam), -9999.0, dtype=np.float32))
        
        self.shared_array_sim_sto   = sharedctypes.RawArray(self.sim_sto._type_   , self.sim_sto   )
        self.shared_array_sim_inflw = sharedctypes.RawArray(self.sim_inflw._type_ , self.sim_inflw )
        self.shared_array_sim_outflw= sharedctypes.RawArray(self.sim_outflw._type_, self.sim_outflw)
        self.shared_array_obs_sto   = sharedctypes.RawArray(self.obs_sto._type_   , self.obs_sto   )
        self.shared_array_obs_inflw = sharedctypes.RawArray(self.obs_inflw._type_ , self.obs_inflw )
        self.shared_array_obs_outflw= sharedctypes.RawArray(self.obs_outflw._type_, self.obs_outflw)

        self.cval              = np.ctypeslib.as_ctypes(np.full((self.ndam, 6), np.nan, dtype=np.float32))
        self.shared_array_cval = sharedctypes.RawArray(self.cval._type_, self.cval)


    def filter_nan(self, s, o):
        """
        Removes data from simulated and observed data wherever the observed data contains NaNs.
        """
        # Combine the simulated and observed data arrays into a single array
        data = np.column_stack((s.flatten(), o.flatten()))

        # Remove any rows that contain NaN or -9999.0 values
        data[data<0] = np.nan
        data = data[~np.isnan(data).any(axis=1)]

        # Split the remaining data into separate simulated and observed arrays
        s_filtered = data[:, 0]
        o_filtered = data[:, 1]

        return s_filtered, o_filtered

    # ========================================
    def NS(self, s, o):
        """
        Nash Sutcliffe efficiency coefficient
        input:
            s: simulated
            o: observed
        output:
            ns: Nash Sutcliffe efficient coefficient
        """
        # s,o = filter_nan(s,o)
        o = ma.masked_where(o <= 0.0, o).filled(0.0)
        s = ma.masked_where(o <= 0.0, s).filled(0.0)
        o = np.compress(o > 0.0, o)
        s = np.compress(o > 0.0, s)
        s, o = self.filter_nan(s, o)

        ns = 1 - sum((s - o) ** 2) / (sum((o - np.mean(o)) ** 2) + 1e-20)

        return ns

    # ========================================
    def NSlog(self, s, o):
        """
        Nash Sutcliffe efficiency coefficient (log-scale)
        input:
            s: simulated
            o: observed
        output:
            ns: Nash Sutcliffe efficient coefficient
        """
        o = ma.masked_where(o <= 0.0, o).filled(0.0)
        s = ma.masked_where(o <= 0.0, s).filled(0.0)
        o = np.compress(o > 0.0, o)
        s = np.compress(o > 0.0, s)
        s, o = self.filter_nan(s, o)

        s = np.log(s)      # warning!!!
        o = np.log(o)

        nslog = 1 - sum((s - o) ** 2) / (sum((o - np.mean(o)) ** 2) + 1e-20)
        return nslog

    # ========================================
    def KGE(self, s, o):
        """
    	Kling Gupta Efficiency (Kling et al., 2012, http://dx.doi.org/10.1016/j.jhydrol.2012.01.011)
    	input:
            s: simulated
            o: observed
        output:
            KGE: Kling Gupta Efficiency
        """
        o = ma.masked_where(o <= 0.0, o).filled(0.0)
        s = ma.masked_where(o <= 0.0, s).filled(0.0)
        o = np.compress(o > 0.0, o)
        s = np.compress(o > 0.0, s)
        s, o = self.filter_nan(s, o)

        if np.mean(o) == 0 or np.mean(s) == 0:
            kge = np.nan
        else:
            B = np.mean(s) / (np.mean(o) + 1e-20)
            y = (np.std(s) / (np.mean(s)+ 1e-20)) / ((np.std(o) / (np.mean(o) + 1e-20)) + 1e-20)
            r = np.corrcoef(o, s)[0, 1]
            kge = 1 - np.sqrt((r - 1) ** 2 + (B - 1) ** 2 + (y - 1) ** 2)
        # print(B)
        # print(y)
        # print(r)
        return kge

    #==================================== 
    ## read simulation 
    def read_sim_data(self, inputlist):

        yyyy = inputlist
        
        self.sim_sto   = np.ctypeslib.as_array(self.shared_array_sim_sto)
        self.sim_inflw = np.ctypeslib.as_array(self.shared_array_sim_inflw)
        self.sim_outflw= np.ctypeslib.as_array(self.shared_array_sim_outflw)

        s_file = self.simdir + '/damtxt-' + str(yyyy) + '.txt'
        if not os.path.exists(s_file):
            print("no file", s_file)
        else:
            with open(s_file, "r", encoding='gb18030', errors='ignore') as f:  
                lines = f.readlines()
     
            head = 1
            # ------ read output dam id  ------   
            for line in lines[0:head]:
                dam_id_out = line.split()
            # dam_id_out = dam_id_out[1:]
            dam_id_out = np.array(dam_id_out[1:]).astype(int)

            # ------ read dam output data ------  
            dam_data = np.loadtxt(s_file, skiprows=1)
            
            # set time
            if calendar.isleap(yyyy):
                dt = 366
            else:
                dt = 365
            target_dt = datetime.date(yyyy, 1, 1)
            st = (target_dt - self.start_dt).days
            et = st + dt
         
            # Find common id
            common_id = set(self.dam_id).intersection(dam_id_out)
            # Find indices of common numbers in dam_id_out
            idx_out = [i for i, x in enumerate(dam_id_out) if x in common_id]
            # index of data ######## revise according to .txt format
            id_sto=  [x * 3 + 1 for x in idx_out]
            id_inflw=  [x * 3 + 2 for x in idx_out]
            id_outflw=  [x * 3 + 3 for x in idx_out]
            # Find indices of common numbers in dam_id
            idx_id = [i for i, x in enumerate(self.dam_id) if x in common_id]

            # data of yyyy
            sto_temp   = dam_data[:, id_sto]
            inflw_temp = dam_data[:, id_inflw]
            outflw_temp= dam_data[:, id_outflw]
            # all the data 
            self.sim_sto[st:et,idx_id]   = sto_temp
            self.sim_inflw[st:et,idx_id] = inflw_temp
            self.sim_outflw[st:et,idx_id]= outflw_temp

            # print(sto_temp)
    #====================================

    ## read observation
    def read_obs_data(self, inputlist):    ## input dam id
        
        o_file = self.obsdir + '/ResOpsUS_' + str(self.dam_id[inputlist]) + '.csv'

        self.obs_sto   = np.ctypeslib.as_array(self.shared_array_obs_sto   )
        self.obs_inflw = np.ctypeslib.as_array(self.shared_array_obs_inflw )
        self.obs_outflw= np.ctypeslib.as_array(self.shared_array_obs_outflw)

        # read observation 
        head = 1  # header lines      # change according to the obs.txt format -- revise by shulei
        if not os.path.exists(o_file):
            # print("no file", o_file)
            sto = np.ones([self.last], np.float32) * -9999.0
            inflw = np.ones([self.last], np.float32) * -9999.0
            outflw = np.ones([self.last], np.float32) * -9999.0
        else:
            with open(o_file, "r", encoding='gb18030', errors='ignore') as f:  
                lines = f.readlines()
            # ------
            sto_temp = {}
            inflw_temp = {}
            outflw_temp = {}
            for line in lines[head::]:
                line_one = line.split(',')
                line_one = [s.strip() for s in line_one]
                line_one = [s.replace('NA', '-9999.0') for s in line_one]

                yyyymmdd = line_one[0].split('-')
                # print(yyyymmdd)
                if len(yyyymmdd) == 3:
                    yyyy = '%04d' % (int(yyyymmdd[0].replace('"', '')))
                    mm = '%02d' % (int(yyyymmdd[1].replace('"', '')))
                    dd = '%02d' % (int(yyyymmdd[2].replace('"', '')))

                    data_dt = datetime.date(int(yyyy), int(mm), int(dd))
                    # ----- read data and mark date
                    if self.start_dt <= data_dt <= self.last_dt:
                        sto_temp[yyyy + mm + dd] = float(line_one[1])       ## 2column- storage
                        inflw_temp[yyyy + mm + dd] = float(line_one[2])     ## 3column- inflw
                        outflw_temp[yyyy + mm + dd] = float(line_one[3])    ## 4column- outflw
                    elif self.last_dt < data_dt:
                        break            
            # ----- check data date with input date
            sto   = []
            inflw = []
            outflw= []
            for day in np.arange(self.last):
                target_dt = self.start_dt + datetime.timedelta(days=int(day))
                
                yyyy = '%04d' % (target_dt.year)
                mm   = '%02d' % (target_dt.month)
                dd   = '%02d' % (target_dt.day)

                if (yyyy + mm + dd) in sto_temp.keys():
                    sto.append(sto_temp[yyyy + mm + dd])
                    inflw.append(inflw_temp[yyyy + mm + dd])
                    outflw.append(outflw_temp[yyyy + mm + dd])
                else:
                    sto.append(-9999.0)
                    inflw.append(-9999.0)
                    outflw.append(-9999.0)

            self.obs_sto[:,inputlist]   = sto
            self.obs_inflw[:,inputlist] = inflw
            self.obs_outflw[:,inputlist]= outflw

    #=======================================    
    def compare_s_and_o(self, dam_i):
        # move to namelist
        # min_len =  50 # Minimum requirements for validation data [days]

        self.cval = np.ctypeslib.as_array(self.shared_array_cval)    
        #---------
        compare_val = np.zeros((1,6)) * np.nan      
        #---------
        s = self.sim_sto[:,dam_i]
        o = self.obs_sto[:,dam_i]
        [s,o] = self.filter_nan(s, o)

        if len(s) > self.min_len:
            compare_val[0,0] = self.NS(s, o)
            compare_val[0,1] = self.KGE(s, o)
        #---------
        s = self.sim_inflw[:,dam_i]
        o = self.obs_inflw[:,dam_i]
        [s,o] = self.filter_nan(s, o)

        if len(s) > self.min_len:
            compare_val[0,2] = self.NS(s, o)
            compare_val[0,3] = self.KGE(s, o)
        #---------
        s = self.sim_outflw[:,dam_i]
        o = self.obs_outflw[:,dam_i]
        [s,o] = self.filter_nan(s, o)

        if len(s) > self.min_len:
            compare_val[0,4] = self.NS(s, o)
            compare_val[0,5] = self.KGE(s, o)
        #---------
        self.cval[dam_i,:] = compare_val

    def main_function(self):
        start_time = time.time()
        # for parallel calculation
        inputlist = np.arange(self.syear, self.eyear + 1)
        # print(inputlist)

        #===================== step[1] read simulaiton  ==========================
        # read sim_data parallel
        print("-----------------------reading simulation---------------------------")
        for inpi in inputlist:
            res            = self.read_sim_data(inpi)
            self.sim_sto   = np.ctypeslib.as_array(self.shared_array_sim_sto   )
            self.sim_inflw = np.ctypeslib.as_array(self.shared_array_sim_inflw )
            self.sim_outflw= np.ctypeslib.as_array(self.shared_array_sim_outflw)

        # print(sim_sto)
        # print(sim_inflw)
        # print(sim_outflw)

        

        # for parallel calculation
        inputlist = np.arange(self.ndam)
        # print(inputlist)

        #===================== step[2] read observation  ==========================
        # read observation
        print("-----------------------reading observation--------------------------")
        for inpi in inputlist:
            res            = self.read_obs_data(inpi)
            self.obs_sto   = np.ctypeslib.as_array(self.shared_array_obs_sto   )
            self.obs_inflw = np.ctypeslib.as_array(self.shared_array_obs_inflw )
            self.obs_outflw= np.ctypeslib.as_array(self.shared_array_obs_outflw)

        # print(obs_sto)
        # print(obs_inflw)
        # print(obs_outflw)

        #===================== step[3] compare sim vs obs  ==========================
        # for parallel calculation
        inputlist = np.arange(self.ndam)
        # print(inputlist)

        # compare
        print("-----------------------  evaluation   --------------------------")      
        for inpi in inputlist:
            # serial BUG!!!!!
            # res = read_obs_data(inpi)
            res      = self.compare_s_and_o(inpi)
            self.cval= np.ctypeslib.as_array(self.shared_array_cval)

        # save validation for all dam
        df = pd.DataFrame({'GRAND_ID':self.dam_id, 'DamName':self.dam_name, 'DamLon':self.dam_lon, 'DamLat':self.dam_lat,
                           'NSval_sto': self.cval[:,0], 'KGEval_sto': self.cval[:,1], 
                           'NSval_inflw': self.cval[:,2], 'KGEval_inflw': self.cval[:,3],
                           'NSval_outflw': self.cval[:,4], 'KGEval_outflw': self.cval[:,5]})

        # delete dam with no data
        df.dropna(subset=df.columns[-6:], how='all', inplace=True)
        # print(df)
        print(self.output)
        df.to_csv(self.output, na_rep='NA',index=False)


        # end timer
        end_time = time.time()

        # calculate and print elapsed time
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")
