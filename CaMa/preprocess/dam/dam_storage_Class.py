
# to estimate flood control voluse from ReGeom data

import multiprocessing
from datetime import datetime
from datetime import date
import os
import numpy as np
import pandas as pd
import sys
from dateutil.relativedelta import relativedelta
import warnings

# ignore FutureWarning messages
warnings.filterwarnings("ignore", category=FutureWarning)

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'


class dam_storage_Class:
   
    def __init__(self, namelist):

        self.name     = 'dam_storage_Class'
        self.version  = '0.1'
        self.release  = '0.1'
        self.date     = 'May 2023'
        self.author   = 'Shulei Zhang, zhangshlei@mail.sysu.edu.cn'   

        #-- read namelist
        self.ptag             =       namelist['General'    ]['Para_Tag'        ]
        self.debug            =       namelist['General'    ]['Debug_Tag'       ]
        self.GRSADdir         =       namelist['dam_storage']['GRSADdir'        ]
        self.ReGeomdir        =       namelist['dam_storage']['ReGeomdir'       ]
        self.ReGeom_ErrorFile =       namelist['dam_storage']['ReGeom_ErrorFile']
        self.savedir          =       namelist['General'    ]['Save_Dir'        ]
        self.ncores           =   int(namelist['General'    ]['Num_Cores'       ])
        self.pc_fld           = float(namelist['dam_storage']['Pc_Fld'          ])
        self.pc_nor           = float(namelist['dam_storage']['Pc_Nor'          ])
        self.pc_con           = float(namelist['dam_storage']['Pc_Con'          ])

        self.damfile          = self.savedir + '/damloc.csv'
        self.outfile          = self.savedir + '/damsto.csv'

        self.grand            = pd.read_csv(self.damfile,header=1)
        self.ndams            = self.grand.shape[0]
        self.error            = pd.read_csv(self.ReGeom_ErrorFile)

  



    def process_dam(self,i):
        gr       = self.grand.iloc[i:i+1]
        nm       = gr['GRAND_ID'].values[0]
        totalsto = gr['CAP_MCM'].values[0] 

        fld_area, fld_sto = np.nan, np.nan
        nor_area, nor_sto = np.nan, np.nan
        con_area, con_sto = np.nan, np.nan

        ## read timeseries file -----
        grsadpath = self.GRSADdir + '/'+ str(nm) + '_intp'
        if not os.path.isfile(grsadpath):
            print('file not found: ' +str(grsadpath)) 
        else:
            df = pd.read_table(grsadpath, index_col=0, parse_dates=True)
            data = df.dropna()
            if np.max(df['3water_enh'].value_counts()) > 12:
                rm_df = df['3water_enh'].value_counts()
                rm_df = rm_df[rm_df > 12]
                rm_df = rm_df.index
                for j in range(len(rm_df)):
                    rm_val = rm_df[j]
                    data['3water_enh'] = data['3water_enh'].replace(rm_val, np.nan)
                data = data.dropna()

            data = data['3water_enh']

            if len(data) < 2:
                print('low data quality: ' +str(grsadpath)) 
            else:
                areamax  = np.max(data)
                fld_area = np.percentile(data, self.pc_fld)
                nor_area = np.minimum(np.percentile(data, self.pc_nor), fld_area*0.9)
                con_area = np.minimum(np.percentile(data, 1),      nor_area*0.9) 
                if self.debug:
                    print('fld_area_org', fld_area,'nor_area_org', nor_area,'con_area_org', con_area)

                ## read reservoir bathymetry data --------------
                regeompath = self.ReGeomdir + '/'+ str(nm) + '.csv'
                if not os.path.isfile(regeompath):
                    print('file not found: ' +str(regeompath)) 
                else:
                    regeom = pd.read_csv(regeompath, header=7)
                    regeom.columns = ['Depth', 'Area', 'Storage']
                    if len(regeom) < 2:
                        print('ReGeom data was empty!!!')
                    else:
                        adj      = regeom['Area'].values[-1] / areamax           
                        fld_area = fld_area * adj
                        nor_area = nor_area * adj
                        con_area = con_area * adj
                        if self.debug:
                            print('fld_area', fld_area, 'areamax', areamax, 'regeom_max', regeom['Area'].values[-1])

                        fld_sto = self.est_sto_by_area(fld_area,regeom,totalsto)                
                        nor_sto = self.est_sto_by_area(nor_area,regeom,totalsto)
                        con_sto = self.est_sto_by_area(con_area,regeom,totalsto) 
    
        ## save data                
        df_i = [nm, totalsto, fld_sto, nor_sto, con_sto, fld_area, nor_area, con_area]

        return df_i




    def est_sto_by_area(self, fld_area,regeom,totalsto):
        # sto_max = 0
        fld_sto = 0
        for i in range(len(regeom)):
            rg = regeom.iloc[i:i+1]
            if rg['Area'].values[0] < fld_area:
                continue
            elif rg['Area'].values[0] == fld_area:
                #fld_sto = rg['Storage'].values[0]
                fld_sto = np.mean(regeom.query('Area == @fld_area')['Storage'])
                sto_max = np.mean(regeom.query('Area == @fld_area')['Storage'])

                adj     = totalsto / regeom['Storage'].values[-1]       
                fld_sto = fld_sto * adj
                # ext_sto = totalsto - fld_sto
                break
            
            elif rg['Area'].values[0] > fld_area:
                rg_p    = regeom.iloc[i-1:i]
                sto_min, area_min = rg_p['Storage'].values[0], rg_p['Area'].values[0]
                sto_max, area_max =   rg['Storage'].values[0],   rg['Area'].values[0]
        
                fld_sto = sto_min + (sto_max - sto_min) * (fld_area - area_min) / (area_max - area_min)
        
                # adj = error_i['V_GRanD_mcm'].values[0] / error_i['V_est_mcm'].values[0]
                adj     = totalsto / regeom['Storage'].values[-1]         
                fld_sto = fld_sto * adj
                # ext_sto = totalsto - fld_sto
                break       

        return fld_sto




    def main_func(self):
        if self.ptag:
            pool = multiprocessing.Pool(self.ncores)
            save_list = pool.map(self.process_dam, range(self.ndams))
        else:
            save_list = []
            for inp in range(self.ndams):
                save_temp = self.process_dam(inp)
                save_list.append(save_temp)
        # save data
        out_vars = ['grand_id', 'totalsto_mcm', 'fldsto_mcm', 'norsto_mcm', 'consto_mcm', 'fldarea', 'norarea','conarea']
        out_data = pd.DataFrame(save_list, columns = out_vars)
        out_data = out_data.sort_values(by=out_data.columns[0])
        print(out_data)
        out_data.to_csv(self.outfile, index=False)
