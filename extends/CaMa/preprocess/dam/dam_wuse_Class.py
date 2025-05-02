'''
Script to idenfity the WaterUse_grids for dams at diverse resolutions
Script to check whether the of different dams overlaps. If so, the grid will be divided equally to each dam.
Apr 25, 2023, written by Shulei Zhang

'''
import glob
import multiprocessing
from multiprocessing import Pool, sharedctypes
import sys
import time
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import os



class dam_wuse_Class:
   
    def __init__(self, namelist):

        self.name     = 'dam_wuse_Class'
        self.version  = '0.1'
        self.release  = '0.1'
        self.date     = 'May 2023'
        self.author   = 'Shulei Zhang, zhangshlei@mail.sysu.edu.cn'

        #-- read namelist
        self.ptag     = namelist['General']['Para_Tag'  ]
        self.ncores   = namelist['General']['Num_Cores' ]
        self.mtag     = namelist['General']['Map_Tag'   ]
        self.mapdir   = namelist['General']['Map_Dir'   ]
        self.outdir   = namelist['General']['Save_Dir'  ]

        self.damlist  = self.outdir + '/damloc.csv'
        self.wuse_dir = self.outdir + '/WaterUse_grids'

        #-- file to save reslts
        if not os.path.isdir(self.wuse_dir):
            os.mkdir(self.wuse_dir)
        self.ix_file  = f'{self.wuse_dir}/ix_{self.mtag}.txt'
        self.iy_file  = f'{self.wuse_dir}/iy_{self.mtag}.txt'
        self.share_file = f'{self.wuse_dir}/grid_share_{self.mtag}.txt'

        #-- read dam data
        self.dam_data = pd.read_csv(self.damlist,skiprows=[0])
        self.dam_data = self.dam_data.sort_values('GRAND_ID')
        self.GRAND_ID= np.array(self.dam_data['GRAND_ID'   ])
        self.dam_name= np.array(self.dam_data['DAM_NAME'   ])
        self.dam_ix  = np.array(self.dam_data['ix'         ]).astype(int)
        self.dam_iy  = np.array(self.dam_data['iy'         ]).astype(int)
        self.dam_vol = np.array(self.dam_data['CAP_MCM'    ])

        self.ndam    = len(self.GRAND_ID)

        self.dam_id     = []
        self.data_area1 = []
        self.data_area2 = []
        self.data_col   = []
        self.data_ix    = []
        self.data_iy    = []
        self.elevtn  = None
        self.grdare  = None
        self.nx      = None
        self.ny      = None
        self.gsize   = None
        self.west    = None
        self.east    = None
        self.south   = None
        self.north   = None
        self.wuse_area = None
        self.max_column  = None
        self.ix_lim  = None
        self.iy_lim  = None



    def volToWUSE_area(self,vol):
        '''
        function to define the relationship between TotalVol_mcm and wateruse_area accodring to 《中国大型水库标准》中库容与灌溉面积的关系
        '''
        # vol = self.dam_vol
        ref_vol = np.array( [0.1,  1,  10,  100, 1000, 10000, 100000, 1000000])     # unit: mcm 
        ref_area = np.array([  2, 20, 200, 2000, 6000, 12000,  18000,   24000])     # unit: km2
        f_area = interp1d(ref_vol, ref_area)
        area = f_area(vol)
        self.wuse_area = area



    def read_map_data(self, map_dir):
        '''
        function to read map data
        '''
        print("Read Map Files: /params.bin") 
        fparam = f'{map_dir}/params.txt'
        with open(fparam, 'r') as f:
            self.nx = int((f.readline()).split()[0])
            self.ny = int((f.readline()).split()[0])
            f.readline()
            self.gsize = float((f.readline()).split()[0])
            self.west = float((f.readline()).split()[0])
            self.east = float((f.readline()).split()[0])
            self.south = float((f.readline()).split()[0])
            self.north = float((f.readline()).split()[0])

        print("Read Map Files: /grdare.bin") 
        finp = f"{map_dir}/grdare.bin"   # unit: m2
        grdare = np.zeros((self.nx, self.ny), dtype=np.float32)
        with open(finp, "rb") as f:
            # read grdare
            f.seek(0)
            grdare_bytes = f.read(4 * self.nx * self.ny)
            grdare = np.frombuffer(grdare_bytes, dtype=np.float32)
            grdare = np.reshape(grdare, (self.nx, self.ny), order='F')
        self.grdare = grdare * 10**(-6) # unit: m2 to km2

        print("Read Map Files: /elevtn.bin") 
        finp = f"{map_dir}/elevtn.bin"     # unit: m
        elevtn = np.zeros((self.nx, self.ny), dtype=np.float32)
        with open(finp, "rb") as f:
            # read elevtn
            f.seek(0)
            elevtn_bytes = f.read(4 * self.nx * self.ny)
            elevtn = np.frombuffer(elevtn_bytes, dtype=np.float32)
            self.elevtn = np.reshape(elevtn, (self.nx, self.ny), order='F')


    def cal_wuse_grid(self, dam_i):
        '''
        function to calculate the wateruse grid (ix,iy) for a given dam
        step[1]: Set the initial search grid
        step[2]: Keep only the downstream area (elevation below that of dam)
        step[3]: accumulate grid area to minimize the difference between acc_area and wuse_area
        step[4]: save ix and iy of the wateruse grids to .txt file    
        '''
        # print(f'dam id: {dam_id[dam_i]}-{dam_i}')
        
        #-- dam information
        ix = int(self.dam_ix[dam_i]-1)  # index_ix of the dam in the map
        iy = int(self.dam_iy[dam_i]-1)  # index_iy of the dam in the map

        dam_elev = self.elevtn[ix,iy]  # elevation of the dam
        dam_area = self.wuse_area[dam_i]  # wateruse area of the dam
        # print(ix)
        # print(iy)
        # print(elevtn[ix,iy])
        # print(grdare[ix,iy])

        #--  
        save_grid_ix = [ix+1]       # to save ix+1
        save_grid_iy = [iy+1]       # to save iy+1
        save_grid_area = [self.grdare[ix,iy]]       # to save grid area 
        save_diff_area = [abs(self.grdare[ix,iy]-dam_area)]       # to save area difference
        acc_area = self.grdare[ix,iy]  # accumulative grid area 
        #--
        di = 0 # number of searches
        ix_search_last = []
        iy_search_last = []
        while (acc_area < dam_area):
            di = di + 1
            #=============== step[1]: Set the initial search grid =============
            #-- Search the outer grids centered on the dam
            ix_search = np.arange(ix-di,ix+di+1)  # ix-di ix ix+di
            iy_search = np.arange(iy-di,iy+di+1)  # iy-di iy iy+di
            grid_search = np.meshgrid(ix_search,iy_search)
            ix_search = grid_search[0]#.flatten()
            iy_search = grid_search[1]#.flatten()
            #-- delete the grids that have been searched
            loc = (~np.isin(ix_search,ix_search_last)+~np.isin(iy_search,iy_search_last))
            #-- save the last outer grids and to delete later
            ix_search_last = ix_search
            iy_search_last = iy_search
            #-- delete the grids that have been searched
            ix_search = ix_search[loc]
            iy_search = iy_search[loc]
            # print('--delete grids that have been searched--')
            # print(ix_search)
            # print(iy_search)
            #-- delete the grids out of the map
            loc = (ix_search<=self.nx-1)*(ix_search>=0)*(iy_search<=self.ny-1)*(iy_search>=0)
            ix_search = ix_search[loc]
            iy_search = iy_search[loc] 
            # print('--delete grids out of the map--')
            # print(ix_search)
            # print(iy_search)

            #=============== step[2]: Keep only the downstream area (elevation below that of dam) =============
            check_ele = self.elevtn[ix_search,iy_search]  # elevation of all the grids
            ix_search = ix_search[np.where(check_ele<dam_elev)]
            iy_search = iy_search[np.where(check_ele<dam_elev)]
            check_ele = check_ele[np.where(check_ele<dam_elev)]
            # print(ix_search)
            # print(iy_search)
            # print(check_ele)
            #-- rearrange ix_search in order of check_ele value from small to large !!! 
            ix_search = ix_search[np.argsort(check_ele)]
            iy_search = iy_search[np.argsort(check_ele)]
            check_ele = check_ele[np.argsort(check_ele)]
            # print(ix_search)
            # print(iy_search)
            # print(check_ele)

            #=============== step[3]: accumulate grid area to minimize the difference between acc_area and wuse_area =============
            num_search = ix_search.shape[0]  #  number of grids to be searched
            for i in range(num_search):
                grid_area = self.grdare[ix_search[i],iy_search[i]]
                if grid_area > 0:
                    acc_area = acc_area + grid_area   # accumulative grid area 
                    diff_area = acc_area - dam_area 
                    # print(f'acc_area = {acc_area}')
                    # print(f'diff_area = {diff_area}')
                    #-- save ix and iy
                    save_grid_ix.append(ix_search[i]+1)       #### !!!!!!!! index + 1
                    save_grid_iy.append(iy_search[i]+1)       #### !!!!!!!! index + 1
                    save_diff_area.append(abs(diff_area))
                    save_grid_area.append(grid_area)              
                    if diff_area > 0 : #and abs(diff_area) < save_diff_area[-1] :# or (diff_area < 0 and abs(diff_area) < grid_area * 1/4):
                        break
        # print(save_diff_area)
        if len(save_diff_area) > 1 : # at least two grids to be compared
            #-- compare the last two area differences：if the last value is greater than the previous value, then delete it
            if save_diff_area[-1] > save_diff_area[-2]:
                acc_area = acc_area - save_grid_area[-1]
                # diff_area = acc_area - dam_area 
                save_grid_ix = save_grid_ix[:-1]
                save_grid_iy = save_grid_iy[:-1]
                # save_diff_area = save_diff_area[:-1]

        len_ixiy = len(save_grid_ix)
        #-- save the number of grids
        # save_len_ixiy[dam_i] = len_ixiy

        ix_list = [dam_i, f'DAM_{self.GRAND_ID[dam_i]}', self.wuse_area[dam_i], acc_area, len_ixiy]
        for i in range(len(save_grid_ix)):   # write grid ix 
            ix_list.append(save_grid_ix[i])

        iy_list = [dam_i, f'DAM_{self.GRAND_ID[dam_i]}', self.wuse_area[dam_i], acc_area, len_ixiy]
        for i in range(len(save_grid_iy)):   # write grid iy
            iy_list.append(save_grid_iy[i])

        
        #-- save the number of grids
        return len_ixiy, ix_list, iy_list
    




    def p01_identify_wuse_grids(self):
        '''
        identify the wateruse_grids
        '''
        self.volToWUSE_area(self.dam_vol)

        self.read_map_data(self.mapdir)

        print('identify the wateruse_grids......')
        len_ixiy = [] # to save the number of grids of each dam

        inputlist = np.arange(self.ndam)    ####### check #######
        open(self.ix_file, 'w')
        open(self.iy_file, 'w')

        if self.ptag :
            #---
            print("Multi core operation......")
            p = Pool(self.ncores) 
            save_list = p.map(self.cal_wuse_grid, inputlist)
            p.close()

            save_len_ixiy_temp = []
            ix_list_temp = []
            iy_list_temp = []
            for temp_list in save_list:
                save_len_ixiy_temp.append(temp_list[0])
                ix_list_temp.append(temp_list[1])
                iy_list_temp.append(temp_list[2])

            max_column = max(save_len_ixiy_temp)
            first_line = '0. serial number // 1. dam id // 2. wateruse area estimated by vol~wuse_area relationship (km2) // 3. wateruse_grid area (km2) // 4. number of wateruse_grid // 5. grid-ix/iy\n'
            with open(self.ix_file, 'a+') as f:
                f.write(first_line)
                f.write(('%10s' % max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for ix in ix_list_temp:   # write grid ix 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (ix[0], ix[1], ix[2], ix[3], ix[4]))  
                    for i in range(5,len(ix)):   # write grid ix 
                        f.write('%8s' % ix[i])
                    f.write('\n')
            f.close()   
                
            with open(self.iy_file, 'a+') as f:
                f.write(first_line)
                f.write(('%10s' % max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for iy in iy_list_temp:   # write grid iy 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (iy[0], iy[1], iy[2], iy[3], iy[4]))  
                    for i in range(5,len(iy)):   # write grid iy 
                        f.write('%8s' % iy[i])
                    f.write('\n')
            f.close()

        else:
            print("Single core operation......")
            save_len_ixiy_temp = []
            ix_list_temp = []
            iy_list_temp = []
            for inpi in inputlist:
                length,ix,iy= self.cal_wuse_grid(inpi)
                save_len_ixiy_temp.append(length)
                ix_list_temp.append(ix)
                iy_list_temp.append(iy)
            
            max_column = max(save_len_ixiy_temp)
            first_line = '0. serial number // 1. dam id // 2. wateruse area estimated by vol~wuse_area relationship (km2) // 3. wateruse_grid area (km2) // 4. number of wateruse_grid // 5. grid-ix/iy\n'
            with open(self.ix_file, 'a+') as f:
                f.write(first_line)
                f.write(('%10s' % max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for ix in ix_list_temp:   # write grid ix 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (ix[0], ix[1], ix[2], ix[3], ix[4]))  
                    for i in range(5,len(ix)):   # write grid ix 
                        f.write('%8s' % ix[i])
                    f.write('\n')
            f.close()
                
            with open(self.iy_file, 'a+') as f:
                f.write(first_line)
                f.write(('%10s' % max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for iy in iy_list_temp:   # write grid iy 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (iy[0], iy[1], iy[2], iy[3], iy[4]))  
                    for i in range(5,len(iy)):   # write grid iy 
                        f.write('%8s' % iy[i])
                    f.write('\n')
            f.close()

        print('--!!! p01_calc_dam_wuse_grids finished !!!--')


       




    def read_grid_ix_iy_data(self):
        print('Reading grid_ix and grid_iy data ...')
        head_row = 2 # skip the first and second lines
        with open(self.ix_file, 'r') as f:
            lines = f.readlines()
        for line in lines[head_row:]:#[2262:2263]:#[head_row:]: #
            # print(line)
            line = line.split()
            self.dam_id.append(line[1]) # read dam_id
            self.data_area1.append(float(line[2]))
            self.data_area2.append(float(line[3]))
            self.data_col.append(int(line[4])) # read the data length of all the dams
            self.data_ix.append(np.array((line[5:])).astype(int))  # read grid_ix of all the dam 
            # print(data_ix)

        #-- read grid_iy data
        with open(self.iy_file, 'r') as f:
            lines = f.readlines()
        for line in lines[head_row:]: 
            line = line.split()
            self.data_iy.append(np.array((line[5:])).astype(int))  # read grid_iy of all the dam 
            # print(data_iy)
        self.max_column = lines[1].split()[0] # read max_column value


    def cal_n_share(self,dam_i):
            '''
            function to calculate the overlapped grids for a given dam
            '''
            # print(f'-- dam id: {self.dam_id[dam_i]}-{dam_i}')
            # print(f'dam id: {GRAND_ID[dam_i]}-{dam_i}')
            dam_ixiy = np.column_stack((self.data_ix[dam_i],self.data_iy[dam_i])) # join two one-dimensional arrays (ix and iy)
            num = self.data_col[dam_i] # data length of the dam 
            
            #=============== step[1]: get the index of dams with overlapping location =============
            #-- find dams with overlapping location
            ix_right = self.ix_lim[:,0] - self.ix_lim[dam_i,1]    # ix_right = left of other dam - right of the give dam : >0 locate at the right of the dam
            ix_left  = self.ix_lim[:,1] - self.ix_lim[dam_i,0]    # ix_left = right of other dam - left of the give dam  : <0 locate at the left of the dam
            iy_up    = self.iy_lim[:,0] - self.iy_lim[dam_i,1]    # iy_up = downward of other dam - upward of the give dam : >0 locate at the upward of the dam
            iy_down  = self.iy_lim[:,1] - self.iy_lim[dam_i,0]    # iY_down = upward of other dam - downward of the give dam : <0 locate at the downward of the dam
            #-- get the index of overlapping dams
            loc_outside = (ix_right>0) + (ix_left<0) + (iy_up>0) + (iy_down<0)  #  True to be deleted, False to be kept
            loc_idx = np.where(loc_outside == False)[0]  # index of the overlapping dams

            #=============== step[2]: search for the overlapped grids for each grid of the given dam  =============
            if len(loc_idx) == 1 : # no overlapping dam 
                save_share = [1 for i in range(num)]
            else: 
                save_share = []
                for i in range(num): # to check each grid of the given dam
                    dam_ixiy_i = dam_ixiy[i,:]

                    n_share = 0
                    for idx in range(len(loc_idx)): # compare the grid to that of the overlapping dams
                        dam_j = loc_idx[idx]
                        check_ixiy = np.column_stack((self.data_ix[dam_j],self.data_iy[dam_j])) - dam_ixiy_i 
                        exist = np.any(np.all(check_ixiy == [0, 0], axis=1)) # [0 0] means the same grid
                        if exist:
                            # index = np.where((check_ixiy == [0, 0]).all(axis=1))[0][0]
                            # print(index)
                            n_share = n_share +1  # number of [0 0]
                    save_share.append(1 / n_share) # the grid is divided equally to each dam
            # print(save_share)

            share_line = [dam_i, f'DAM_{self.GRAND_ID[dam_i]}', self.data_area1[dam_i], self.data_area2[dam_i], self.data_col[dam_i]]
            for i in range(len(save_share)):  
                share_line.append(save_share[i])
            
            return share_line


    def p02_identify_overlapping_grids(self):
        self.read_grid_ix_iy_data()

        print('Searching overlapping grids ...')
        self.ix_lim = np.zeros([self.ndam,2])
        self.iy_lim = np.zeros([self.ndam,2])
        for dam_i in range(self.ndam):
            self.ix_lim[dam_i,:] = min(self.data_ix[dam_i]), max(self.data_ix[dam_i])
            self.iy_lim[dam_i,:] = min(self.data_iy[dam_i]), max(self.data_iy[dam_i])

        #-- calculate the overlapped grids for a given dam
        inputlist = np.arange(self.ndam)
        open(self.share_file, 'w')

        if self.ptag :
            print("Multi core operation......")
            p = Pool(self.ncores)
            share_list=p.map(self.cal_n_share, inputlist)
            p.close()
            
            first_line = '0. serial number // 1. dam id // 2. wateruse area estimated by vol~wuse_area relationship (km2) // 3. wateruse_grid area (km2) // 4. number of wateruse_grid // 5. grid proportion (0, 1]\n'
            with open(self.share_file, 'a+') as f:
                f.write(first_line)
                f.write( ('%10s' % self.max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for sl in share_list:   # write grid sl 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (sl[0], sl[1], sl[2], sl[3], sl[4]))  
                    for i in range(5,len(sl)):   # write grid sl 
                        f.write('%8.4f' % sl[i])
                    f.write('\n')

        else:
            print("Single core operation......")
            share_list = []
            for inpi in inputlist:
                share_list.append(self.cal_n_share(inpi))

            first_line = '0. serial number // 1. dam id // 2. wateruse area estimated by vol~wuse_area relationship (km2) // 3. wateruse_grid area (km2) // 4. number of wateruse_grid // 5. grid proportion (0, 1]\n'
            with open(self.share_file, 'a+') as f:
                f.write(first_line)
                f.write( ('%10s' % self.max_column) + ('%32s' % 'max_number_of_wateruse_grid') + '\n')
                for sl in share_list:   # write grid sl 
                    f.write('%10s %10s %10.3f %10.3f %10s' % (sl[0], sl[1], sl[2], sl[3], sl[4]))  
                    for i in range(5,len(sl)):   # write grid sl 
                        f.write('%8.4f' % sl[i])
                    f.write('\n')
        print('--!!! p02_calc_overlapping_grids finished !!!--')




    def main_function(self):
        self.p01_identify_wuse_grids()
        self.p02_identify_overlapping_grids()


