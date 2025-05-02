import pandas as pd
import csv
import os
import numpy as np
import pandas as pd
import multiprocessing
from collections import defaultdict
import time




class dam_basicInfo_Class:
   
    def __init__(self, namelist):

        self.name       = 'dam_basicInfo_Class'
        self.version    = '0.1'
        self.release    = '0.1'
        self.date       = 'May 2023'
        self.author     = 'Shulei Zhang, zhangshlei@mail.sysu.edu.cn'   

        #-- read namelist
        self.ptag       = namelist['General'      ]['Para_Tag'  ]
        self.debug      = namelist['General'      ]['Debug_Tag' ]
        self.ncores     = int(namelist['General'  ]['Num_Cores' ])
        self.mtag       = namelist['General'      ]['Map_Tag'   ]
        self.mapdir     = namelist['General'      ]['Map_Dir'   ]
        self.savedir    = namelist['General'      ]['Save_Dir'  ]
        self.GRanD_if   = namelist['dam_basicInfo']['GRanD_If'  ]
        self.minerror   = float(namelist['dam_basicInfo']['Min_Error'  ])
        self.minuparea  = float(namelist['dam_basicInfo']['Min_Uparea'  ])
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)

        # save file
        self.GRanD_of   = self.savedir + '/dam_inplist.csv'
        self.damfile    = self.savedir + '/damloc.csv'


    def p01_creat_damlist(self):
        # read original file
        GRanD_data = pd.read_csv(self.GRanD_if)
        GRanD_column = GRanD_data.columns.tolist()

        # extract selected data
        column_names = ['GRAND_ID','DAM_NAME','LONG_DD','LAT_DD','CAP_MCM','CATCH_SKM','MAIN_USE','YEAR']
        dam_inp = GRanD_data[column_names]

        # delete dam with area = 0
        dam_inp = dam_inp[dam_inp['CATCH_SKM'] > 0]

        # update dam_name
        dam_inp['DAM_NAME'] = dam_inp['DAM_NAME'].str.replace(" ", "-")
        nan_rows = dam_inp['DAM_NAME'].isna()
        dam_inp.loc[nan_rows, 'DAM_NAME'] = "GRanD-" + dam_inp.loc[nan_rows, 'GRAND_ID'].astype(str)
        if self.debug :
            print("dam_inp.loc[nan_rows, 'DAM_NAME']",dam_inp.loc[nan_rows, 'DAM_NAME'])

        # update main_use
        dam_inp['MAIN_USE'] = dam_inp['MAIN_USE'].fillna('Other-NA')
        dam_inp['MAIN_USE'] = dam_inp['MAIN_USE'].str.replace(" ", "-")
        if self.debug :
            print("dam_inp['MAIN_USE'].unique()",dam_inp['MAIN_USE'].unique())

        # save file
        dam_inp.to_csv(self.GRanD_of, index=False)

        with open(self.GRanD_of, 'r') as f:
            reader = csv.reader(f)
            data = [row for row in reader]
        new_row = [f"{len(data) - 1}", "NDAMS"]
        data.insert(0, new_row)
        with open(self.GRanD_of, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(data)




    def p02_identify_damloc(self):
        #  Code to allocate GRanD reservoir on CaMa-Flood river map
        #  -- input:  GRanD reservoir list      
        #  -- output: tentative allocation file 
        #
        #  written by Risa Hanazaki, 2021Mar (get_rivinfo_glb.F90)
        #  converted to python and updated by Shulei Zhang (zhangshlei@mail.sysu.edu.cn), May 2023

        self.outdir = self.savedir
        self.damtmpfile = os.path.join(self.outdir, 'tmp_damloc.csv')
        if os.path.exists(self.damtmpfile):
            os.remove(self.damtmpfile)
        self.check_dir(self.outdir)
        
        out_vars = ['GRAND_ID','ix', 'iy', 'uparea_cama']

        # read dam info
        self.dam_info = pd.read_csv(self.GRanD_of, header=1)
        ndams    = self.dam_info.shape[0]
        if self.debug :
            print(self.dam_info)
            print("---------------------------")
            print("")
        
        # read map files
        fparam = os.path.join(self.mapdir, 'params.txt')
        with open(fparam, 'r') as f:
            lines = f.readlines()
            self.nx = int(lines[0].split()[0])
            self.ny = int(lines[1].split()[0])
            self.gsize = float(lines[3].split()[0])
            self.west = float(lines[4].split()[0])
            self.east = float(lines[5].split()[0])
            self.south = float(lines[6].split()[0])
            self.north = float(lines[7].split()[0])
        f.close()
        if self.debug :
            print('Map Domain W-E-S-N: ', self.west, self.east, self.south, self.north)
            print('Map Resolution    : ', self.gsize)
            print('Map NX,NY         : ', self.nx, self.ny)
            print("---------------------------")
            print("")

        # read map bin files
        self.read_bin_data()
        save_list = []
        # calculate ix iy for each dam
        if self.ptag:
            print("parallel computing..................................")
            pool = multiprocessing.Pool(self.ncores)
            save_list = pool.map(self.process_dam, range(ndams))

        else:
            print("serial computing....................................")
            for inp in range(ndams):
                save_temp = self.process_dam(inp)
                save_list.append(save_temp)
        # print("save_list",save_list)

        # save data
        out_data = pd.DataFrame(save_list, columns = out_vars)
        out_data = out_data.sort_values(by=out_data.columns[0])
        out_data = pd.merge(self.dam_info,out_data,on='GRAND_ID')
        out_data.to_csv(self.damtmpfile, index=False)




    def p03_complete_damcsv(self):
        # save file
        self.damfile = os.path.join(self.outdir, 'damloc.csv')
        
        self.damcsv  = pd.read_csv(self.damtmpfile)

        # treat multiple dams in one grid
        if self.debug :
            print('')
            print('treat multiple dams in one grid')
        # count grids with multiple dams allocated
        cnt = defaultdict(int)
        for index, row in self.damcsv.iterrows():
            grandid, ix, iy = row['GRAND_ID'], row['ix'], row['iy']
            key = tuple([ix,iy])
            cnt[key] += 1
        if self.debug :
            print('cnt=', cnt)

        # remove smaller dams in such grids
        damcsv_update = self.damcsv.copy()
        for k, v in cnt.items():
            if v > 1:
                if self.debug :
                    print('')
                    print('multiple dams on one grid!!:', k, v)
                ix, iy = k
                dams = self.damcsv.query('ix == @ix & iy == @iy')
                maxsto = np.max(dams['CAP_MCM'])
                rmdams = dams.query('CAP_MCM != @maxsto')
                # if len(rmdams) == 0:
                #     maxfsto = np.max(dams['fldsto_mcm'])
                #     rmdams = dams.query('fldsto_mcm != @maxfsto')
                if self.debug :
                    print('remove:', rmdams)
                damcsv_update.drop(index=rmdams.index, inplace=True)
                if self.debug :
                    print(damcsv_update.query('ix==@ix & iy==@iy'))
        if self.debug :
            print('-----------------------------------')

        # remove dams with small drainage area
        if self.debug :
            print('')
            print('remove dams with small drainage area')
            print('-----------------------------------')
        damcsv_update = damcsv_update.query('uparea_cama >= @self.minuparea')
        damcsv_update = damcsv_update.dropna()


        # save output
        damcsv_update.to_csv(self.damfile, index=None)
        if self.debug :
            print(' ')
            print('dam locations:', self.damfile)

        # add the first row
        with open(self.damfile, 'r') as f:
            reader = csv.reader(f)
            data = [row for row in reader]
        new_row = [f"{len(data) - 1}", "NDAMS"]
        data.insert(0, new_row)
        with open(self.damfile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(data)




    def read_bin_data(self):

        finp = os.path.join(self.mapdir, 'nextxy.bin')
        with open(finp, 'rb') as f:
            self.nextx = np.fromfile(f, dtype=np.int32, count=self.nx * self.ny, sep='')
            self.nexty = np.fromfile(f, dtype=np.int32, count=self.nx * self.ny, sep='')
        self.nextx = np.reshape(self.nextx, (self.nx, self.ny), order='F')
        self.nexty = np.reshape(self.nexty, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read nextxy.bin successfully')

        finp = os.path.join(self.mapdir, 'basin.bin')
        with open(finp, 'rb') as f:
            self.basin = np.fromfile(f, dtype=np.int32, count=self.nx * self.ny, sep='')
        self.basin = np.reshape(self.basin, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read basin.bin successfully')

        finp = os.path.join(self.mapdir, 'upgrid.bin')
        with open(finp, 'rb') as f:
            self.upgrid = np.fromfile(f, dtype=np.int32, count=self.nx * self.ny, sep='')
        self.upgrid = np.reshape(self.upgrid, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read upgrid.bin successfully')

        finp = os.path.join(self.mapdir, 'uparea.bin')
        with open(finp, 'rb') as f:
            self.uparea = np.fromfile(f, dtype=np.float32, count=self.nx * self.ny, sep='')
        self.uparea = np.reshape(self.uparea, (self.nx, self.ny), order='F')         
        self.uparea = self.uparea * 1e-6      ##### unit conversion !!!!!!!!!
        if self.debug :
            print('Read uparea.bin successfully')

        finp = os.path.join(self.mapdir, 'ctmare.bin')
        with open(finp, 'rb') as f:
            self.ctmare = np.fromfile(f, dtype=np.float32, count=self.nx * self.ny, sep='')
        self.ctmare = np.reshape(self.ctmare, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read ctmare.bin successfully')

        finp = os.path.join(self.mapdir, 'elevtn.bin')
        with open(finp, 'rb') as f:
            self.elevtn = np.fromfile(f, dtype=np.float32, count=self.nx * self.ny, sep='')
        self.elevtn = np.reshape(self.elevtn, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read elevtn.bin successfully')

        finp = os.path.join(self.mapdir, 'lonlat.bin')
        with open(finp, 'rb') as f:
            self.outlon = np.fromfile(f, dtype=np.float32, count=self.nx * self.ny, sep='')
            self.outlat = np.fromfile(f, dtype=np.float32, count=self.nx * self.ny, sep='')
        self.outlon = np.reshape(self.outlon, (self.nx, self.ny), order='F')
        self.outlat = np.reshape(self.outlat, (self.nx, self.ny), order='F')
        if self.debug :
            print('Read lonlat.bin successfully')




    def process_dam(self, dam):
        grandid = self.dam_info['GRAND_ID'][dam]
        damname = self.dam_info['DAM_NAME'][dam]
        lon = np.float32(self.dam_info['LONG_DD'][dam])
        lat = np.float32(self.dam_info['LAT_DD'][dam])
        totalsto = np.float32(self.dam_info['CAP_MCM'][dam])
        upreal   = np.float32(self.dam_info['CATCH_SKM'][dam])
        ix = 0
        iy = 0
        if self.debug :
            print('')
            print("grandid:", grandid, " damname:", damname, " uparea:", upreal," lon:", lon, " lat:", lat, " totalsto:", totalsto)
            print("nx:", self.nx, " ny:", self.ny,"west:", self.west, " north:", self.north,  ix, iy)
        ix,iy = self.calc_ixiy(lon, lat, ix, iy)
        if self.debug : 
            print("ix,iy", ix, iy)

        if ix < 0 or iy < 0:
            for cnt in range(1, 5):
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        if dx == 0 and dy == 0:
                            continue
                        lon = lon + self.gsize * cnt * dy
                        lat = lat + self.gsize * cnt * dx
                        ix, iy = self.calc_ixiy(lon, lat, ix, iy)
                        if ix > 0 and iy > 0:
                            break
                    else:
                        continue
                    break

        if (ix < 0 or iy < 0):
            print(grandid, damname, lon, lat, ix, iy, upreal, "undefined", totalsto)

            #-- save the identified ix iy
            # return grandid, ix, iy, -9999

        elif (ix > 0 and iy > 0):
            print(grandid,damname, lon, lat, ix, iy, upreal, self.uparea[ix-1,iy-1], totalsto)

            #-- check area error
            error = abs(self.uparea[ix-1,iy-1] - upreal)
            if error > self.minerror * upreal:
                ix, iy, error = self.modify_damloc(ix, iy, error, upreal, self.uparea)

            #-- save the identified ix iy
            return grandid, ix, iy, self.uparea[ix-1,iy-1]


    def calc_ixiy(self, lon, lat, ix, iy):
        isHires = 0
        # print('')

        # check 15sec hires map availability
        floc = self.mapdir + '/15sec/location.txt'
        if os.path.isfile(floc):
            with open(floc, 'r') as f:
                f.readline()
                f.readline()
                line = f.readline()
                buf, tag, _, _, _, _, mx, my, csize = line.split()
                mx=int(mx)
                my=int(my)
                csize=np.float32(csize)
            if tag == '15sec':
            #print('USE 15sec hires map')
                isHires = 15

                catmx = np.zeros((mx, my), dtype=np.int32)
                catmy = np.zeros((mx, my), dtype=np.int32)
                finp = self.mapdir + '/15sec/15sec.catmxy.bin'
                print('finp', finp)
                with open(finp, 'rb') as f:
                    catmx_bytes = f.read(2*mx*my)
                    catmx = np.frombuffer(catmx_bytes, dtype=np.int16).astype(np.int32).reshape((mx, my))
                    catmy_bytes = f.read(2*mx*my)
                    catmy = np.frombuffer(catmy_bytes, dtype=np.int16).astype(np.int32).reshape((mx, my))

                jx = int((lon-self.west)/csize) + 1
                jy = int((self.north-lat)/csize) + 1
                if jx <= 0 or jx > mx or jy <= 0 or jy > my:
                    # print('ix,iy cannot be defined')
                    ix = -99
                    iy = -99
                    # stop
                ix = catmx[jx-1, jy-1]
                iy = catmy[jx-1, jy-1]
                if ix <= 0 or ix > self.nx or iy <= 0 or iy > self.ny:
                    # print('ix,iy cannot be defined')
                    ix = -99
                    iy = -99
                    # stop

        if isHires != 15:
            # check 1min hires map availability
            floc = self.mapdir + '/1min/location.txt'
            if os.path.isfile(floc):
                with open(floc, 'r') as f:
                    f.readline()
                    f.readline()
                    line = f.readline()
                    buf, tag, _, _, _, _, mx, my, csize = line.split()
                    mx=int(mx)
                    my=int(my)
                    csize=np.float32(csize)
                if tag == '1min':
                    isHires = 60
                    finp = self.mapdir + '/1min/1min.catmxy.bin'
                    with open(finp, 'rb') as f:
                        catmx = np.fromfile(f, dtype=np.int16, count=mx*my).reshape((mx, my), order='F')
                        catmy = np.fromfile(f, dtype=np.int16, count=mx*my).reshape((mx, my), order='F')
                    jx = int((lon - self.west) / csize) + 1
                    jy = int((self.north - lat) / csize) + 1
                    if jx <= 0 or jx > mx or jy <= 0 or jy > my:
                        # print('ix,iy cannot be defined')
                        ix = -99
                        iy = -99
                        # raise ValueError('ix,iy cannot be defined')
                    ix = catmx[jx-1, jy-1]
                    iy = catmy[jx-1, jy-1]
                    if ix <= 0 or ix > self.nx or ix <= 0 or iy > self.ny:
                        # print('ix,iy cannot be defined')
                        ix = -99
                        iy = -99
                        # raise ValueError('ix,iy cannot be defined')
        if isHires == 0:
            ix = int((lon-self.west) / self.gsize) + 1
            iy = int((self.north-lat) / self.gsize) + 1

            glon = self.west + self.gsize * (ix - 0.5)
            glat = self.north - self.gsize * (iy - 0.5)

        if ix > 0 and iy > 0:
        # print("nextx(ix,iy):", nextx[ix-1, iy-1])
            if self.nextx[ix-1, iy-1] == -9999:
                print('NOT LAND GRID')
                ix = -99
                iy = -99
        # print(ix,iy)
        return ix,iy




    


    def modify_damloc(self, ix, iy, error, upreal, uparea):
        if self.debug :
            print("error >= uparea_real*minerror ; modify dam location")
            print("uparea_real=", upreal, " uparea=", uparea[ix-1,iy-1], " error=", error)  

        # searching -------------------------------------------
        ix_m,iy_m = ix,iy
        error_m   = error
        for j_x in range(-1, 2, 1):
            for j_y in range(-1, 2, 1):
                ix_tmp = ix + j_x
                iy_tmp = iy + j_y
                error_tmp = abs(uparea[ix_tmp-1,iy_tmp-1] - upreal)

                if error_tmp >= error_m:    #if there is no better grid????
                    continue
                else: #error_tmp<error_m
                    ix_m, iy_m = ix_tmp, iy_tmp
                    error_m = error_tmp
                if self.debug :
                    print("modified location:", ix_m, iy_m, 'up_real=', uparea[ix_m-1,iy_m-1], 'error=', error_m)

        if error_m >= self.minerror * upreal:
            if self.debug :
                print('')
                print("still have error >= minerror!!!!!!")
            ix_m,iy_m = ix,iy
            error_m   = error

            for j_x in range(-2, 3, 1):
                for j_y in range(-2, 3, 1):
                    ix_tmp = ix + j_x
                    iy_tmp = iy + j_y
                    error_tmp = abs(uparea[ix_tmp-1,iy_tmp-1] - upreal)
    
                    if error_tmp >= error_m:    #if there is no better grid????
                        continue
                    else: #error_tmp<error_m
                        ix_m, iy_m = ix_tmp, iy_tmp
                        error_m = error_tmp
                        if self.debug :
                            print("modified location:", ix_m, iy_m, 'up_cama=', uparea[ix_m-1,iy_m-1], 'error=', error_m)
        
        if self.debug :
            print("final modified location:", ix_m, iy_m, 'up_cama=' ,uparea[ix_m-1,iy_m-1], 'error=', error_m)

        return ix_m,iy_m,error_m




    def check_dir(self,dir):
        if os.path.exists(dir):
            if self.debug :
                print (dir+' exists')
        else:
            os.makedirs(dir)
            if self.debug :
                print (dir+' created')



    def main_func(self):
        start_time = time.time()

        print('--- Start Process Dam Basic Info ---')
        self.p01_creat_damlist()
        self.p02_identify_damloc()
        self.p03_complete_damcsv()
        print('--- Done ---')





