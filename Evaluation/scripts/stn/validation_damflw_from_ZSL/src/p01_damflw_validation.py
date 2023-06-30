import calendar
import datetime
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
def filter_nan(s, o):
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
def NS(s, o):
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
    s, o = filter_nan(s, o)

    ns = 1 - sum((s - o) ** 2) / (sum((o - np.mean(o)) ** 2) + 1e-20)

    return ns

# ========================================
def NSlog(s, o):
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
    s, o = filter_nan(s, o)

    s = np.log(s)      # warning!!!
    o = np.log(o)

    nslog = 1 - sum((s - o) ** 2) / (sum((o - np.mean(o)) ** 2) + 1e-20)
    return nslog

# ========================================
def KGE(s, o):
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
    s, o = filter_nan(s, o)

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
def read_sim_data(inputlist):    # input: year of data
    yyyy = inputlist

    sim_sto = np.ctypeslib.as_array(shared_array_sim_sto)
    sim_inflw = np.ctypeslib.as_array(shared_array_sim_inflw)
    sim_outflw = np.ctypeslib.as_array(shared_array_sim_outflw)

    s_file = s_dir + '/damtxt-' + str(yyyy) + '.txt'
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
        st = (target_dt - start_dt).days
        et = st + dt
     
        # Find common id
        common_id = set(dam_id).intersection(dam_id_out)
        # Find indices of common numbers in dam_id_out
        idx_out = [i for i, x in enumerate(dam_id_out) if x in common_id]
        # index of data ######## revise according to .txt format
        id_sto=  [x * 3 + 1 for x in idx_out]
        id_inflw=  [x * 3 + 2 for x in idx_out]
        id_outflw=  [x * 3 + 3 for x in idx_out]
        # Find indices of common numbers in dam_id
        idx_id = [i for i, x in enumerate(dam_id) if x in common_id]

        # data of yyyy
        sto_temp = dam_data[:, id_sto]
        inflw_temp = dam_data[:, id_inflw]
        outflw_temp = dam_data[:, id_outflw]
        # all the data 
        sim_sto[st:et,idx_id] = sto_temp
        sim_inflw[st:et,idx_id] = inflw_temp
        sim_outflw[st:et,idx_id] = outflw_temp

        # print(sto_temp)

        # # Read storage, inflow, and outflow data
        # for i in range(len(idx_out)):
        #     # print(dam_id_out[idx_out[i]])
        #     id_data = idx_out[i] * 3 +1
        #     # data at idx_out[i] from damtxt-xxxx.txt
        #     sto_temp = [float(line.split()[id_data]) for line in lines[head::]]
        #     # inflw_temp = [float(line.split()[id_data+1]) for line in lines[head::]]
        #     # outflw_temp = [float(line.split()[id_data+2]) for line in lines[head::]]
        #     # write data into idx_id[i] 
        #     # print(dam_id[idx_id[i]])
        #     sim_sto[st:et,idx_id[i]] = sto_temp
        #     # sim_inflw[st:et,idx_id[i]] = inflw_temp
        #     # sim_outflw[st:et,idx_id[i]] = outflw_temp

        # sto_temp = []
        # inflw_temp = []
        # outflw_temp = []
        # for i in range(ndam):
        #     dam_id_str = str(dam_id[i])
        #     # print(dam_id_str)
        #     # check if the output exists
        #     if dam_id_str in dam_id_out: 
        #         idx_id = dam_id_out.index(dam_id_str)
        #         idx_id = idx_id * 3 + 1    ##### revise according to output format
        #         # read storage
        #         sto_temp = [float(line.split()[idx_id]) for line in lines[head::]]
        #         # print(sto_temp)
        #         # read inflow
        #         inflw_temp = [float(line.split()[idx_id+1]) for line in lines[head::]]
        #         # print(inflw_temp)
        #         # read outflw
        #         outflw_temp = [float(line.split()[idx_id+2]) for line in lines[head::]]
        #         # print(outflw_temp)

        #         # each year
        #         sim_sto[st:et,i] = sto_temp
        #         # inflw[st:et,i] = inflw_temp
        #         # outflw[st:et,i] = outflw_temp

#==================================== 
## read observation
def read_obs_data(inputlist):    ## input dam id
    o_file = o_dir + '/ResOpsUS_' + str(dam_id[inputlist]) + '.csv'

    obs_sto = np.ctypeslib.as_array(shared_array_obs_sto)
    obs_inflw = np.ctypeslib.as_array(shared_array_obs_inflw)
    obs_outflw = np.ctypeslib.as_array(shared_array_obs_outflw)

    # read observation 
    head = 1  # header lines      # change according to the obs.txt format -- revise by shulei
    if not os.path.exists(o_file):
        # print("no file", o_file)
        sto = np.ones([last], np.float32) * -9999.0
        inflw = np.ones([last], np.float32) * -9999.0
        outflw = np.ones([last], np.float32) * -9999.0
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
                if start_dt <= data_dt <= last_dt:
                    sto_temp[yyyy + mm + dd] = float(line_one[1])       ## 2column- storage
                    inflw_temp[yyyy + mm + dd] = float(line_one[2])     ## 3column- inflw
                    outflw_temp[yyyy + mm + dd] = float(line_one[3])    ## 4column- outflw
                elif last_dt < data_dt:
                    break            
        # ----- check data date with input date
        sto = []
        inflw = []
        outflw = []
        for day in np.arange(last):
            target_dt = start_dt + datetime.timedelta(days=int(day))
            yyyy = '%04d' % (target_dt.year)
            mm = '%02d' % (target_dt.month)
            dd = '%02d' % (target_dt.day)
            if (yyyy + mm + dd) in sto_temp.keys():
                sto.append(sto_temp[yyyy + mm + dd])
                inflw.append(inflw_temp[yyyy + mm + dd])
                outflw.append(outflw_temp[yyyy + mm + dd])
            else:
                sto.append(-9999.0)
                inflw.append(-9999.0)
                outflw.append(-9999.0)

        obs_sto[:,inputlist] = sto
        obs_inflw[:,inputlist] = inflw
        obs_outflw[:,inputlist] = outflw

# =====================================
# def write_result(dam_i, ):

    # s = sim_sto
    # [s,o] = filter_nan(sim, obs)

    # NSval = NS(sim, obs)
    # # NSLval = NSlog(sim, obs)
    # KGEval = KGE(sim, obs)


    # fname = res_file
    # with open(fname, "a") as f:
    #     vali = "%12s %12.4f %12.4f %12.2f %12.2f %12.4f %12.4f" \
    #         % (pnames[dam_i], latlist[dam_i], lonlist[dam_i], NSval, KGEval)
    #     f.write(vali + '\n')
    # f.close()
    # return 0

#=======================================    
def compare_s_and_o(dam_i):
    min_len =  50 # Minimum requirements for validation data [days]

    cval = np.ctypeslib.as_array(shared_array_cval)    
    #---------
    compare_val = np.zeros((1,6)) * np.nan      
    #---------
    s = sim_sto[:,dam_i]
    o = obs_sto[:,dam_i]
    [s,o] = filter_nan(s, o)
    if len(s) > min_len:
        compare_val[0,0] = NS(s, o)
        compare_val[0,1] = KGE(s, o)
    #---------
    s = sim_inflw[:,dam_i]
    o = obs_inflw[:,dam_i]
    [s,o] = filter_nan(s, o)
    if len(s) > min_len:
        compare_val[0,2] = NS(s, o)
        compare_val[0,3] = KGE(s, o)
    #---------
    s = sim_outflw[:,dam_i]
    o = obs_outflw[:,dam_i]
    [s,o] = filter_nan(s, o)
    if len(s) > min_len:
        compare_val[0,4] = NS(s, o)
        compare_val[0,5] = KGE(s, o)
    #---------
    cval[dam_i,:] = compare_val


start_time = time.time()
#====================================  setting  ====================================
# syear, smon, sday = 1979, 1, 1
# eyear, emon, eday = 1980, 12, 31
# tag = 'src_v411_15min_h06_test'

# ## read dam list
# dam_dir = '/tera02/zhangsl/cama-flood/cama_v4.10/map/dam_para'
# dam_file = dam_dir + '/dam_nw_06min/dam_params_nw_06min.csv'

# ## read observation
# o_dir = '/tera02/zhangsl/cama-flood/obs/ResOpsUS/time_series_all'

# ## read simluation
# s_dir = '/tera02/zhangsl/cama-flood/cama_v4.10/out/' + tag #src_nw_06min_dam_h06_test'

# ## save results
# res_dir = '/tera02/zhangsl/cama-flood/cama_v4.10/etc/validation_damflw/validation/' + tag 
# if not os.path.exists(res_dir):
#     os.mkdir(res_dir)
# res_file = res_dir + '/validation_' + str(syear) + '-' + str(eyear) + '.csv'

# ## parallel calculation set
# para_flag = 1
# num_cores = 72

syear, smon, sday = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
eyear, emon, eday = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
tag = sys.argv[7]

## read dam list
dam_file = f'./inp/damlist.csv'

## read observation
o_dir = f'./inp/obs'

## read simluation
s_dir = f'./inp/sim'

## save results
res_file = f'./validation_{syear}-{eyear}.csv'

## parallel calculation set
para_flag = int(sys.argv[8])
num_cores = int(sys.argv[9])


#====================================  setting  ====================================
dam_data = pd.read_csv(dam_file,skiprows=[0])
dam_data = dam_data.sort_values('GRAND_ID')

# dam_data = dam_data.iloc[360:361]
dam_id = np.array(dam_data['GRAND_ID'])
dam_name = np.array(dam_data['DamName'])
dam_lon = np.array(dam_data['DamLon'])
dam_lat = np.array(dam_data['DamLat'])

ndam = len(dam_id)
# print(np.where(dam_id==1207)[0])
# print(dam_id)

start_dt = datetime.date(syear, smon, sday)
last_dt = datetime.date(eyear, emon, eday)
last = (last_dt - start_dt).days + 1

#====================================  step[1] read simulaiton  ====================================
# multiprocessing array
sim_sto = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_sim_sto = sharedctypes.RawArray(sim_sto._type_, sim_sto)
sim_inflw = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_sim_inflw = sharedctypes.RawArray(sim_inflw._type_, sim_inflw)
sim_outflw = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_sim_outflw = sharedctypes.RawArray(sim_outflw._type_, sim_outflw)

# for parallel calculation
inputlist = np.arange(syear, eyear + 1)
# print(inputlist)

# read sim_data parallel
if para_flag == 1:
    if __name__ == "__main__":            
        p = Pool(num_cores)
        res = list(p.map(read_sim_data, inputlist))
        sim_sto = np.ctypeslib.as_array(shared_array_sim_sto)
        sim_inflw = np.ctypeslib.as_array(shared_array_sim_inflw)
        sim_outflw = np.ctypeslib.as_array(shared_array_sim_outflw)
        p.terminate()
else:
    for inpi in inputlist:
        res = read_sim_data(inpi)
        sim_sto = np.ctypeslib.as_array(shared_array_sim_sto)
        sim_inflw = np.ctypeslib.as_array(shared_array_sim_inflw)
        sim_outflw = np.ctypeslib.as_array(shared_array_sim_outflw)

# print(sim_sto)
# print(sim_inflw)
# print(sim_outflw)

#====================================  step[2] read observation  ====================================
# multiprocessing array
obs_sto = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_obs_sto = sharedctypes.RawArray(obs_sto._type_, obs_sto)
obs_inflw = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_obs_inflw = sharedctypes.RawArray(obs_inflw._type_, obs_inflw)
obs_outflw = np.ctypeslib.as_ctypes(np.full((last, ndam), -9999.0, dtype=np.float32))
shared_array_obs_outflw = sharedctypes.RawArray(obs_outflw._type_, obs_outflw)

# for parallel calculation
inputlist = np.arange(ndam)
# print(inputlist)

# read observation
if para_flag == 1:
    if __name__ == "__main__":            
        p = Pool(num_cores)
        res = list(p.map(read_obs_data, inputlist))
        obs_sto = np.ctypeslib.as_array(shared_array_obs_sto)
        obs_inflw = np.ctypeslib.as_array(shared_array_obs_inflw)
        obs_outflw = np.ctypeslib.as_array(shared_array_obs_outflw)
        p.terminate()
else:
    for inpi in inputlist:
        res = read_obs_data(inpi)
        obs_sto = np.ctypeslib.as_array(shared_array_obs_sto)
        obs_inflw = np.ctypeslib.as_array(shared_array_obs_inflw)
        obs_outflw = np.ctypeslib.as_array(shared_array_obs_outflw)

# print(obs_sto)
# print(obs_inflw)
# print(obs_outflw)

#====================================  step[3] compare sim vs obs  ====================================
# multiprocessing array
cval = np.ctypeslib.as_ctypes(np.full((ndam, 6), np.nan, dtype=np.float32))
shared_array_cval = sharedctypes.RawArray(cval._type_, cval)

# for parallel calculation
inputlist = np.arange(ndam)
# print(inputlist)

# compare
if para_flag == 1:
    if __name__ == "__main__":            
        p = Pool(num_cores)
        res = list(p.map(compare_s_and_o, inputlist))
        cval = np.ctypeslib.as_array(shared_array_cval)
        p.terminate()
else:
    for inpi in inputlist:
        res = compare_s_and_o(inpi)
        cval = np.ctypeslib.as_array(shared_array_cval)

# save validation for all dam
df = pd.DataFrame({'GRAND_ID':dam_id, 'DamName':dam_name, 'DamLon':dam_lon, 'DamLat':dam_lat,
                   'NSval_sto': cval[:,0], 'KGEval_sto': cval[:,1], 
                   'NSval_inflw': cval[:,2], 'KGEval_inflw': cval[:,3],
                   'NSval_outflw': cval[:,4], 'KGEval_outflw': cval[:,5]})

# delete dam with no data
df.dropna(subset=df.columns[-6:], how='all', inplace=True)
print(df)
df.to_csv(res_file, na_rep='NA',index=False)


# end timer
end_time = time.time()

# calculate and print elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")