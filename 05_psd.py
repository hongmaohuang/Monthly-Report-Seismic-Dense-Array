# ===================================
# 在project資料夾下mkdir output_psd
# ===================================
# %%  
import obspy
import glob 
import pandas as pd
from obspy import UTCDateTime
from obspy.io.sac import attach_paz
from obspy.signal.invsim import corn_freq_2_paz
from matplotlib.dates import date2num
from matplotlib import pyplot as plt
import os
import shutil
import numpy as np
from obspy.signal import PPSD
from obspy.io.xseed import Parser

startdate = '20230103'
project_path = '/home/hmhuang/Work/Research_Assistant/Ilan_geothermal/'
sac_path = '/raid2/ILAN2022/Sorted_oneday/' 
resp_path = '/raid2/ILAN2022/PZs/'
stations_file = '/raid2/ILAN2022/stations.csv'
output_path = project_path + 'output_psd/'
sensitivity_resp = 306846
corner_freq = 5
damping = 0.7
gain = 1
component = ['DPZ']
pre_filt = (0.05, 0.1, 30.0, 50.0)

all_date = sorted(os.listdir(sac_path))
del all_date[0:all_date.index(startdate)]

sta = sorted(pd.read_csv(stations_file).Station)
st = obspy.Stream()

for k in sta:
    if os.path.isfile(output_path + str(k) + '_psd.png'):
        print('file: '  + str(k) + '_psd.png is exist !')
        continue
    else:
        print('file: ' + str(k) + '_psd.png is processing !')
        all_path = []
        for i in range(len(all_date)):
            all_path.insert(i, glob.glob(f'{sac_path}{all_date[i]}/*{k}*{component[0]}*.SAC'))
        for j in range(len(all_path)):
            try:
                st += obspy.read(all_path[j][0], dtype='float16')
                st.merge(method=1, fill_value='interpolate')                              
            except:
                continue
        if len(st) == 0:
            print('No Data in Station: ' + str(k))
            continue
        else:
            tr = st[0]
            attach_paz(tr, resp_path + component[0])
            tr.detrend(type='demean')  
            tr.taper(max_percentage=0.05, type='cosine', max_length=len(tr[0].data), side='both')
            paz = {'gain': gain,
            'poles': tr.stats.paz.poles,
            'sensitivity': sensitivity_resp,
            'zeros': tr.stats.paz.zeros}
            ppsd = PPSD(tr.stats, paz)
            ppsd.add(st)
            ppsd.plot(filename = output_path + str(k) + '_psd.png', cmap=obspy.imaging.cm.pqlx, period_lim=(0.5,100), xaxis_frequency=True)
