# ================================
# 畫15天的時頻圖（所有站點）
# 請先在project資料夾下：
# mkdir -r output_15days
# ================================
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


startdate = '20221101'
part_plot = '1st'
# 1st or 2nd >>上半月or下半月
# 202210月更新 >> 可直接畫30天ㄉ
project_path = '/home/hmhuang/Work/Research_Assistant/Ilan_geothermal/'
sac_path = '/raid2/ILAN2022/Sorted_oneday/' 
resp_path = '/raid2/ILAN2022/PZs/'
stations_file = '/raid2/ILAN2022/stations.csv'
output_path = project_path + 'output_Spec/'
# ======================
# 可改output資料夾
#output_path = project_path + 'output_active/'
# ======================
sensitivity_resp = 306846
corner_freq = 5
damping = 0.7
gain = 1
component = ['DPZ']

paz_1hz = corn_freq_2_paz(corner_freq, damp=damping)  
paz_1hz['sensitivity'] = sensitivity_resp
paz_1hz['gain'] = gain
pre_filt = (0.05, 0.1, 110.0, 125.0)


all_date = sorted(os.listdir(sac_path))
del all_date[0:all_date.index(startdate)]
sta = sorted(pd.read_csv(stations_file).Station)
# ======================
# 要改特定站點從這邊改
#del sta[3:len(sta)]
# ======================
st = obspy.Stream()
# %%
for k in sta:
    if os.path.isfile(output_path + str(k) + ".png"):
        print('file: ' + str(k) + '.png is exist !')
        continue
    else:
        all_path = []
        for i in range(len(all_date)):
            all_path.insert(i, glob.glob(f'{sac_path}{all_date[i]}/*{k}*{component[0]}*.SAC'))
        if part_plot == '1st':
            start_day = 0
            end_day = int(len(all_path)/2)
        else:
            start_day = int(len(all_path)/2) +1
            end_day = len(all_path)
        start_day = 0
        end_day = len(all_path)
        # ======================
        # 要改特定日期從這邊改
        #start_day = 1
        #end_day = start_day + 10
        # ======================
        for j in range(start_day, end_day):
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
            tr.decimate(3)
            attach_paz(tr, resp_path + component[0])
            #st.simulate(paz_remove=tr.stats.paz, pre_filt=pre_filt)
            st.simulate(paz_remove=paz_1hz, pre_filt=pre_filt)
            tr = st[0]
            tr.detrend(type='demean')  
            tr.taper(max_percentage=0.05, type='cosine', max_length=len(tr[0].data), side='both')
            # =======
            # 時頻圖
            # =======
            sig = tr.data
            time = tr.times("matplotlib")
            labelUTCtime1 = tr.stats.starttime
            time1 = date2num(labelUTCtime1.datetime)
            labelUTCtime2 = tr.stats.endtime
            time2 = date2num(labelUTCtime2.datetime)
            NFFT = 1024*8
            num_overlap = 512*8
            Fs = tr.stats.sampling_rate
            hspace = 0.7
            fig = plt.figure(figsize=(10,7))
            ax1 = plt.subplot2grid((3,1), (0,0), rowspan=1)
            ax1.plot(time, sig, 'k', linewidth=1)
            ax1.set_xlim(time1, time2)
            ax1.xaxis_date()
            fig.autofmt_xdate()
            ax1.grid(True)
            plt.title(tr.stats.station + '.' + tr.stats.channel)
            plt.subplots_adjust(hspace=hspace)  
            plt.jet()
            ax2 = plt.subplot2grid((3,1), (1,0), rowspan=2)
            ax2.specgram(sig, Fs=Fs, NFFT=NFFT, noverlap=num_overlap,vmin=-250,vmax=-100, scale='dB')
            ax2.set_yscale('symlog')
            cbar_x = ax2.get_position().x1 + 0.01 
            cbar_y = ax2.get_position().y0
            cbar_h = ax2.get_position().height
            cbar = fig.add_axes([cbar_x, cbar_y, 0.02, cbar_h])
            ax2.set_ylim([0, 50])
            ax2.set_xlabel("Time [sec]")
            ax2.set_ylabel("Frequency [Hz]")
            plt.colorbar(ax2.images[0], label='(dB)', cax=cbar) 
            plt.rcParams['font.sans-serif'] =  'Nimbus Roman'
            fig.savefig(output_path + tr.stats.station + ".png")
            print(tr.stats.station + ".png is finished")
            plt.close()
