# ================================
# 畫一個指定地震的時頻圖（所有站點）
# 1st Step for the Routine!!
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
# ========
# 可改參數
# ========

eq_date = '20230113'
eq_time = '185601'
lon_epic = 122.01
lat_epic = 24.92
th_line = 0.003
# 選站門檻值

# ==========
# 比較不用改
# ==========
project_path = '/home/hmhuang/Work/Research_Assistant/Ilan_geothermal/'
sac_path = '/raid2/ILAN2022/Sorted_oneday/' + eq_date 
resp_path = '/raid2/ILAN2022/PZs/'
stations_file = '/raid2/ILAN2022/stations.csv'
'''
sensitivity_resp = 1
corner_freq = 5
damping = 0.7
gain = 3000
'''
sensitivity_resp = 306846
corner_freq = 5
damping = 0.7
gain = 1
component = ['DPZ']
pre_filt = (0.05, 0.1, 30.0, 50.0)

if os.path.exists(project_path + eq_date):
    shutil.rmtree(project_path + eq_date)
    os.mkdir(project_path + eq_date)
else:
    os.mkdir(project_path + eq_date)



# ============================================
# 找要跑ㄉ站
# 1. 先找最遠的，並把最遠的一站跟震央連起一線
# 2. 找這條線上所有的站
# ============================================
lon_sta = pd.read_csv(stations_file).Lon
lat_sta = pd.read_csv(stations_file).Lat
name_sta = pd.read_csv(stations_file).Station
dist_sta = []
for i in range(len(lon_sta)):   
    dist_sta.insert(i, (((lon_sta[i]-lon_epic)**2) + ((lat_sta[i]-lat_epic)**2))**0.5)
ind_max = np.argmax(dist_sta)
a = np.array([[lon_epic, 1],[lon_sta[ind_max],1]])
b = np.array([lat_epic, lat_sta[ind_max]])
p = np.linalg.solve(a, b)
diff_sta = []
for i in range(len(lon_sta)):  
    diff_sta.insert(i, p[0]*lon_sta[i]+p[1]-lat_sta[i])
index_sta_all = np.where(abs(np.array(diff_sta)) < th_line)
inthecase_sta = list(name_sta[index_sta_all[0]])

# ======
# 讀資料
# ======
for j in component:
    all_sta_that_day = []
    for k in range(len(inthecase_sta)):
        if len(glob.glob(f'{sac_path}/*{inthecase_sta[k]}*{j}*.SAC')) == 0:
            print('no data in ' + str(inthecase_sta[k]))
            continue
        else:
            all_sta_that_day.insert(k, glob.glob(f'{sac_path}/*{inthecase_sta[k]}*{j}*.SAC'))
    for i in all_sta_that_day:
        st = obspy.read(i[0])
        eq_time_comb = UTCDateTime(eq_date[0:4] + '-' + eq_date[4:6] + '-' + eq_date[6:8] + 'T' + eq_time[0:2] + ':' + eq_time[2:4] + ':' + eq_time[4:6])
        st.trim(eq_time_comb - 300, eq_time_comb  + 300)
        tr = st[0]
        # ============
        # 去除儀器響應
        # ============
        attach_paz(tr, resp_path + j)
        paz_1hz = corn_freq_2_paz(corner_freq, damp=damping)  
        paz_1hz['sensitivity'] = sensitivity_resp
        paz_1hz['gain'] = gain
        st.simulate(paz_remove=paz_1hz, pre_filt=pre_filt)
        # ===============================
        # dtrend, dmean, taper, filter
        # ===============================
        tr.detrend(type='demean')  
        tr.detrend(type='linear')
        tr.taper(max_percentage=0.05, type='cosine', max_length=len(tr[0].data), side='both')
        #tr.filter("highpass", freq=0.5) 
        # =======
        # 時頻圖
        # =======
        sig = tr.data
        time = tr.times("matplotlib")
        labelUTCtime1 = tr.stats.starttime
        time1 = date2num(labelUTCtime1.datetime)
        labelUTCtime2 = tr.stats.endtime
        time2 = date2num(labelUTCtime2.datetime)
        NFFT = 1024
        num_overlap = 512
        Fs = tr.stats.sampling_rate
        hspace = 0.7
        fig = plt.figure(figsize=(10,7))
        ax1 = plt.subplot2grid((3,1), (0,0), rowspan=1)
        ax1.plot(time, sig, 'k', linewidth=1)
        ax1.set_xlim(time1, time2)
        ax1.xaxis_date()
        fig.autofmt_xdate()
        ax1.grid(True)
        plt.title(eq_date[0:6] + '.' +  i[0][42:46] + "."  + j)
        plt.subplots_adjust(hspace=hspace)  
        plt.jet()
        ax2 = plt.subplot2grid((3,1), (1,0), rowspan=2)
        ax2.specgram(sig, Fs=Fs, NFFT=NFFT, noverlap=num_overlap,vmin=-250,vmax=-100, scale='dB')
        #ax2.set_yscale('symlog')
        cbar_x = ax2.get_position().x1 + 0.01 
        cbar_y = ax2.get_position().y0
        cbar_h = ax2.get_position().height
        cbar = fig.add_axes([cbar_x, cbar_y, 0.02, cbar_h])
        ax2.set_ylim([0, 50])
        ax2.set_xlabel("Time [sec]")
        ax2.set_ylabel("Frequency [Hz]")
        plt.colorbar(ax2.images[0], label='(dB)', cax=cbar) 
        fig.savefig(project_path + eq_date + '/' + i[0][42:46] + "."  + j + ".png")
        plt.close()
        plt.rcParams['font.sans-serif'] =  'Nimbus Roman'
        plt.show()
# %%
