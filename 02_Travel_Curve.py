# =====================================================
# 畫一個指定地震的走時曲線（地震與最遠站連線之經過的站點）
# 2nd Step for the Routine!!
# =====================================================

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
import math
# ========v
# 可改參數
# ========

eq_date = '20230105'
eq_time = '171138'
lon_epic = 122.01
lat_epic = 24.87
th_line = 0.003

th_multi = 0
trim_start = 0
trim_end = 30
xlim_plot = [0, 30]
text_loc = xlim_plot[1] + 1 
minimize_parm = 10
# ==========
# 比較不用改
# ==========
project_path = '/home/hmhuang/Work/Research_Assistant/Ilan_geothermal/'
sac_path = '/raid2/ILAN2022/Sorted_oneday/' + eq_date 
resp_path = '/raid2/ILAN2022/PZs/'
stations_file = '/raid2/ILAN2022/stations.csv'
sensitivity_resp = 306846
corner_freq = 5
damping = 0.7
gain = 1
component = ['DPZ']
pre_filt = (0.05, 0.1, 30.0, 50.0)



# ============================================
#                   找站
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
inthercase_sta_dist_km = 2*math.pi*6378.137*(np.array(dist_sta)[index_sta_all[0]])/360

d_LL = {'Lon': lon_sta[index_sta_all[0]], 'Lat': lat_sta[index_sta_all[0]]}
df_LL = pd.DataFrame(data=d_LL)
df_LL.to_csv(eq_date + '_event.txt', index=None, sep=' ')


# ============
# 抓資料路徑
# ============
d = {'sta': inthecase_sta, 'dist': inthercase_sta_dist_km}
df = pd.DataFrame(data=d)
df_resorted = df.sort_values('dist', ignore_index=True)
inthecase_sta = df_resorted.sta
inthercase_sta_dist_km = df_resorted.dist
#重新整理地震站及其與震央距離
fig = plt.figure(figsize=(8, 10))

inthecase_sta_new = []
inthercase_sta_dist_km_new = []

for j in component:
    all_sta_that_day = []
    for k in range(len(inthecase_sta)):
        if k == len(inthecase_sta)-1:
            continue
        else:
            if abs(inthercase_sta_dist_km[k+1]-inthercase_sta_dist_km[k])> th_multi:
                all_sta_that_day.insert(k, glob.glob(f'{sac_path}/*{inthecase_sta[k]}*{j}*.SAC'))
                inthecase_sta_new.insert(k, inthecase_sta[k])
                inthercase_sta_dist_km_new.insert(k, inthercase_sta_dist_km[k])
            else:
                continue
        #若距離差不多就不抓了
# ======
# 讀資料
# ======     
    for idx, i in enumerate(all_sta_that_day):
        #print(i)
        eq_time_comb = UTCDateTime(eq_date[0:4] + '-' + eq_date[4:6] + '-' + eq_date[6:8] + 'T' + eq_time[0:2] + ':' + eq_time[2:4] + ':' + eq_time[4:6])
        try:
            st = obspy.read(i[0])
            tr = st[0]
            st.trim(eq_time_comb + trim_start , eq_time_comb  + trim_end)
            #print(eq_time_comb)
            #st.write('/home/hmhuang/Work/Research_Assistant/Ilan_geothermal/outeq_' + eq_date + '/' + eq_date + '_' + i[0][42:46] + '.sac')
        except:
            continue    
        #print(st[0].stats.starttime)
        
        # ============
        # 去除儀器響應
        # ============
        pre_filt = (0.05, 0.1, 30.0, 50.0)     
        attach_paz(tr, resp_path + j)
        paz_1hz = corn_freq_2_paz(corner_freq, damp=damping)  
        paz_1hz['sensitivity'] = sensitivity_resp
        paz_1hz['gain'] = gain
        st.simulate(paz_remove=paz_1hz, pre_filt=pre_filt)
        # =====================
        # dtrend, dmean, taper
        # =====================
        tr.detrend(type='demean')  
        tr.detrend(type='linear')
        try:
            tr.taper(max_percentage=0.05, type='cosine', max_length=len(tr.data), side='both')
            tr.normalize()
            tr.filter('bandpass', freqmin=1, freqmax=10)
        except:
            pass
        
        y = tr.data/minimize_parm
        x = tr.times()
        plt.plot(x, (y+inthercase_sta_dist_km[idx]), 'k', linewidth=0.7)
        plt.text(text_loc, np.mean((y+inthercase_sta_dist_km[idx])), '(' + str(inthecase_sta[idx]) + ')', color='k', fontsize=10)

plt.rcParams['font.sans-serif'] =  'Nimbus Roman'
plt.ylabel('Distance [km]', fontsize=12)
plt.xlabel('Time [sec]' , fontsize=12)
plt.title('Event in ' + eq_date + ' (' + eq_time[0:2] + ':' + eq_time[2:4] + ':' + str(int(eq_time[4:6])+trim_start) + ')')
plt.xlim(xlim_plot)
plt.savefig('Event_' + eq_date + '.png', dpi = 300)
plt.show()
 # %%