# %% 
#
#  3rd Step for the Routine!!
# 
# %% 
import pygmt
import pandas as pd 

river_file = "riverpoly/riverpoly.shp"
project_path = "/home/hmhuang/Work/Research_Project/Shihmen/"

stations_file = '/raid2/ILAN2022/stations.csv'



depth_epic = 4
focal_mechanism = dict(strike=231, dip=31, rake=-95, magnitude=4.11)

eq_date = '20230105'
eq_time = '171138'
lon_epic = 122.01
lat_epic = 24.87


lon = pd.read_csv(stations_file).Lon
lat = pd.read_csv(stations_file).Lat
sta = pd.read_csv(stations_file).Station	
lon_used = pd.read_table(eq_date + '_event.txt', sep=' ').Lon
#print(lon_used)

lat_used = pd.read_table(eq_date + '_event.txt', sep=' ').Lat
#print(lat_used)

Region = [min(lon)-0.003, max(lon)+0.003, min(lat)-0.003, max(lat)+0.003]
grid = pygmt.datasets.load_earth_relief(resolution="01s", region=Region)
dgrid = pygmt.grdgradient(grid=grid, azimuth=300)

with pygmt.config(FONT="Times-Roman"):
    fig = pygmt.Figure()
    fig.grdimage(grid=grid, projection="M15c", frame='a', cmap="gray",transparency=30)
    fig.plot(x = lon, y = lat, pen='faint', style='t7p', color = 'yellow', label='Stations')
    fig.plot(x = lon_used, y = lat_used, pen='faint', style='t7p', color = 'green', label='Stations in Plot')
    #fig.plot(x = lon_epic, y = lat_epic, pen='faint', style='c10p', color = 'red', label='Source')


    with fig.inset(
    position="jBL+o0.5c/0.2c",
    box="+pblack",
    #[121.397, 122.1, 24.459, 24.9]
    #[121.397, 122.7, 24.459, 25.1]
    region=[121.397, 122.7, 24.459, 25.1],
    projection="M3c",):
        fig.coast(
            land="gray",
            shorelines=True,
            water="white",
        )
        rectangle = [[Region[0], Region[2], Region[1], Region[3]]]
        fig.plot(data=rectangle, style="r+s", pen="0.02c,red")
        fig.meca(focal_mechanism, scale="0.3c", longitude=lon_epic, latitude=lat_epic, depth=depth_epic)
    fig.text(x=121.699, y=24.7035, text='Lanyang River', font='20p,Times-BoldItalic, white')
    fig.legend(position='JBR+jBR+o0.2c+w3/1', box='+gwhite+p1p+r')
    #fig.legend(position='JBR+jBR+o0.2c+w3/1.5', box='+gwhite+p1p+r')
    fig.show()
    #fig.savefig("my-figure.png")
# %%
