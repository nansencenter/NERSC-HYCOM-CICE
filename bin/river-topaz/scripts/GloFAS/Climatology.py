import matplotlib.pyplot as plt
import netCDF4
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
mpl.rcParams.update({'font.size': 12}) #Pour changer taille des
import pandas as pd

months=["Annual mean","January", "February", "March",
       "April", "May", "June", "July", "August", "September",
       "October", "November", "December"]

#Creating empty dataframe
df=pd.DataFrame({'lon':[], 'lat':[], 'month':[], 'flux':[]})
# df.insert(0,"Name",[0,0,months[0]])
# new_row = {'lon':0, 'lat':0, 'month':months[0], 'flux':0}
# df.loc[len(df)] = new_row

cntr=0
with open('river/rivers_ahype-ehype_clim_rev2.dat','r') as f:

    for line in f:
        if line[0]=='T':
            L=line.split('  ')
            name=L[1]
            values=f.readline()
            values=values[5:-1].split(' ')
            coord=f.readline()
            coord=coord[4:-1].split(' ')
            lat,lon=coord[0].strip(),coord[-1].strip()
            cntr+=1

            for j in range(len(values)):
                new_row = {'lon': lon, 'lat': lat, 'month': months[j], 'flux': values[j]}
                df.loc[len(df)] = new_row
            print(cntr)
df.to_csv('river/climatology.csv', index=False)