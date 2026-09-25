import xarray as xr
import pandas as pd
from functions import *
import os
import sys


#ArcticGRO info
NAME=['Lena','Ob','Yenisei']

path = '/cluster/work/users/adracc/ArcticGRO' #Path to ArcticGRO discharge data (csv files taken from https://arcticgreatrivers.org/data/ )

def compute_ArtcicGRO_Dataframes(path_ArcticGRO_data,output_path):
    """
    Creates pandas Dataframes with ArcticGRO data.
    The first dataframe simply contains daily discharge in ArcticGRO observations for each river. With a column for each river and every row being different date
    The second dataframe contains monthly average for each river
    The third dataframe contains yearly average for each river

    These daily dataframe is used to compute the GloFAS correction
    """
    for river in NAME:
        try:
            river_file = [i for i in os.listdir(path_ArcticGRO_data) if os.path.isfile(os.path.join(path_ArcticGRO_data,i)) and \
                    i.startswith(river)][0]

            df_GRO=pd.read_excel(path_ArcticGRO_data+'/'+river_file).drop(columns=['agency', 'station_code', 'station_name','flag'])
            # print(df_GRO['date'].str.slice(stop=4))

            df_GRO['year']=df_GRO['date'].str.slice(stop=4)
            df_GRO['month']=df_GRO['date'].str.slice(start=5,stop=7)

            tmp_dis=df_GRO['discharge'] #saving it in m3/S
            df_GRO['discharge']=df_GRO['discharge']*24*3600*10**(-9) #discharge in km3/day


            aggregation_functions = {'discharge':'sum',  'date':'first', 'year':'first','river':'first'}
            df_GRO_yearly=df_GRO.groupby('year').aggregate(aggregation_functions)

            df_GRO['discharge']=df_GRO['discharge']*365/12 #discharge in km3/year
            aggregation_functions = {'discharge':'mean',  'date':'first', 'month':'first','river':'first'}
            df_GRO_monthly=df_GRO.groupby('month').aggregate(aggregation_functions)

            print(df_GRO_monthly)
            df_GRO['discharge']=tmp_dis #back in m3/s

            df_GRO_yearly.to_csv(output_path+'df_GRO_yearly_'+river)
            df_GRO_monthly.to_csv(output_path+'df_GRO_monthly_'+river)
            df_GRO.to_csv(output_path+'df_GRO_daily_'+river)
        except:
            print("No river named "+river)

compute_ArtcicGRO_Dataframes(path_ArcticGRO_data=path,output_path='Dataframes/')