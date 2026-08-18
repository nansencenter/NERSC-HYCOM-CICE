import numpy as np
from netCDF4 import Dataset as NetCDFFile
import pickle
import abplot
import matplotlib.pyplot as plt
import calendar
import datetime
from dateutil.relativedelta import relativedelta
import argparse
import warnings
import getpass
warnings.filterwarnings('ignore')
user = getpass.getuser()

# python model_integrates.py TP2 070 1990 2020 --depth=200

def main(region,experiment,date1,date2,depth,workdir):

    workdir = workdir + "/" + user + "/"
    expt = 'expt_'+experiment[0:2]+'.'+experiment[2]
    
    if region == "TP2":
        region = "TP2a0.10" 
        f=open('/cluster/home/cagyum/codes/VALIDATION/maskTP2.pckl','rb')
    if region == "TP5":
        region = "TP5a0.06" 
        f=open('/cluster/home/cagyum/codes/VALIDATION/maskTP5.pckl','rb')
    if region == "TP0":
        region = "TP0a1.00"
        f=open('/cluster/home/cagyum/codes/VALIDATION/maskTP0.pckl','rb')

    REGmask = pickle.load(f,encoding="latin1")
    f.close()

    depthm = abplot.getdepth(region=region)
    dataf = workdir + "/" + region + "/" + expt + "/data/"

#    NIT={}; PHO={}; SIL={}; CHL={}; OXY={}
    depths = ['0-10','10-30','30-50','50-100']
#    for d in depths:
#      for m in range(17):
#        NIT[d][m]=[np.zeros((depthm.shape))]; PHO[d][m]=[np.zeros((depthm.shape))]; SIL[d][m]=[np.zeros((depthm.shape))]; CHL[d][m]=[np.zeros((depthm.shape))]; OXY[d][m]=[np.zeros((depthm.shape))]; 

    NIT = {d: {m: np.zeros(depthm.shape) for m in range(17)} for d in depths}
    PHO = {d: {m: np.zeros(depthm.shape) for m in range(17)} for d in depths}
    SIL = {d: {m: np.zeros(depthm.shape) for m in range(17)} for d in depths}
    CHL = {d: {m: np.zeros(depthm.shape) for m in range(17)} for d in depths}
    OXY = {d: {m: np.zeros(depthm.shape) for m in range(17)} for d in depths}


    nyear = int(date2) - int(date1) + 1
    datenow = datetime.datetime(date1,1,1)
    dateend = datetime.datetime(date2,12,31)  
    while datenow <= dateend :
        print(datenow)
        nyearday = 366 if calendar.isleap(datenow.year) else 365
        _,number_of_days = calendar.monthrange(datenow.year, datenow.month)
        nday = datenow.timetuple().tm_yday
        season_days = 0.0
        if int(datenow.month) in [12, 1, 2]: 
            season = 13
            _,s_days = calendar.monthrange(datenow.year-1, 12)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 1)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 2)
            season_days = season_days + s_days
        if int(datenow.month) in [3, 4, 5]:
            season = 14
            _,s_days = calendar.monthrange(datenow.year, 3)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 4)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 5)
            season_days = season_days + s_days
        if int(datenow.month) in [6, 7, 8]:
            season = 15
            _,s_days = calendar.monthrange(datenow.year, 6)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 7)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 8)
            season_days = season_days + s_days
        if int(datenow.month) in [9, 10, 11]:
            season = 16
            _,s_days = calendar.monthrange(datenow.year, 9)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 10)
            season_days = season_days + s_days
            _,s_days = calendar.monthrange(datenow.year, 11)
            season_days = season_days + s_days

        arcf = dataf + "archm." + str(datenow.year) + "_" + str(nday).zfill(3) + "_12.a"
        arc  = abplot.openfile(arcf)


        intg10   = abplot.integr(arc,'no3',10)
        intg30   = abplot.integr(arc,'no3',30)
        intg50   = abplot.integr(arc,'no3',50)
        intg100  = abplot.integr(arc,'no3',100)
        intg100 = (intg100 - intg50)/50.0
        intg50 = (intg50 - intg30)/20.0
        intg30 = (intg30 - intg10)/20.0
        intg10 = intg10/10.0

        NIT['0-10'][datenow.month - 1] = NIT['0-10'][datenow.month - 1] + intg10/number_of_days/nyear
        NIT['0-10'][12] = NIT['0-10'][12] + intg10/nyearday/nyear
        NIT['0-10'][season] = NIT['0-10'][season] + intg10/season_days/nyear
        NIT['10-30'][datenow.month - 1] = NIT['10-30'][datenow.month - 1] + intg30/number_of_days/nyear
        NIT['10-30'][12] = NIT['10-30'][12] + intg30/nyearday/nyear
        NIT['10-30'][season] = NIT['10-30'][season] + intg30/season_days/nyear
        NIT['30-50'][datenow.month - 1] = NIT['30-50'][datenow.month - 1] + intg50/number_of_days/nyear
        NIT['30-50'][12] = NIT['30-50'][12] + intg50/nyearday/nyear
        NIT['30-50'][season] = NIT['30-50'][season] + intg50/season_days/nyear
        NIT['50-100'][datenow.month - 1] = NIT['50-100'][datenow.month - 1] + intg100/number_of_days/nyear
        NIT['50-100'][12] = NIT['50-100'][12] + intg100/nyearday/nyear
        NIT['50-100'][season] = NIT['50-100'][season] + intg100/season_days/nyear

        intg10   = abplot.integr(arc,'pho',10)
        intg30   = abplot.integr(arc,'pho',30)
        intg50   = abplot.integr(arc,'pho',50)
        intg100  = abplot.integr(arc,'pho',100)
        intg100 = (intg100 - intg50)/50.0
        intg50 = (intg50 - intg30)/20.0
        intg30 = (intg30 - intg10)/20.0
        intg10 = intg10/10.0

        PHO['0-10'][datenow.month - 1] = PHO['0-10'][datenow.month - 1] + intg10/number_of_days/nyear
        PHO['0-10'][12] = PHO['0-10'][12] + intg10/nyearday/nyear
        PHO['0-10'][season] = PHO['0-10'][season] + intg10/season_days/nyear
        PHO['10-30'][datenow.month - 1] = PHO['10-30'][datenow.month - 1] + intg30/number_of_days/nyear
        PHO['10-30'][12] = PHO['10-30'][12] + intg30/nyearday/nyear
        PHO['10-30'][season] = PHO['10-30'][season] + intg30/season_days/nyear
        PHO['30-50'][datenow.month - 1] = PHO['30-50'][datenow.month - 1] + intg50/number_of_days/nyear
        PHO['30-50'][12] = PHO['30-50'][12] + intg50/nyearday/nyear
        PHO['30-50'][season] = PHO['30-50'][season] + intg50/season_days/nyear
        PHO['50-100'][datenow.month - 1] = PHO['50-100'][datenow.month - 1] + intg100/number_of_days/nyear
        PHO['50-100'][12] = PHO['50-100'][12] + intg100/nyearday/nyear
        PHO['50-100'][season] = PHO['50-100'][season] + intg100/season_days/nyear

        intg10   = abplot.integr(arc,'sil',10)
        intg30   = abplot.integr(arc,'sil',30)
        intg50   = abplot.integr(arc,'sil',50)
        intg100  = abplot.integr(arc,'sil',100)
        intg100 = (intg100 - intg50)/50.0
        intg50 = (intg50 - intg30)/20.0
        intg30 = (intg30 - intg10)/20.0
        intg10 = intg10/10.0

        SIL['0-10'][datenow.month - 1] = SIL['0-10'][datenow.month - 1] + intg10/number_of_days/nyear
        SIL['0-10'][12] = SIL['0-10'][12] + intg10/nyearday/nyear
        SIL['0-10'][season] = SIL['0-10'][season] + intg10/season_days/nyear
        SIL['10-30'][datenow.month - 1] = SIL['10-30'][datenow.month - 1] + intg30/number_of_days/nyear
        SIL['10-30'][12] = SIL['10-30'][12] + intg30/nyearday/nyear
        SIL['10-30'][season] = SIL['10-30'][season] + intg30/season_days/nyear
        SIL['30-50'][datenow.month - 1] = SIL['30-50'][datenow.month - 1] + intg50/number_of_days/nyear
        SIL['30-50'][12] = SIL['30-50'][12] + intg50/nyearday/nyear
        SIL['30-50'][season] = SIL['30-50'][season] + intg50/season_days/nyear
        SIL['50-100'][datenow.month - 1] = SIL['50-100'][datenow.month - 1] + intg100/number_of_days/nyear
        SIL['50-100'][12] = SIL['50-100'][12] + intg100/nyearday/nyear
        SIL['50-100'][season] = SIL['50-100'][season] + intg100/season_days/nyear

        intg10   = abplot.integr(arc,'ECO_oxy',10)
        intg30   = abplot.integr(arc,'ECO_oxy',30)
        intg50   = abplot.integr(arc,'ECO_oxy',50)
        intg100  = abplot.integr(arc,'ECO_oxy',100)
        intg100 = (intg100 - intg50)/50.0
        intg50 = (intg50 - intg30)/20.0
        intg30 = (intg30 - intg10)/20.0
        intg10 = intg10/10.0

        OXY['0-10'][datenow.month - 1] = OXY['0-10'][datenow.month - 1] + intg10/number_of_days/nyear
        OXY['0-10'][12] = OXY['0-10'][12] + intg10/nyearday/nyear
        OXY['0-10'][season] = OXY['0-10'][season] + intg10/season_days/nyear
        OXY['10-30'][datenow.month - 1] = OXY['10-30'][datenow.month - 1] + intg30/number_of_days/nyear
        OXY['10-30'][12] = OXY['10-30'][12] + intg30/nyearday/nyear
        OXY['10-30'][season] = OXY['10-30'][season] + intg30/season_days/nyear
        OXY['30-50'][datenow.month - 1] = OXY['30-50'][datenow.month - 1] + intg50/number_of_days/nyear
        OXY['30-50'][12] = OXY['30-50'][12] + intg50/nyearday/nyear
        OXY['30-50'][season] = OXY['30-50'][season] + intg50/season_days/nyear
        OXY['50-100'][datenow.month - 1] = OXY['50-100'][datenow.month - 1] + intg100/number_of_days/nyear
        OXY['50-100'][12] = OXY['50-100'][12] + intg100/nyearday/nyear
        OXY['50-100'][season] = OXY['50-100'][season] + intg100/season_days/nyear

        intg10   = abplot.integr(arc,'chl',10)
        intg30   = abplot.integr(arc,'chl',30)
        intg50   = abplot.integr(arc,'chl',50)
        intg100  = abplot.integr(arc,'chl',100)
        intg100 = (intg100 - intg50)/50.0
        intg50 = (intg50 - intg30)/20.0
        intg30 = (intg30 - intg10)/20.0
        intg10 = intg10/10.0

        CHL['0-10'][datenow.month - 1] = CHL['0-10'][datenow.month - 1] + intg10/number_of_days/nyear
        CHL['0-10'][12] = CHL['0-10'][12] + intg10/nyearday/nyear
        CHL['0-10'][season] = CHL['0-10'][season] + intg10/season_days/nyear
        CHL['10-30'][datenow.month - 1] = CHL['10-30'][datenow.month - 1] + intg30/number_of_days/nyear
        CHL['10-30'][12] = CHL['10-30'][12] + intg30/nyearday/nyear
        CHL['10-30'][season] = CHL['10-30'][season] + intg30/season_days/nyear
        CHL['30-50'][datenow.month - 1] = CHL['30-50'][datenow.month - 1] + intg50/number_of_days/nyear
        CHL['30-50'][12] = CHL['30-50'][12] + intg50/nyearday/nyear
        CHL['30-50'][season] = CHL['30-50'][season] + intg50/season_days/nyear
        CHL['50-100'][datenow.month - 1] = CHL['50-100'][datenow.month - 1] + intg100/number_of_days/nyear
        CHL['50-100'][12] = CHL['50-100'][12] + intg100/nyearday/nyear
        CHL['50-100'][season] = CHL['50-100'][season] + intg100/season_days/nyear

        datenow = datenow + relativedelta(days=1)

    savef = dataf + "python_model_averages_"+str(date1)+"_"+str(date2)+".pckl"
    openf=open(savef,"wb")
    pickle.dump([NIT,PHO,SIL,OXY,CHL],openf)
    openf.close()
#     

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('date1',  type=int, help='e.g 2010')
    parser.add_argument('date2',  type=int, help='e.g 2010')
    parser.add_argument('--depth',  type=float, default=200)
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.date1,args.date2,depth=args.depth,workdir=args.workdir)
