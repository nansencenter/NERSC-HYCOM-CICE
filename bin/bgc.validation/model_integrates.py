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

    DIA={}; FLA={}; MIC={}; MES={}; MES750={}; PP={}; CHL={}; CCL={}; DOC={}; DET={}; DOC100={}; DET100={}
    for m in range(17):
        DIA[m]=[np.zeros((depthm.shape))]; FLA[m]=[np.zeros((depthm.shape))]; MIC[m]=[np.zeros((depthm.shape))]; MES[m]=[np.zeros((depthm.shape))]; MES750[m]=[np.zeros((depthm.shape))]; PP[m]=[np.zeros((depthm.shape))]; CHL[m]=[np.zeros((depthm.shape))]; CCL[m]=[np.zeros((depthm.shape))]; DOC[m]=[np.zeros((depthm.shape))]; DOC100[m]=[np.zeros((depthm.shape))]; DET[m]=[np.zeros((depthm.shape))]; DET100[m]=[np.zeros((depthm.shape))]; 

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

        intg = abplot.integr(arc,'ECO_dia',depth)
        DIA[datenow.month - 1] = DIA[datenow.month - 1] + intg/number_of_days/nyear
        DIA[12] = DIA[12] + intg/nyearday/nyear
        DIA[season] = DIA[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_fla',depth)
        FLA[datenow.month - 1] = FLA[datenow.month - 1] + intg/number_of_days/nyear
        FLA[12] = FLA[12] + intg/nyearday/nyear
        FLA[season] = FLA[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_ccl',depth)
        CCL[datenow.month - 1] = CCL[datenow.month - 1] + intg/number_of_days/nyear
        CCL[12] = CCL[12] + intg/nyearday/nyear
        CCL[season] = CCL[season] + intg/season_days/nyear

#        intg = abplot.getvarib(arc,'ECO_prim',depth) # getvarib already integrates PP within 200 meters 
#        PP[datenow.month - 1] = PP[datenow.month - 1] + (intg*86400.)/number_of_days/nyear
#        PP[12] = PP[12] + intg/nyearday/nyear
#        PP[season] = PP[season] + intg/season_days/nyear
        intg = abplot.integr(arc,'ECO_prim',depth) # getvarib already integrates PP within 200 meters 
        PP[datenow.month - 1] = PP[datenow.month - 1] + (intg*86400.)/number_of_days/nyear
        PP[12] = PP[12] + intg/nyearday/nyear
        PP[season] = PP[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'total_ch',10)
        CHL[datenow.month - 1] = CHL[datenow.month - 1] + intg/10./number_of_days/nyear # I average upper 10 meters
        CHL[12] = CHL[12] + intg/nyearday/nyear
        CHL[season] = CHL[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_meso',depth)
        MES[datenow.month - 1] = MES[datenow.month - 1] + intg/number_of_days/nyear
        MES[12] = MES[12] + intg/nyearday/nyear
        MES[season] = MES[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_meso',750)
        MES750[datenow.month - 1] = MES750[datenow.month - 1] + intg/number_of_days/nyear
        MES750[12] = MES750[12] + intg/nyearday/nyear
        MES750[season] = MES750[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_micr',depth)
        MIC[datenow.month - 1] = MIC[datenow.month - 1] + intg/number_of_days/nyear
        MIC[12] = MIC[12] + intg/nyearday/nyear 
        MIC[season] = MIC[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_dom',depth)
        DOC[datenow.month - 1] = DOC[datenow.month - 1] + intg/number_of_days/nyear
        DOC[12] = DOC[12] + intg/nyearday/nyear
        DOC[season] = DOC[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_dom',100)
        DOC100[datenow.month - 1] = DOC100[datenow.month - 1] + intg/number_of_days/nyear
        DOC100[12] = DOC100[12] + intg/nyearday/nyear
        DOC100[season] = DOC100[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_det',depth)
        DET[datenow.month - 1] = DET[datenow.month - 1] + intg/number_of_days/nyear
        DET[12] = DET[12] + intg/nyearday/nyear
        DET[season] = DET[season] + intg/season_days/nyear

        intg = abplot.integr(arc,'ECO_det',100)
        DET100[datenow.month - 1] = DET100[datenow.month - 1] + intg/number_of_days/nyear
        DET100[12] = DET100[12] + intg/nyearday/nyear
        DET100[season] = DET100[season] + intg/season_days/nyear

        datenow = datenow + relativedelta(days=1)

    savef = dataf + "python_model_integrates_"+str(date1)+"_"+str(date2)+".pckl"
    openf=open(savef,"wb")
    pickle.dump([DIA,FLA,CCL,MIC,MES,MES750,PP,CHL,DOC,DOC100,DET,DET100],openf)
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
