from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date
import datetime
import warnings
import numpy as np
import argparse
import getpass
import os
from make_mask_for_generic_use import get_mask
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import abfile
import pickle
########################
warnings.filterwarnings('ignore')
user = getpass.getuser()
global REGmask
global REGmmask
global time
global dates
global yearmin
global yearmax
global model
global modelts
global workdir
global experiment
global region

# SOME STATISTICS #
def pc_bias(s,o):
         return 100.0*sum(s-o)/sum(o)
def bias(s,o):
         return np.mean(s-o)
def correlation(s,o):
         return np.corrcoef(o, s)[0,1]
def rmse(s,o):
         return np.sqrt(np.mean((s-o)**2))
def nstd(s,o):
         return np.std(s)/np.std(o)


def main(region_in,experiment_in,yearmin,yearmax,workdir_in):
    global REGmask
    global REGmmask
#    global time
#    global dates
#    global model
#    global modelts
#    global workdir
#    workdir="/nird/projects/NS9481K/caglar/model_archives"
    global variables
    global regions
    global workdir, region, experiment # Tell main to use the global ones
    workdir = workdir_in
    region = region_in
    experiment = experiment_in

    if region == "TP2":
        f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskTP2.pckl','rb')
        REGmmask = pickle.load(f,encoding="latin1")
        f.close()
        region = "TP2a0.10"
    if region == "TP5":
        f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskTP5.pckl','rb')
        REGmmask = pickle.load(f,encoding="latin1")
        f.close()
        region = "TP5a0.06"
    if region == "NAT":
        f=open('/cluster/projects/nn9481k/BGC.Validation/OMmaskNAT.pckl','rb')
        REGmmask = pickle.load(f,encoding="latin1")
        f.close()       
        region = "NATa1.00"

    experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]

    # https://doi.org/10.5194/bg-19-1087-2022
    obsfile = "/cluster/projects/nn9481k/BGC.Validation/surface.carbon.data.nc"
    nc = NetCDFFile(obsfile)

    lat = nc.variables["latitude"][:]
    lon = nc.variables["longitude"][:]
    nctime = nc['time']
    time = num2date(nctime[:], units=nctime.units)
    indexmin = np.abs(time -datetime.datetime(yearmin,1,15) ).argmin()
    indexmax = np.abs(time -datetime.datetime(yearmax,12,15) ).argmin()
    time = time[indexmin:indexmax+1]
    dates = mdates.date2num(time)

    fgco2 = nc.variables["fgco2"][indexmin:indexmax+1,:,:]
    std_fgco2 = nc.variables["fgco2_uncertainty"][indexmin:indexmax+1,:,:]
    ph = nc.variables["ph"][indexmin:indexmax+1,:,:]
    std_ph = nc.variables["ph_uncertainty"][indexmin:indexmax+1,:,:]
    talk = nc.variables["talk"][indexmin:indexmax+1,:,:]
    std_talk = nc.variables["talk_uncertainty"][indexmin:indexmax+1,:,:]
    tco2 = nc.variables["tco2"][indexmin:indexmax+1,:,:]
    std_tco2 = nc.variables["tco2_uncertainty"][indexmin:indexmax+1,:,:]

    xx,yy = np.meshgrid(lon,lat) 
    dummy = get_mask(xx,yy)
    REGmask = dummy

    # obsts = {
    #           "FGCO2": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "stdFGCO2": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "ALK": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "stdALK": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "DIC": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "stdDIC": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "PH": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]},
    #           "stdPH": {"LFTB":[],"NRWB":[],"BRTS":[],"IBIB":[],"GSIS":[]}
    # }

    variables = ["FGCO2", "stdFGCO2", "ALK", "stdALK", "DIC", "stdDIC", "PH", "stdPH"]
    regions = ["LFTB", "NRWB", "BRTS", "IBIB", "GSIS"]

    obsts = {v: {r: [] for r in regions} for v in variables}

    # 1. Map your dictionary keys to the actual data variables
    data_map = {
      "FGCO2": fgco2, "PH": ph, "ALK": talk, "DIC": tco2,
      "stdFGCO2": std_fgco2, "stdPH": std_ph, "stdALK": std_talk, "stdDIC": std_tco2
    }

    # 2. Define the regions you want to process
    target_regions = regions

    # 3. Use a nested loop to populate the dictionary
    for key, data_var in data_map.items():
      for regionn in target_regions:
        obsts[key][regionn] = get_timeseries(data_var, regionn)


    # obsts["FGCO2"]["Barents"] = get_timeseries(fgco2,"Barents")
    # obsts["FGCO2"]["NorwegianN"] = get_timeseries(fgco2,"NorwegianN")
    # obsts["FGCO2"]["NorwegianS"] = get_timeseries(fgco2,"NorwegianS")
    # obsts["stdFGCO2"]["Barents"] = get_timeseries(std_fgco2,"Barents")
    # obsts["stdFGCO2"]["NorwegianN"] = get_timeseries(std_fgco2,"NorwegianN")
    # obsts["stdFGCO2"]["NorwegianS"] = get_timeseries(std_fgco2,"NorwegianS")

    # obsts["PH"]["Barents"] = get_timeseries(ph,"Barents")
    # obsts["PH"]["NorwegianN"] = get_timeseries(ph,"NorwegianN")
    # obsts["PH"]["NorwegianS"] = get_timeseries(ph,"NorwegianS")
    # obsts["stdPH"]["Barents"] = get_timeseries(std_ph,"Barents")
    # obsts["stdPH"]["NorwegianN"] = get_timeseries(std_ph,"NorwegianN")
    # obsts["stdPH"]["NorwegianS"] = get_timeseries(std_ph,"NorwegianS")

    # obsts["ALK"]["Barents"] = get_timeseries(talk,"Barents")
    # obsts["ALK"]["NorwegianN"] = get_timeseries(talk,"NorwegianN")
    # obsts["ALK"]["NorwegianS"] = get_timeseries(talk,"NorwegianS")
    # obsts["stdALK"]["Barents"] = get_timeseries(std_talk,"Barents")
    # obsts["stdALK"]["NorwegianN"] = get_timeseries(std_talk,"NorwegianN")
    # obsts["stdALK"]["NorwegianS"] = get_timeseries(std_talk,"NorwegianS")

    # obsts["DIC"]["Barents"] = get_timeseries(tco2,"Barents")
    # obsts["DIC"]["NorwegianN"] = get_timeseries(tco2,"NorwegianN")
    # obsts["DIC"]["NorwegianS"] = get_timeseries(tco2,"NorwegianS")
    # obsts["stdDIC"]["Barents"] = get_timeseries(std_tco2,"Barents")
    # obsts["stdDIC"]["NorwegianN"] = get_timeseries(std_tco2,"NorwegianN")
    # obsts["stdDIC"]["NorwegianS"] = get_timeseries(std_tco2,"NorwegianS")

    # Creates the dictionary with your predefined variables and regions
    modelts = {v: {r: [] for r in regions} for v in variables}

    # modelts = {
    #           "FGCO2": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "stdFGCO2": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "ALK": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "stdALK": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "DIC": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "stdDIC": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "PH": {"Barents":[],"NorwegianN":[],"NorwegianS":[]},
    #           "stdPH": {"Barents":[],"NorwegianN":[],"NorwegianS":[]}
    # }

    nameout = workdir + user + "/" + \
        region + "/" + experiment + "/data/carbon_"+str(yearmin)+"_"+str(yearmax)+".pckl"

    modeldata = get_model_monthly_averages(modelts,yearmin,yearmax)
    ff = open(nameout,"wb")
    pickle.dump([modeldata,obsts,dates],ff)
    ff.close()    

#    for var in ["FGCO2","ALK","DIC","PH"]:
#	for subregion in ["Barents","NorwegianN","NorwegianS"]:
#		plot_timeseries(obsts[var][subregion],obsts["std"+var][subregion],modeldata[var][subregion],modeldata["std"+var][subregion],varname=subregion+"_"+var)	
#		fig = plt.gcf()
#	        fig.canvas.print_figure("/nird/projects/NS9481K/caglar/model_archives/TP5a0.06/expt_03.1/data/CO2_"+var+"_"+subregion+".png",dpi=300)	
#                plt.close(fig)

def clear_model():

    model = {v: {r: [] for r in regions} for v in variables}

    return model

#def take_model_average(modelts,model):
#
#    for variables in ["FGCO2","ALK","DIC","PH"]:
#        for subregion in ["Barents","NorwegianN","NorwegianS"]:
#            modelts[variables][subregion].append( np.array(model[variables][subregion]).mean() )
#            if variables == "PH":
#               m = np.array( model[variables][subregion] )
#               m = np.ma.masked_where(m>1E6,m)
#               print(subregion,m.mean())
#            modelts["std"+variables][subregion].append( np.std( np.array(model[variables][subregion])) )
#    return modelts

def take_model_average(modelts,model):

    for variables in ["FGCO2","ALK","DIC","PH"]:
        for subregion in regions:
            m = np.array( model[variables][subregion] )
            m = np.ma.masked_where(np.abs(m)>1E6,m)
            if variables == "FGCO2": m = m/1000.*365.
            modelts[variables][subregion].append( m.mean() )
            modelts["std"+variables][subregion].append( np.std(m)) 
    return modelts

def get_model_monthly_averages(modelts,yearmin,yearmax):
    datem = datetime.datetime(yearmin,1,1)
    month = 1
    model = clear_model()
    while datem <= datetime.datetime(yearmax,12,31):
        nday = datem.strftime('%j')
        if month == datem.month:
            if datem == datetime.datetime(yearmax,12,31):
               modelts = take_model_average(modelts,model)
            else:        
               model = model_append(datem.year,nday,model)
        else:
            month = datem.month
            modelts = take_model_average(modelts,model)
            model = clear_model()
            model = model_append(datem.year,nday,model)
#        print(nday,month,datem.month)
        datem = datem + datetime.timedelta(days=1)
    return modelts

def model_append(year,nday,model):

    fair, ph, dic, alk = get_model_day(year,nday)
    for subregion in regions:
        model["FGCO2"][subregion].append( fair[ REGmmask[subregion] == 1 ])
        model["DIC"][subregion].append( dic[ REGmmask[subregion] == 1 ])
        model["ALK"][subregion].append( alk[ REGmmask[subregion] == 1 ])
        model["PH"][subregion].append( ph[ REGmmask[subregion] == 1 ])
    return model

def get_model_day(year,nday):
    global workdir, region, experiment
    print("archm."+str(year)+"_"+str(nday).zfill(3)+"_12.b")
#    f = abfile.ABFileArchv(workdir + user + "/" + \
#                  region + "/" + experiment + \
#                    "/data/"+"archm."+str(year)+"_"+str(nday).zfill(3)+"_12.b","r")
    f = abfile.ABFileArchv(workdir + user + "/" + \
                  region + "/" + experiment + \
                    "/data/"+"archm."+str(year)+"_"+str(nday).zfill(3)+"_12.b","r")

    fair = f.read_field('CO2_fair',0)
    ph = f.read_field('CO2_pH',1)
    dic = f.read_field('CO2_c',1)
    alk = f.read_field('CO2_TA',1)

    return fair, ph, dic, alk

def get_timeseries(varib,region):
    outvarib = []
    for t in range(varib.shape[0]):
        outvarib.append( varib[t,:,:][ REGmask[region]==1 ].mean() ) 
    return outvarib

def get_title(varname):
    if "FGCO2" in varname:
       title = varname.split("_")[0]+" surface downward mass flux of CO2 expressed as carbon" 
       ytitle = r"molC m$^{-2}$ yr$^{-1}$"
    if "PH" in varname:
       title = varname.split("_")[0]+" pH"
       ytitle = r" "
    if "ALK" in varname:
       title = varname.split("_")[0]+" total alkalinity in surface seawater"
       ytitle = r"mmolC m$^{-3}$"
    if "DIC" in varname:
       title = varname.split("_")[0]+" surface ocean dissolved inorganic carbon"
       ytitle = r"mmolC m$^{-3}$"
    else:
       title=""
       ytitle=""
    return ytitle,title

# def plot_timeseries(varib,std,mvarib,mstd, *args, **kwargs):
#     varib = np.array(varib)
#     std = np.array(std)
#     mvarib = np.array(mvarib)
#     mstd = np.array(mstd)
#     fig = plt.figure(figsize=(10,3.5),facecolor='w')
#     ax  = fig.add_subplot(1,1,1)
#     plt.plot(dates,varib,color='#1f77b4',label="OBS REP L4")
#     plt.fill_between(dates,varib+std,varib-std,color='#1f77b4',alpha=0.35)
#     plt.plot(dates,mvarib,color='#ff7f0e',label="model")
#     plt.fill_between(dates,mvarib+mstd,mvarib-mstd,color='#ff7f0e',alpha=0.35)
#     ax.xaxis.axis_date()
#     ax.set_xlim(datetime.datetime(time[0].year,1,1),datetime.datetime(time[-1].year+1,1,1) )
#     ax.yaxis.set_tick_params(labelsize='large')
#     ax.xaxis.set_tick_params(labelsize='large')
#     plt.xticks(rotation = 25)
#     varname = kwargs.get('varname',None)
#     if varname is not None:
#         ytitle,title = get_title(varname) 
#         ax.set_ylabel(ytitle,fontsize=13)
#         plt.title(title,fontsize=14)
#         printname = varname+"_"+str(yearmin)+"_"+str(yearmax)+".png"
#     ax.grid(color='k', linestyle=':', linewidth=1)
#     b = bias(mvarib,varib)
#     e = rmse(mvarib,varib)
#     #ax.text(0.025, 0.95, 'bias: '+str(round(b,4)), fontsize=10, transform=ax.transAxes)
#     #ax.text(0.025, 0.9, 'rmse : '+str(round(e,4)), fontsize=10, transform=ax.transAxes)
#     print(varname+"  bias: "+str(round(b,4)))
#     print(varname+"  error: "+str(round(e,4)))
#     plt.legend()
#     plt.tight_layout()


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('yearmin',  type=int, help='e.g 2010')
    parser.add_argument('yearmax',  type=int, help='e.g 2010')
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.yearmin,args.yearmax,args.workdir)
