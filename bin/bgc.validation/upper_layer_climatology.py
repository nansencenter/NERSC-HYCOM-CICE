import abplot
import argparse
import warnings
import getpass
import pickle
warnings.filterwarnings('ignore')
user = getpass.getuser()

def main(region,experiment,date1,date2,depth,workdir):

    workdir = workdir + "/" + user + "/"
    expt = 'expt_'+experiment[0:2]+'.'+experiment[2]
    
    if region == "TP2":
        region = "TP2a0.10" 
    if region == "TP5":
        region = "TP5a0.06" 
    if region == "TP0":
        region = "TP0a1.00"
    if region == "NAT":
        region = "NATa1.00"

    dataf = workdir + "/" + region + "/" + expt + "/data/"
    print(dataf)

    NIT = {"0-10":{},'10-30':{},"30-50":{},"50-100":{}}
    SIL = {"0-10":{},'10-30':{},"30-50":{},"50-100":{}}
    PHO = {"0-10":{},'10-30':{},"30-50":{},"50-100":{}}
    CHL = {"0-10":{},'10-30':{},"30-50":{},"50-100":{}}
    OXY = {"0-10":{},'10-30':{},"30-50":{},"50-100":{}}

    for m in range(1,13):

      arcf = dataf + "CLIM." + str(m).zfill(2) + "." + str(date1) + "." + str(date2) + ".a"
      arc  = abplot.openfile(arcf)
      print("month:",m,"file:",arcf)
      # these below first integrate the total depth range, then remove the depth range above, and divide by the remaining depth to take average

      for var in ['no3', 'pho', 'sil', 'ECO_oxy', 'total_ch']:

         intg10   = abplot.integr(arc,var,10)
         intg30   = abplot.integr(arc,var,30)
         intg50   = abplot.integr(arc,var,50)
         intg100  = abplot.integr(arc,var,100)
         intg100 = (intg100 - intg50)/50.0
         intg50 = (intg50 - intg30)/20.0
         intg30 = (intg30 - intg10)/20.0
         intg10 = intg10/10.0

         if var == 'no3'     :  NIT["0-10"][m] = intg10; NIT["10-30"][m] = intg30; NIT["30-50"][m] = intg50; NIT["50-100"][m] = intg100;
         if var == 'sil'     :  SIL["0-10"][m] = intg10; SIL["10-30"][m] = intg30; SIL["30-50"][m] = intg50; SIL["50-100"][m] = intg100;
         if var == 'pho'     :  PHO["0-10"][m] = intg10; PHO["10-30"][m] = intg30; PHO["30-50"][m] = intg50; PHO["50-100"][m] = intg100;
         if var == 'total_ch'     :  CHL["0-10"][m] = intg10; CHL["10-30"][m] = intg30; CHL["30-50"][m] = intg50; CHL["50-100"][m] = intg100;
         if var == 'ECO_oxy' :  OXY["0-10"][m] = intg10; OXY["10-30"][m] = intg30; OXY["30-50"][m] = intg50; OXY["50-100"][m] = intg100;

    arcf = dataf + "CLIM." + str(date1) + "." + str(date2) + ".a"
    arc  = abplot.openfile(arcf)
    print("annual climatology file:", arcf)

    intg10   = abplot.integr(arc,var,10)
    intg30   = abplot.integr(arc,var,30)
    intg50   = abplot.integr(arc,var,50)
    intg100  = abplot.integr(arc,var,100)
    intg100 = (intg100 - intg50)/50.0
    intg50 = (intg50 - intg30)/20.0
    intg30 = (intg30 - intg10)/20.0
    intg10 = intg10/10.0

    if var == 'no3'     :  NIT["0-10"][0] = intg10; NIT["10-30"][0] = intg30; NIT["30-50"][0] = intg50; NIT["50-100"][0] = intg100;
    if var == 'sil'     :  SIL["0-10"][0] = intg10; SIL["10-30"][0] = intg30; SIL["30-50"][0] = intg50; SIL["50-100"][0] = intg100;
    if var == 'pho'     :  PHO["0-10"][0] = intg10; PHO["10-30"][0] = intg30; PHO["30-50"][0] = intg50; PHO["50-100"][0] = intg100;
    if var == 'total_ch'     :  CHL["0-10"][0] = intg10; CHL["10-30"][0] = intg30; CHL["30-50"][0] = intg50; CHL["50-100"][0] = intg100;
    if var == 'ECO_oxy' :  OXY["0-10"][0] = intg10; OXY["10-30"][0] = intg30; OXY["30-50"][0] = intg50; OXY["50-100"][0] = intg100;


    savef = dataf + "upper_layer_climatology_"+str(date1)+"_"+str(date2)+".pckl"
    openf=open(savef,"wb")
    pickle.dump([NIT,PHO,SIL,OXY,CHL],openf)
    openf.close()

    print("")
    print("now take integrates")
    print("")
    DIA={}; FLA={}; MIC={}; MES={}; MES750={}; PP={}; CHL={}; CCL={}; DOC={}; DET={}; DOC100={}; DET100={}
    for m in range(1,13):

        arcf = dataf + "CLIM." + str(m).zfill(2) + "." + str(date1) + "." + str(date2) + ".a"
        arc  = abplot.openfile(arcf)
        print("month:",m,"file:",arcf)

        DIA[m] = abplot.integr(arc,'ECO_dia',depth)
        FLA[m] = abplot.integr(arc,'ECO_fla',depth)
        CCL[m] = abplot.integr(arc,'ECO_ccl',depth)  
        MIC[m] = abplot.integr(arc,'ECO_micr',depth)
        MES[m] = abplot.integr(arc,'ECO_meso',depth)
        PP[m]  = abplot.integr(arc,'ECO_prim',depth)*86400.
        CHL[m] = abplot.integr(arc,'total_ch',10)/10.0 # CHL is integrated to surface / and then averaged
        DOC[m] = abplot.integr(arc,'ECO_dom',depth)
        DET[m] = abplot.integr(arc,'ECO_det',depth)
        DOC100[m] = abplot.integr(arc,'ECO_dom',100)
        DET100[m] = abplot.integr(arc,'ECO_det',100) 
        MES750[m] = abplot.integr(arc,'ECO_meso',750)

    arcf = dataf + "CLIM." + str(date1) + "." + str(date2) + ".a"
    arc  = abplot.openfile(arcf)
    print("annual climatology file:", arcf)

    DIA[0] = abplot.integr(arc,'ECO_dia',depth)
    FLA[0] = abplot.integr(arc,'ECO_fla',depth)
    CCL[0] = abplot.integr(arc,'ECO_ccl',depth)
    MIC[0] = abplot.integr(arc,'ECO_micr',depth)
    MES[0] = abplot.integr(arc,'ECO_meso',depth)
    PP[0]  = abplot.integr(arc,'ECO_prim',depth)*86400.
    CHL[0] = abplot.integr(arc,'total_ch',10)/10.0 # CHL is integrated to surface / and then averaged
    DOC[0] = abplot.integr(arc,'ECO_dom',depth)
    DET[0] = abplot.integr(arc,'ECO_det',depth)
    DOC100[0] = abplot.integr(arc,'ECO_dom',100)
    DET100[0] = abplot.integr(arc,'ECO_det',100) 
    MES750[0] = abplot.integr(arc,'ECO_meso',750)

    savef = dataf + "climatology_integrates_"+str(date1)+"_"+str(date2)+".pckl"
    openf=open(savef,"wb")
    pickle.dump([DIA,FLA,CCL,MIC,MES,MES750,PP,CHL,DOC,DOC100,DET,DET100],openf)
    openf.close()

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('date1',  type=int, help='e.g 2010')
    parser.add_argument('date2',  type=int, help='e.g 2010')
    parser.add_argument('--depth',  type=float, default=200, help='depth to integrate towards surface')
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.date1,args.date2,depth=args.depth,workdir=args.workdir)
