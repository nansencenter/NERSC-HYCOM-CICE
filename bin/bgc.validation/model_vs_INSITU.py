import pickle
import warnings
import numpy as np
import datetime
import abfile
import abplot
import getpass
import argparse
warnings.filterwarnings('ignore')
user = getpass.getuser()

def bias(s,o):
        return np.mean(s-o)
def rmse(s,o):
        return np.sqrt(np.mean((s-o)**2))


f=open("/nird/projects/NS9481K/BGC.Validation/INSITU.OMmask.pckl","rb")
INSITU = pickle.load(f)

def main(region,experiment,date1,date2,workdir):

   subregion = "subregion"+region
   iindex = "II"+region
   jindex = "JJ"+region
   if region == "TP0":
        region = "TP0a1.00"
   if region == "TP2":
        region = "TP2a0.10"
   if region == "TP5":
        region = "TP5a0.06"
   if region == "NAT":
        region = "NATa1.00"


   experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]


   abgrid = abfile.ABFileGrid(workdir + user + "/" + \
              region + "/topo/regional.grid","r")
   plon=abgrid.read_field("plon")
   plat=abgrid.read_field("plat")    


   NO3 = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   SIL = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   OXY = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   PHO = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   CHL = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   TEM = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}
   SAL = {"lat":[], "lon":[], "model":[], "obs":[], "depth":[], "date":[], "region":[]}

   index = []
   datenow = datetime.datetime(date1,1,1)
   while datenow < datetime.datetime(date2+1,1,1):

    try:
       y = datenow.year
       m = datenow.month
       d = datenow.day
       nday = datetime.datetime( y, m, d ).strftime('%j')
       index = [ index for index, (element1, element2, element3) \
               in enumerate(zip(INSITU["year"],INSITU["month"],INSITU["day"])) \
                    if element1 == y and element2 == m and element3 == d ]

       if len(index) > 0:
          print(datetime.datetime(y,m,d))
          f = workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"archm."+str( y )+"_"+str(nday).zfill(3)+"_12.b"
          ab = abplot.openfile(f)
          no3 = abplot.get3d(ab,"no3")
          pho = abplot.get3d(ab,"pho")
          sil = abplot.get3d(ab,"sil")
          oxy = abplot.get3d(ab,"ECO_oxy")
          chl = abplot.get3d(ab,"total_ch")
          tem = abplot.get3d(ab,"temp") 
          sal = abplot.get3d(ab,"salin")
          dep = abplot.depth3d(ab) 

          no3_found = False 
          sil_found = False
          pho_found = False
          oxy_found = False
          chl_found = False
          tem_found = False
          sal_found = False
          for t in index:
              i = INSITU[iindex][ int(t) ]
              j = INSITU[jindex][ int(t) ]
              intpNO3 = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], no3[:,j,i] )
              intpPHO = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], pho[:,j,i] )
              intpSIL = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], sil[:,j,i] )
              intpOXY = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], oxy[:,j,i] )
              intpCHL = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], chl[:,j,i] )
              intpTEM = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], tem[:,j,i] )
              intpSAL = np.interp( INSITU["depth"][ int(t) ], dep[:,j,i], sal[:,j,i] )
#              print( INSITU["depth"][ int(t) ], INSITU["latitude"][ int(t) ],INSITU["longitude"][ int(t) ],\
#                     plon[j,i],plat[j,i],intpNO3)
#              print(INSITU["year"][ int(t) ],INSITU["month"][ int(t) ],INSITU["day"][ int(t) ])
              lat = INSITU["latitude"][ int(t) ]
              lon = INSITU["longitude"][ int(t) ]

              if np.logical_and( np.logical_and( intpTEM <40. , not np.isnan(INSITU["temp"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not tem_found:
                    print(' --- TEM present this day --- ')
                    tem_found = True
                 TEM["lat"].append(lat); TEM["lon"].append(lon); TEM["date"].append(datenow);
                 TEM["depth"].append(INSITU["depth"][ int(t) ]); TEM["model"].append(intpTEM);
                 TEM["region"].append(INSITU[subregion][ int(t) ]); TEM["obs"].append(INSITU["temp"][ int(t) ]);

              if np.logical_and( np.logical_and( intpSAL <50. , not np.isnan(INSITU["psal"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not sal_found:
                    print(' --- SAL present this day --- ')
                    sal_found = True
                 SAL["lat"].append(lat); SAL["lon"].append(lon); SAL["date"].append(datenow);
                 SAL["depth"].append(INSITU["depth"][ int(t) ]); SAL["model"].append(intpSAL);
                 SAL["region"].append(INSITU[subregion][ int(t) ]); SAL["obs"].append(INSITU["psal"][ int(t) ]);

              if np.logical_and( np.logical_and( intpNO3 <1000. , not np.isnan(INSITU["ntra"][ int(t) ] )), 
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not no3_found:
                    print(' --- NO3 present this day --- ')
                    no3_found = True 
                 NO3["lat"].append(lat); NO3["lon"].append(lon); NO3["date"].append(datenow);
                 NO3["depth"].append(INSITU["depth"][ int(t) ]); NO3["model"].append(intpNO3);
                 NO3["region"].append(INSITU[subregion][ int(t) ]); NO3["obs"].append(INSITU["ntra"][ int(t) ]);

              if np.logical_and( np.logical_and( intpPHO <1000. , not np.isnan(INSITU["phos"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not pho_found:
                    print(' --- PO4 present this day --- ')
                    pho_found = True
                 PHO["lat"].append(lat); PHO["lon"].append(lon); PHO["date"].append(datenow);
                 PHO["depth"].append(INSITU["depth"][ int(t) ]); PHO["model"].append(intpPHO);
                 PHO["region"].append(INSITU[subregion][ int(t) ]); PHO["obs"].append(INSITU["phos"][ int(t) ]);

              if np.logical_and( np.logical_and( intpOXY <1000. , not np.isnan(INSITU["doxy"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not oxy_found:
                    print(' --- OXY present this day --- ')
                    oxy_found = True
                 OXY["lat"].append(lat); OXY["lon"].append(lon); OXY["date"].append(datenow);
                 OXY["depth"].append(INSITU["depth"][ int(t) ]); OXY["model"].append(intpOXY);
                 OXY["region"].append(INSITU[subregion][ int(t) ]); OXY["obs"].append(INSITU["doxy"][ int(t) ]);

              if np.logical_and( np.logical_and( intpSIL <1000. , not np.isnan(INSITU["slca"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not sil_found:
                    print(' --- SIL present this day --- ')
                    sil_found = True
                 SIL["lat"].append(lat); SIL["lon"].append(lon); SIL["date"].append(datenow);
                 SIL["depth"].append(INSITU["depth"][ int(t) ]); SIL["model"].append(intpSIL);
                 SIL["region"].append(INSITU[subregion][ int(t) ]); SIL["obs"].append(INSITU["slca"][ int(t) ]);

              if np.logical_and( np.logical_and( intpCHL <1000. , not np.isnan(INSITU["cphl"][ int(t) ] )),
                    not INSITU[subregion][ int(t) ] == 'OUTSIDE') :
                 if not chl_found:
                    print(' --- CHL present this day --- ')
                    chl_found = True
                 CHL["lat"].append(lat); CHL["lon"].append(lon); CHL["date"].append(datenow);
                 CHL["depth"].append(INSITU["depth"][ int(t) ]); CHL["model"].append(intpCHL);
                 CHL["region"].append(INSITU[subregion][ int(t) ]); CHL["obs"].append(INSITU["cphl"][ int(t) ]);

       datenow = datenow + datetime.timedelta(days=1) 
    except:

       print(' ')
       print('script did not work this day:',datenow)
       datenow = datenow + datetime.timedelta(days=1)

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.NO3."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(NO3,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.PHO."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(PHO,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.SIL."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(SIL,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.OXY."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(OXY,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.CHL."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(CHL,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.TEM."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(TEM,fout)
   fout.close()

   fout = open(workdir + user + "/" + \
              region + "/" + experiment + \
              "/data/"+"colocated.SAL."+str(date1)+"."+str(date2)+".pckl","wb")
   pickle.dump(SAL,fout)
   fout.close()

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('date1',  type=int, help='e.g 2010')
    parser.add_argument('date2',  type=int, help='e.g 2010')
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.date1,args.date2,workdir=args.workdir)
