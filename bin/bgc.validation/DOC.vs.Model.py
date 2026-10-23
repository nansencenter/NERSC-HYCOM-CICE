import pickle
#import matplotlib
#matplotlib.use('WebAgg')
import matplotlib.pyplot as plt
import abplot
import argparse
import getpass
user = getpass.getuser()
import numpy as np

# python DOC.vs.Model.py TP2 070 1990 2020
# Even though there is no Argo in 90's, the code still works

def main(region,experiment,date1,date2,workdir):

    experiment = 'expt_'+experiment[0:2]+'.'+experiment[2]
    if region == "TP0":
        region = "TP0a1.00"
    if region == "TP2":
        region = "TP2a0.10"
    if region == "TP4":
        region = "TP4a0.12"
    if region == "TP5":
        region = "TP5a0.06"
    if region == "NAT":
        region = "NATa1.00"

    plon,plat = abplot.getgrid(region=region)
    depthm = abplot.getdepth(region=region)
    
    DATAf = open("./doc/doc.pckl","rb")
    dDATE, dLAT, dLON, dPRES, dDOC = pickle.load(DATAf)

    COLOCATED = {'model':[],'obs':[],\
                 'depth':[],'lat':[],'lon':[],'date':[]}

    for p in range(len(dDATE)):
        year = dDATE[p].year
        doy = dDATE[p].timetuple().tm_yday
        is_it_true = ( np.logical_and( \
            np.logical_and( int(year) >= int(date1) , int(year) <= int(date2) ), \
             dLAT[p] > plat.min() )
        ) 
        if is_it_true:
            constructed_name = workdir + user + "/" + \
                               region + "/" + experiment + "/data/archm." + \
                               str(year) + "_" + str(doy).zfill(3) + "_12.a"  
            #print('constructing for ',constructed_name)
            ab = abplot.openfile(constructed_name)
            model = abplot.get3d(ab,'ECO_dom')
            depth = abplot.depth3d(ab)
                
            cooINDEX = abs( plat-dLAT[p] ) + abs( plon-dLON[p] )
            JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

            if not np.ma.is_masked(depthm[JJ,II]):
                try: 
                    model1d = model[:,JJ,II]
                    depth1d = depth[:,JJ,II]
                    interpolated = np.interp(dPRES[p], depth1d, model1d)

#                    index = interpolated < 0.0
#                    interpolated = np.ma.masked_where(index, interpolated)

                    if interpolated >= 0.0 and dPRES[p] <= depth1d.max():

                        COLOCATED['lat'].append( dLAT[p] )
                        COLOCATED['lon'].append( dLON[p] )
                        COLOCATED['depth'].append( dPRES[p] )
                        COLOCATED['date'].append( dDATE[p] )
                        COLOCATED['obs'].append( dDOC[p] )
                        COLOCATED['model'].append( interpolated )

                        print(dDATE[p], dLAT[p], dLON[p], dPRES[p], dDOC[p], interpolated)

                except:
                      pass


    fOUT = open(workdir + user + "/" + \
        region + "/" + experiment + "/data/DOC.Model.Colocated."+date1+"_"+date2+".pckl",'wb')
    pickle.dump(COLOCATED,fOUT)
    fOUT.close()



if __name__ == "__main__" :
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('date1',  type=str, default='1900', help='e.g 2010')
    parser.add_argument('date2',  type=str, default='2100', help='e.g 2010')
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")

    args = parser.parse_args()

    main(args.region,args.experiment,args.date1,args.date2,workdir=args.workdir)
    


