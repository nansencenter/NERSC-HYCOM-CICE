import pickle
#import matplotlib
#matplotlib.use('WebAgg')
import matplotlib.pyplot as plt
import abplot
import argparse
import getpass
user = getpass.getuser()
import numpy as np

# python Argo.POC.vs.Model.py main TP2 070 1990 2020
# Even though there is no Argo in 90's, the code still works

def main(region,experiment,date1,date2,workdir):

    expout = experiment
    regionout = region

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
    
    DATAf = open("/cluster/projects/nn9481k/BGC.Validation/POC.argo.processed.pckl","rb")
    DATA = pickle.load(DATAf)

    COLOCATED = {'model1d':[],'model1m':[],'argo1d':[],'argo1m':[],\
                 'depth1d':[],'lat1d':[],'lon1d':[],'date1d':[],'lat1m':[],'lon1m':[],'date1m':[]}

    for argo in range(len(DATA)):
        print('Argo: ',DATA[argo]['number'])
        argol = DATA[argo]['binned'].shape[0]
        argod = DATA[argo]['binned'].shape[1]
        for loc in range(argol):
            #print('Date: ', DATA[argo]['date'])
            year = DATA[argo]['date'][loc].year
            doy = DATA[argo]['date'][loc].timetuple().tm_yday
            #print('Argo dates ',DATA[argo]['date'][0]," to ",DATA[argo]['date'][-1])
            if np.logical_and( int(year) >= int(date1) , int(year) <= int(date2) ):
                constructed_name = workdir + user + "/" + \
                               region + "/" + experiment + "/data/archm." + \
                               str(year) + "_" + str(doy).zfill(3) + "_12.a"  
                print('constructing for Argo',DATA[argo]['number'],constructed_name)
                ab = abplot.openfile(constructed_name)
                model = abplot.get3d(ab,'poc')
                depth = abplot.depth3d(ab)
                
                cooINDEX = abs( plat-DATA[argo]['lat'][loc] ) + abs( plon-DATA[argo]['lon'][loc] )
                JJ,II = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

                if not np.ma.is_masked(depthm[JJ,II]):
                  try: 
                    model1d = model[:,JJ,II]
                    depth1d = depth[:,JJ,II]
                    interpolated = np.interp(depth1d, np.arange(argod), DATA[argo]['binned'][loc])

                    index = interpolated < 0.0
                    interpolated = np.ma.masked_where(index, interpolated)
                    model1d = np.ma.masked_where(index, model1d)
                    model1m = np.interp(np.arange(argod),depth1d,model1d)
                    argo1m = np.ma.masked_where(DATA[argo]['binned'][loc]<0.0,DATA[argo]['binned'][loc])
                    model1m = np.ma.masked_where(DATA[argo]['binned'][loc]<0.0,model1m)

                    dlim = np.arange(argod)[~DATA[argo]['binned'][loc].mask].max()
                    print('argo depth:',dlim,'| model depth:',depthm[JJ,II])
                    print(DATA[argo]['lat'][loc], DATA[argo]['lon'][loc])
                    print(plat[JJ,II],plon[JJ,II])
                    for d in range(50):
                        #if not np.ma.is_masked(interpolated[d]) : print(depth1d[d],model1d[d],interpolated[d])
                        COLOCATED['argo1d'].append( interpolated[d] ) 
                        COLOCATED['model1d'].append( model1d[d] )
                        COLOCATED['depth1d'].append( depth1d[d] )
                        COLOCATED['date1d'].append( DATA[argo]['date'][loc] )
                        COLOCATED['lat1d'].append( DATA[argo]['lat'][loc] )
                        COLOCATED['lon1d'].append( DATA[argo]['lon'][loc] )

                    for d in range(argod):
                        #if not np.ma.is_masked(model1m[d]) : print(d,model1m[d],argo1m[d])
                        COLOCATED['argo1m'].append( argo1m[d] ) 
                        COLOCATED['model1m'].append( model1m[d] )
                        COLOCATED['date1m'].append( DATA[argo]['date'][loc] )
                        COLOCATED['lat1m'].append( DATA[argo]['lat'][loc] )
                        COLOCATED['lon1m'].append( DATA[argo]['lon'][loc] )

                  except:
                      pass

                print('')

    fOUT = open(workdir + user + "/" + \
        region + "/" + experiment + "/data/POC.Model.Argo.Colocated."+date1+"_"+date2+".pckl","wb")
    pickle.dump(COLOCATED,fOUT)
    fOUT.close()

def plot(workdir):
    print('plotting')
    fIN = open(workdir + user + "/POC.Model.Argo.Colocated.pckl",'rb')
    PLOT = pickle.load(fIN)
    fIN.close()

    plt.plot(PLOT['model1d'],PLOT['depth1d'],linestyle=':',color='k')

    plt.show()

if __name__ == "__main__" :
    '''
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
    '''

    parser = argparse.ArgumentParser(description="Run different modes")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subparser for 'main'
    parser_main = subparsers.add_parser("main", help="Run main analysis")
    parser_main.add_argument("region", type=str)
    parser_main.add_argument("experiment", type=str)
    parser_main.add_argument("date1", type=str)
    parser_main.add_argument("date2", type=str)
    parser_main.add_argument("--workdir", type=str, default="/cluster/work/users/")

    # Subparser for 'plot'
    parser_plot = subparsers.add_parser("plot", help="Generate plots")
#    parser_plot.add_argument("file", type=str)
#    parser_plot.add_argument("variable", type=str)
#    parser_plot.add_argument("--output", type=str, default="plot.png")
    parser_plot.add_argument("--workdir", type=str, default="/cluster/work/users/")

    args = parser.parse_args()

    if args.command == "main":
        main(args.region, args.experiment, args.date1, args.date2, args.workdir)
    elif args.command == "plot":
        plot(args.workdir)

