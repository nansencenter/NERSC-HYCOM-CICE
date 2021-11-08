import numpy as np
import abfile.abfile as abf
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import argparse
import matplotlib.pyplot as plt
import getpass
'''
CAGLAR - Aug28,2020
how to use:

the code assumes by default you are using FRAM
thus looks for the directory '/cluster/work/users/'
If you are using another machine, use the argument:
--workdir='XXX'
The code also assumes that the region is located as follows:
i.e '/cluster/work/users/cagyum/TPxxx' - this script will retrieve your username

region options are:
TP0, TP2, TP4, TP5, NAT and NA2
the script will assume NA2 and TP4 are old hycom codes

provide the experiment as follows: 020

month as integer, e.g. January is 1

accepts only these tracers: nitrate, phosphate, silicate, volume

execute anywhere:
python hycom_plot_river.py TP0 020 nitrate 6
or
python hycom_plot_river.py TP0 020 nitrate 6 --clim="0,0.05"
'''
def main(region,experiment,tracer,month,cmap,clim,workdir):

    user = getpass.getuser()
    if region == 'TP0' or region == 'TP2' or region == 'TP5' or region == 'NAT' :
        version = 'new'
    elif region == 'TP4' or region == 'NA2':
        version = 'old'
    else :
        print('wrong or not implemented domain is provided')
        quit()

    if tracer == "nitrate" or tracer == 'phosphate' or \
        tracer == 'silicate' or tracer == "volume" :
        pass
    else :
        print('tracer name not correct')
        print('choose: nitrate, phosphate, silicate, volume')
        quit()

    units = "(mgC m$^{-2}$ s$^{-1}$), nutrient converted to C"
    if version == "new":
        if tracer == "nitrate":
            name = "ECO_no3"
        elif tracer == "phosphate":
            name = "ECO_pho"
        elif tracer == "silicate":
            name = "ECO_sil"
        else:
            name = "rivers"
            units = "(m s$^{-1}$)"

    if version == "old":
        if tracer == "nitrate":
            name = "rivnitr"
        elif tracer == "phosphate":
            name = "rivphos"
        elif tracer == "silicate":
            name = "rivsili"
        else:
            name = "rivers"
            units = "(m s$^{-1}$)"

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
    if region == "NA2":
        region = "NA2a0.80"


    abgrid = abf.ABFileGrid(workdir + user + "/" + \
        region + "/topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    jdm,idm=plon.shape

    abriver = abf.AFile(idm,jdm,workdir + user + "/" + \
        region + "/force/rivers/" + experiment + \
            "/" + name + ".a","r")
    river  = abriver.read_record(int(month)-1)
    abriver.close()

    if version == "old":
        if tracer == "nitrate":
            river = river * 12.01 * 6.625 / 14.007
        elif tracer == "phosphate":
            river = river * 12.01 * 106. / 30.97
        elif tracer == "silicate":
            river = river * 12.01 * 106. / 28.09

    abdepth = abf.ABFileBathy(workdir + user + "/" + \
        region + "/force/rivers/SCRATCH/regional.depth.b","r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")
    river = np.ma.masked_where(depthm.mask,river)

    # now plot
    fig=plt.figure(figsize=(6,5),facecolor='w')
    ax = fig.add_subplot(1,1,1)
    ax.set_position([0.01,0.02,0.865,0.96])
    cmap = plt.get_cmap(cmap)
    ax.set_facecolor('xkcd:gray')
    pmesh = plt.pcolormesh(river,cmap=cmap)
    cb=ax.figure.colorbar(pmesh)
    if clim is not None : pmesh.set_clim(clim)
    plt.text(0.015,1.05, "%s %s %s %s" %("river", tracer, "flux month:",month.zfill(2)),transform=ax.transAxes,fontsize=13)
    plt.text(0.015,1.005, units,transform=ax.transAxes,fontsize=8)

    # save figure
    fig.canvas.print_figure(workdir+user+"/"+region+"/force/rivers/"+experiment+"/"+name+".png",dpi=180)

if __name__ == "__main__" :

    class ClimParseAction(argparse.Action) :
        def __call__(self, parser, args, values, option_string=None):
           tmp = values.split(",")
           tmp = [float(elem) for elem in tmp[0:2]]
           setattr(args, self.dest, tmp)

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('region',  type=str)
    parser.add_argument('experiment',  type=str)
    parser.add_argument('tracer',  type=str)
    parser.add_argument('month', type=str)
    parser.add_argument('--cmap',  type=str,default="Spectral_r",\
        help="matplotlib colormap to use")
    parser.add_argument('--clim',  \
        action=ClimParseAction,default=None,help="range of colormap")
    parser.add_argument('--workdir', type=str, \
        default="/cluster/work/users/", \
            help="machine specific work directory (above user folder)")
    args = parser.parse_args()

    main(args.region,args.experiment,args.tracer,args.month,\
        cmap=args.cmap,clim=args.clim,workdir=args.workdir)
