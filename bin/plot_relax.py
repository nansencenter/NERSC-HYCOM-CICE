import abfile
import matplotlib.pyplot as plt
import getpass
import argparse
import numpy as np
import gridxsec
import modeltools.hycom
import os,fnmatch
import os.path
'''
CAGLAR - Sep03,2020
how to use:
python hycom_plot_relax.py

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

accepts only these tracers: nitrate, phosphate, silicate,oxygen,dic,alkalinity,temperature,salinity

execute anywhere

some options:
plots either vertical or horizontal. You need to provide one of them. See below.
if vertical, you need to provide:
--section='-70,10,50,50' # for geographical coordinates
--section='0,100,50,50' --ij # for model native indexes
  sequence: lon1,lon2,lat1,lat2

if horizontal, you need to provide layer number
--layer=1   

if you want to limit color range
--clim="0,2500"

plots temperature vertical section along geographical coordinates:
python plot_relax.py TP0 020 temperature 6 vertical --section='-70,10,50,50'

plots nitrate vertical section along model native grid indexes
python plot_relax.py TP0 020 nitrate 6 vertical --section='0,100,50,50' --ij

plots nitrate at level 1
python plot_relax.py TP0 020 nitrate 6 horizontal --layer=1
'''
def main(region,experiment,tracer,month,axis,layer,cmap,clim,workdir,\
    section,ijspace):
    user = getpass.getuser()

    if region == 'TP0' or region == 'TP2' or region == 'TP5' or region == 'NAT' :
        version = 'new'
    elif region == 'TP4' or region == 'NA2':
        version = 'old'
    else :
        print('wrong or not implemented domain is provided')
        quit()
    
    if axis == 'horizontal' or axis == 'vertical':
        pass
    else:
        print('provide arguments "horizontal" or "vertical"')
        quit()

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

    if tracer == "nitrate" or tracer == 'phosphate' or \
        tracer == 'silicate' or tracer == "oxygen" or \
           tracer == 'dic' or tracer == "alkalinity" or \
              tracer == 'temperature' or tracer == "salinity" :
        pass
    else :
        print('tracer name not correct')
        print('choose: nitrate, phosphate, silicate, \
            oxygen, dic, alkalinity, temperature, salinity')
        quit()

    units = "(mgC m$^{-2}$ s$^{-1}$), tracer converted to C"
    key="trc" # trc for tracers,  below tem and sal is set if necessary
    if version == "new":
        if tracer == "nitrate":
            name = "relax.ECO_no3"
        elif tracer == "phosphate":
            name = "relax.ECO_pho"
        elif tracer == "silicate":
            name = "relax.ECO_sil"
        elif tracer == "oxygen":
            name = "relax.ECO_oxy"
            units = "mmol m$^{-3}$"
        elif tracer == "dic":
            name = "relax.CO2_dic"
            units = "mmol m$^{-3}$"
        elif tracer == "alkalinity":
            name = "relax.CO2_alk"
            units = "mEq m$^{-3}$"
        elif tracer == "temperature":
            name = "relax_tem"
            key="tem"
            units = "$^{o}C$"
        elif tracer == "salinity":
            name = "relax_sal"
            key="sal"
            units = "psu"

    if version == "old":
        if tracer == "nitrate":
            name = "relax_nit"
        elif tracer == "phosphate":
            name = "relax_pho"
        elif tracer == "silicate":
            name = "relax_sil"
        elif tracer == "oxygen":
            name = "relax_oxy"
            units = "mmol m$^{-3}$"
        elif tracer == "dic":
            name = "relax_dic"
            units = "mmol m$^{-3}$"
        elif tracer == "alkalinity":
            name = "relax_alk"
            units = "mEq m$^{-3}$"
        elif tracer == "temperature":
            name = "relax_tem"
            key="tem"
            units = "$^{o}C$"
        elif tracer == "salinity":
            name = "relax_sal"
            key="sal"
            units = "psu"

    abgrid = abfile.ABFileGrid(workdir + user + "/" + \
        region + "/topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    jdm,idm=plon.shape    

    abdepth = abfile.ABFileBathy(workdir + user + "/" + \
        region + "/relax/" + experiment + "/SCRATCH/regional.depth.b", \
            "r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")

    abrelax = abfile.ABFileRelax(workdir + user + "/" + \
        region + "/relax/" + experiment + \
            "/" + name + ".a","r")

    # now plot
    fig=plt.figure(figsize=(6,5),facecolor='w')
    ax = fig.add_subplot(1,1,1)
    ax.set_position([0.01,0.02,0.865,0.96])
    cmap = plt.get_cmap(cmap)
    ax.set_facecolor('xkcd:gray')
    if axis == 'horizontal':
       if layer is None:
          print(" ")
          print("provide the layer number to be plotted, e.g. --layer=1") 
          print("quitting ...")
          print(" ")
          quit()
       else:
          relax  = abrelax.read_field(key,np.int(layer),np.int(month)-1)
          abrelax.close()
          pmesh = plt.pcolormesh(relax,cmap=cmap)
          cb=ax.figure.colorbar(pmesh)
          if clim is not None : pmesh.set_clim(clim)

    if axis == 'vertical':
        if section is None:
          print(" ")
          print("provide the section to be plotted")
          print("--section='lon1,lon2,lat1,lat2'")
          print("--section='-30.1,2.5,50.0,75.5'") 
          print("quitting ...")
          print(" ")
          quit()    
        else:
          lon1 = section[0]
          lon2 = section[1]
          lat1 = section[2]
          lat2 = section[3]        

        # pick up indexes
        if ijspace:
           sec = gridxsec.SectionIJSpace([lon1,lon2],[lat1,lat2],plon,plat) 
        else:
           sec = gridxsec.Section([lon1,lon2],[lat1,lat2],plon,plat)
        I,J=sec.grid_indexes
        dist=sec.distance
        slon=sec.longitude
        slat=sec.latitude                

        dpname = modeltools.hycom.layer_thickness_variable["archive"]
        # get arbitrary relaxation thicknesses
        dummy = sorted(fnmatch.filter(\
            os.listdir(workdir + user + "/" + \
                region + "/relax/" + experiment), \
                    'relax.0000_[012]*_00.a'))

        dummyarch = workdir + user + "/" + \
            region + "/relax/" + experiment + \
               "/" + dummy[np.int(month)-1]
        dummyfile = abfile.ABFileArchv(dummyarch,"r")
        kdm = max(dummyfile.fieldlevels)
        intfsec=np.zeros((kdm+1,I.size))
        datasec=np.zeros((kdm+1,I.size))
        for k in range(kdm):
            dp2d = dummyfile.read_field(dpname,k+1)
            data2d = abrelax.read_field(key,k+1,np.int(month)-1)
            dp2d=np.ma.filled(dp2d,0.)
            dp2d = dp2d/modeltools.hycom.onem
            data2d=np.ma.filled(data2d,1e30)

            intfsec[k+1,:] = intfsec[k,:] + dp2d[J,I]
            if k==0 : datasec[k,:] = data2d[J,I]
            datasec[k+1,:] = data2d[J,I]
        datasec = np.ma.masked_where(datasec>0.5*1e30,datasec)

        x = dist/1000. # km
        pmesh=ax.pcolormesh(x,-intfsec,datasec,cmap=cmap)
        cb=ax.figure.colorbar(pmesh)
        if clim is not None : pmesh.set_clim(clim)

        axm=fig.add_subplot(212)
        axm.set_position([0.85,0.85,0.15,0.15])
        axm.set_facecolor('xkcd:gray')
        pmesh2 = plt.pcolormesh(dummyfile.read_field(dpname,1)*0.,cmap=cmap)
        pltsec = plt.plot(I,J,'r',lw=2)
        plt.xticks([])
        plt.yticks([])

    plt.text(0.015,1.05, "%s %s %s" %(tracer, "relaxation month:",\
        month.zfill(2)),transform=ax.transAxes,FontSize=13)
    if axis == 'horizontal' :
       plt.text(0.015,1.005, units + " @layer=" + layer,transform=ax.transAxes,FontSize=8)
    else:
       plt.text(0.015,1.005, units,transform=ax.transAxes,FontSize=8)
    # save figure

    counter = 0
    plottemp = workdir + user + "/" + \
        region + "/relax/" + experiment + \
            "/" + name + "_" + axis + "_%s" + ".png"
    plotname = plottemp %(np.str(counter).zfill(3))

    while os.path.isfile(plotname):
          counter += 1
          plotname = plottemp %(np.str(counter).zfill(3))
          
    fig.canvas.print_figure(plotname,dpi=180)
    print(" ")
    print("figure: " + plotname)
    print(" ")

if __name__ == "__main__" :

    class ClimParseAction(argparse.Action) :
        def __call__(self, parser, args, values, option_string=None):
           tmp = values.split(",")
           tmp = [float(elem) for elem in tmp[0:2]]
           setattr(args, self.dest, tmp)
    class CoordParseAction(argparse.Action) :
        def __call__(self, parser, args, values, option_string=None):
           tmp = values.split(",")
           tmp = [float(elem) for elem in tmp[0:4]]
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
    parser.add_argument('axis',  type=str, help='horizontal or vertical')
    parser.add_argument('--layer', type=str, default=None, help="start with 1")
    parser.add_argument('--section',  \
        action=CoordParseAction,default=None,\
            help="plotting section, lon1,lon2,lat1,lat2")
    parser.add_argument('--ij', action="store_true",default=False)
    args = parser.parse_args()

    main(args.region,args.experiment,args.tracer,args.month,args.axis,\
        cmap=args.cmap,clim=args.clim,workdir=args.workdir,\
            layer=args.layer,section=args.section,ijspace=args.ij)
