import sys,os,fnmatch
#fram = "/cluster/work/users/cagyum/Argo-GMF"
#adana = '/Users/cagyum/Desktop/Argo-GMF'
#sys.path.append(adana)
#sys.path.append(fram)
import argoinput
import pickle
import numpy as np

current_dir = os.getcwd()

dirPOC = "/nird/datapeak/NS9481K/BGCDATA/ARGO_DAC_Sprof/ARCMFC/BGC/"
dirOXY = "/nird/datapeak/NS9481K/BGCDATA/ARGO_DAC_Sprof/ARCMFC/BIO/"
dirNIT = "/nird/datapeak/NS9481K/BGCDATA/ARGO_DAC_Sprof/ARCMFC/BGC+/"

nPOCargo = np.size(sorted(fnmatch.filter(os.listdir(dirPOC), '*.nc')))
POC = {}
argocount=0
for f in sorted(fnmatch.filter(os.listdir(dirPOC), '*.nc')):
    print(" ")
    print(f)
    number = f[0:7]
    print("The code is processing POC from number: ",number)

    try:
        target = dirPOC + number + "_Sprof.nc"
        link_name = current_dir + "/GL_PR_PF_" + number + ".nc"

        if os.path.islink(link_name):
            print(f"Symlink already exists: {link_name}")
        elif os.path.exists(link_name):
            print(f"A file or directory already exists at: {link_name}")
        else:
            os.symlink(target, link_name)
            print(f"Created symlink: {link_name} -> {target}")

        pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = argoinput.readARGO(number, 'TEMP_ADJUSTED',all=True)
        poc, pocqc = argoinput.readARGO(number, 'pocK')
        abinnedt,abinnedd,abinned = argoinput.bindata(pres,dates,poc,dtop=0,dbot=2500)

        POC[argocount] = {"number": number, "lat": lat, "lon": lon, "date": dates, "raw": poc, "rawd": pres, "binned": abinned, "binnedt": abinnedt }

        os.unlink(link_name)

        argocount = argocount + 1
    except:
        pass

outf = open("./POC.argo.processed.pckl","wb")
pickle.dump(POC,outf)
outf.close()

nCHLargo = np.size(sorted(fnmatch.filter(os.listdir(dirPOC), '*.nc')))
CHL = {}
argocount=0
for f in sorted(fnmatch.filter(os.listdir(dirPOC), '*.nc')):
    print(" ")
    print(f)
    number = f[0:7]
    print("The code is processing CHL from number: ",number)

    try:
        target = dirPOC + number + "_Sprof.nc"
        link_name = current_dir + "/GL_PR_PF_" + number + ".nc"

        if os.path.islink(link_name):
            print(f"Symlink already exists: {link_name}")
        elif os.path.exists(link_name):
            print(f"A file or directory already exists at: {link_name}")
        else:
            os.symlink(target, link_name)
            print(f"Created symlink: {link_name} -> {target}")
        pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = argoinput.readARGO(number, 'TEMP_ADJUSTED',all=True)
        chl, chlqc = argoinput.readARGO(number, 'CHLA_ADJUSTED')
        abinnedt,abinnedd,abinned = argoinput.bindata(pres,dates,chl,dtop=0,dbot=250)

        CHL[argocount] = {"number": number, "lat": lat, "lon": lon, "date": dates, "raw": poc, "rawd": pres, "binned": abinned, "binnedt": abinnedt }

        os.unlink(link_name)

        argocount = argocount + 1
    except:
        pass

outf = open("./CHL.argo.processed.pckl","wb")
pickle.dump(CHL,outf)
outf.close()


nOXYargo = np.size(sorted(fnmatch.filter(os.listdir(dirOXY), '*.nc')))
OXY = {}
directories = {}
directories[0] = {'d':dirOXY}
directories[1] = {'d':dirPOC}
argocount=0
for dd in [0,1]:
 for f in sorted(fnmatch.filter(os.listdir(directories[dd]['d']), '*.nc')):
    print(" ")
    print(f)
    number = f[0:7]
    print("The code is processing OXY from number: ",number)

    try:
        target = directories[dd]['d'] + number + "_Sprof.nc"
        link_name = current_dir + "/GL_PR_PF_" + number + ".nc"

        if os.path.islink(link_name):
            print(f"Symlink already exists: {link_name}")
        elif os.path.exists(link_name):
            print(f"A file or directory already exists at: {link_name}")
        else:
            os.symlink(target, link_name)
            print(f"Created symlink: {link_name} -> {target}")
        pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = argoinput.readARGO(number, 'TEMP_ADJUSTED',all=True)
        oxy, oxyqc = argoinput.readARGO(number, 'DOXY_ADJUSTED')
        abinnedt,abinnedd,abinned = argoinput.bindata(pres,dates,oxy,dtop=0,dbot=1000)

        OXY[argocount] = {"number": number, "lat": lat, "lon": lon, "date": dates, "raw": poc, "rawd": pres, "binned": abinned, "binnedt": abinnedt }

        os.unlink(link_name)

        argocount = argocount + 1
    except:
        pass

outf = open("./OXY.argo.processed.pckl","wb")
pickle.dump(OXY,outf)
outf.close()


nNITargo = np.size(sorted(fnmatch.filter(os.listdir(dirNIT), '*.nc')))
NIT = {}
argocount=0
for f in sorted(fnmatch.filter(os.listdir(dirNIT), '*.nc')):
    print(" ")
    print(f)
    number = f[0:7]
    print("The code is processing NIT from number: ",number)

    try:
        target = dirNIT + number + "_Sprof.nc"
        link_name = current_dir + "/GL_PR_PF_" + number + ".nc"

        if os.path.islink(link_name):
            print(f"Symlink already exists: {link_name}")
        elif os.path.exists(link_name):
            print(f"A file or directory already exists at: {link_name}")
        else:
            os.symlink(target, link_name)
            print(f"Created symlink: {link_name} -> {target}")
        pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = argoinput.readARGO(number, 'TEMP_ADJUSTED',all=True)
        nit, nitqc = argoinput.readARGO(number, 'NITRATE_ADJUSTED')
        abinnedt,abinnedd,abinned = argoinput.bindata(pres,dates,nit,dtop=0,dbot=1000)

        NIT[argocount] = {"number": number, "lat": lat, "lon": lon, "date": dates, "raw": poc, "rawd": pres, "binned": abinned, "binnedt": abinnedt }

        os.unlink(link_name)

        argocount = argocount + 1
    except:
        pass

outf = open("./NIT.argo.processed.pckl","wb")
pickle.dump(NIT,outf)
outf.close()
