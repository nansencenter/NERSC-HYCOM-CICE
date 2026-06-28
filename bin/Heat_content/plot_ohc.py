import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
from datetime import date,timedelta

from matplotlib.dates import drange

import matplotlib.pyplot as plt


def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days+1)):
        yield start_date + timedelta(n)

def drange(start_date, end_date, ndays):
    numdays=int((end_date - start_date).days+1)
    for n in range(0,numdays,ndays):
        yield start_date + timedelta(n)

def plot_2Dfld(fld2,climP,varnam,str0,str1):
   fig=plt.figure(figsize=(5.5,5.5),facecolor='w')
   ax = fig.add_subplot(111)
   ax.set_position([0.01,0.02,0.865,0.96])
   cmap="Paired"
   #cmap="summer"
   cmap="swirly"
   cmap="viridis"
   cmap="nipy_spectral"
   cmap="cubehelix"
   cmap="rainbow"
   
   pmesh = ax.pcolormesh(fld2,cmap=cmap)
   cb=ax.figure.colorbar(pmesh)
   if climP is not None : pmesh.set_clim(climP)
   if varnam is not None:
      #plt.text(0.025,0.02,"%s"%(varnam),transform=ax.transAxes,fontsize=13)
      plt.text(0.025,0.035,"%s %s %s"%(varnam,"\n",str1),transform=ax.transAxes,fontsize=13)
      if str0 is not None:
         if len(str0)>8:
             plt.text(0.7,0.85,"%s %s %s"%(str0[0:8],"\n",str0[8:]),transform=ax.transAxes,fontsize=13)
         else:
             plt.text(0.7,0.85,"%s %s %s"%(str0[0:4],"\n",str0[4:]),transform=ax.transAxes,fontsize=13)
         plotname=varnam+'_'+str0+'.jpg'
      else:
         plotname=varnam+'.jpg'
   else:
      plotname='test.jpg'
   fig.canvas.print_figure(plotname,dpi=180)
   plt.close(fig)


def plot_3lines_anom(tdate,vart1,vart2,vart3,var1nam,var2nam,var3nam,strtext):
   import matplotlib.pyplot as plt
   import matplotlib.dates as mdates
   arr1=np.array(vart1,dtype=float)
   arr2=np.array(vart2,dtype=float)
   arr3=np.array(vart3,dtype=float)

   # pick a scale based on the largest anomaly amplitude
   amp=max(np.nanmax(np.abs(arr1)),np.nanmax(np.abs(arr2)),np.nanmax(np.abs(arr3)))
   exp=int(np.floor(np.log10(amp))) if amp>0 else 0
   scale=10**exp
   slabel="x10^%i"%exp

   if len(vart1)>20:
      fig=plt.figure(figsize=(12,5),facecolor='w')
   elif len(vart1)>10:
      fig=plt.figure(figsize=(10,5),facecolor='w')
   else:
      fig=plt.figure(figsize=(6,5),facecolor='w')

   ax=fig.add_subplot(111)
   ax.patch.set_alpha(0.8)
   ax.plot(tdate,arr1/scale,label=var1nam,linestyle='-',linewidth=1.5)
   ax.plot(tdate,arr2/scale,label=var2nam,linestyle='--',linewidth=1.5)
   ax.plot(tdate,arr3/scale,label=var3nam,linestyle=':',linewidth=1.5)
   ax.axhline(0,color='k',linewidth=0.8,linestyle='-')

   ax.xaxis.set_major_locator(mdates.YearLocator(base=3))
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   ax.set_xlabel('Year',fontsize=14)
   ax.set_ylabel('OHC anomaly (%s J/m$^2$)'%slabel,fontsize=14)
   ax.set_xlim(tdate[0],max(tdate))
   plt.legend()
   plt.grid(visible=True,color='gray',linestyle=':',linewidth=1.0)

   plotname='OHC_anom_'+strtext+'.jpg' if strtext else 'OHC_anom.jpg'
   fig.canvas.print_figure(plotname,dpi=200)
   plt.close(fig)


def plot_3lines(tdate,vart1,vart2,vart3,var1nam,var2nam,var3nam,strin,strtext):
   import matplotlib.pyplot as plt
   import matplotlib.dates as mdates
   arr1=np.array(vart1,dtype=float)
   arr2=np.array(vart2,dtype=float)
   arr3=np.array(vart3,dtype=float)

   amp=max(np.nanmax(np.abs(arr1)),np.nanmax(np.abs(arr2)),np.nanmax(np.abs(arr3)))
   exp=int(np.floor(np.log10(amp))) if amp>0 else 0
   scale=10**exp

   ave=np.array([np.nanmean(arr1),np.nanmean(arr2),np.nanmean(arr3)])/scale
   str0="Mean (x10^%i): %.2f / %.2f / %.2f (J/m2)"%(exp,ave[0],ave[1],ave[2])
   print(str0)

   if len(vart1)>20:
      fig=plt.figure(figsize=(12,5),facecolor='w')
   elif len(vart1)>10:
      fig=plt.figure(figsize=(10,5),facecolor='w')
   else:
      fig=plt.figure(figsize=(6,5),facecolor='w')

   ax=fig.add_subplot(111)
   ax.patch.set_alpha(0.8)
   ax.plot(tdate,arr1/scale,label=var1nam,linestyle='-',linewidth=1.5)
   ax.plot(tdate,arr2/scale,label=var2nam,linestyle='--',linewidth=1.5)
   ax.plot(tdate,arr3/scale,label=var3nam,linestyle=':',linewidth=1.5)

   y1,y2=plt.ylim()
   yy1=0.9*y1+0.1*y2
   yy2=0.2*y1+0.8*y2

   ax.text(tdate[5],yy1,str0,fontsize=12)
   if len(strin)>1:
      for k in range(len(strin)):
         ax.text(tdate[-8],yy2*0.98**(k),strin[k],fontsize=14)

   ax.xaxis.set_major_locator(mdates.YearLocator(base=3))
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   ax.set_xlabel('Year',fontsize=14)
   ax.set_ylabel('OHC (x10^%i J/m$^2$)'%exp,fontsize=14)
   ax.set_xlim(tdate[0],max(tdate))
   plt.legend()
   plt.grid(visible=True,color='gray',linestyle=':',linewidth=1.0)

   if strtext is not None:
      plotname='OHC_time_'+strtext+'.jpg'
   else:
      plotname='OHC_time_year.jpg'
   fig.canvas.print_figure(plotname,dpi=200)
   plt.close(fig)


def plot_2lines(tdate,vart1,vart2,var1nam,var2nam,strin,strtext):
   import matplotlib.pyplot as plt
   import matplotlib.dates as mdates
   ave=np.zeros(2)
   ave[0]=np.nanmean(np.array(vart1))
   ave[1]=np.nanmean(np.array(vart2))

   print(ave)
   tmp1=ave[0].astype(float)/10**23
   tmp2=ave[1].astype(float)/10**23
   str0="Mean (x10^23): %.2f / %.2f (J)"%(tmp1,tmp2) 
   print(str0)
   if len(vart1)>20:
      fig=plt.figure(figsize=(12,5),facecolor='w')
   elif len(vart1)>10:
      fig=plt.figure(figsize=(10,5),facecolor='w')
   else:
      fig=plt.figure(figsize=(6,5),facecolor='w')
   print(ave)
   ax = fig.add_subplot(111)
   ax.patch.set_alpha(0.8)
   cmap="bwr"
   ax.plot(tdate,vart1,label=var1nam,linestyle='-',linewidth=1.5)
   ax.plot(tdate,vart2,label=var2nam,linestyle='--',linewidth=1.5)


   y1,y2=plt.ylim()
   yy1=0.9*y1+0.1*y2
   yy2=0.2*y1+0.8*y2
   
   ax.text(tdate[5],yy1,str0,fontsize=14)
   if len(strin)>1:
      for k in range(len(strin)):
         ax.text(tdate[-8],yy2*0.98**(k),strin[k],fontsize=14)
   else:    
      ax.text(tdate[-10],ave[1],strin,fontsize=14)
   # yearly on X-axis
   ax.xaxis.set_major_locator(mdates.YearLocator(base=3))
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

   ax.set_xlabel('Year',fontsize=14)
   ax.set_xlim(tdate[0],max(tdate))
   plt.legend()
   plt.grid(visible=True,color='gray',linestyle=':',linewidth=1.0)

   if strtext is not None:
      plotname='OHC_time_'+strtext+'.jpg'
   else:
      plotname='OHC_time_year.jpg'
   fig.canvas.print_figure(plotname,dpi=200)
   plt.close(fig)


   plt.close()


def plot_myline(tdate,vart1,climP,varnam,str0):
   import matplotlib.pyplot as plt
   import matplotlib.dates as mdates
   ave=np.nanmean(np.array(vart1))

   print(vart1)
   if len(vart1)>20:
      fig=plt.figure(figsize=(12,5),facecolor='w')
   elif len(vart1)>10:
      fig=plt.figure(figsize=(10,5),facecolor='w')
   else:
      fig=plt.figure(figsize=(6,5),facecolor='w')

   ax = fig.add_subplot(111)
   #ax.set_position([0.01,0.02,0.865,0.96])
   #ax.patch.set_facecolor('gray')
   ax.patch.set_alpha(0.8)
   #cmap="rainbow"
   cmap="bwr"

   mloc=int((len(tdate)-1)/2)
   nloc=int((len(tdate)-1)/9)
   ax.plot(tdate,vart1,label=str0,linestyle='-',linewidth=1.5)

   # yearly on X-axis
   ax.xaxis.set_major_locator(mdates.YearLocator(base=3))
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

   # Optional: Add minor ticks for months for better visual reference
   #ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=6)) # Every 6 months
   fig.autofmt_xdate()
   
   if climP is not None:
      if len(climP)>0:
         ax.set_ylim(climP)
   y1,y2=plt.ylim()
   yy=0.9*y1+0.1*y2
   if vart1=='SSH' or vart1 =='ssh':
      ax.set_ylabel('Diff. '+ varnam + ' relative to E0 (cm)')
      ax.text(tdate[mloc],-.2,str1)
   elif vart1=='sst' or vart1 == 'SST':
      ax.set_ylabel('Diff. '+ varnam + " relative to E0 ({chr(176)}C)")
      ax.text(tdate[nloc],.13,str1)
   else:
      if varnam.find("BIAS")>=0:
         ax.set_ylabel(varnam + " in OHC (J)")
         str1="Mean: %.2fx10^17"%(ave/1e17)
         ax.text(tdate[mloc],yy,str1,size=13)
      elif varnam.find("RMSE")>=0:
         ax.set_ylabel(varnam + " in OHC (J)")
         str1="Mean: %.2fx10^18"%(ave/1e18)
         ax.text(tdate[mloc],yy,str1,size=13)
      else:
         ax.set_ylabel(varnam + " in OHC",fontsize=14)
         str1="Mean: %.2f"%(ave)
         ax.text(tdate[mloc],yy,str1,size=13)
      print('Default ylabel')
   #ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=np.arange(1, 13, 2)))
   ax.set_xlabel('Year',fontsize=14)
   ax.set_xlim(tdate[0],max(tdate))
   plt.legend()
   plt.grid(visible=True,color='gray',linestyle=':',linewidth=1.0)
   if varnam is not None:
      if str0 is not None:
         plotname=varnam+'_'+str0+'.jpg'
      else:
         plotname=varnam+'_year.jpg'
   else:
      plotname='timetest.jpg'
   fig.canvas.print_figure(plotname,dpi=200)
   plt.close(fig)


