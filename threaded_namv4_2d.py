#!/usr/bin/env python
#PBS -N fv3py
#PBS -l walltime=0:05:00
#PBS -l nodes=1:ppn=8
#PBS -q debug
#PBS -A fv3-cpu
#PBS -o fv3py.out
#PBS -j oe

from timeit import default_timer as timer
tic=timer()
import sys,os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap, cm, maskoceans
import multiprocessing
import numpy as np
import math
#from netCDF4 import Dataset
import numpy.ma as ma
#from netcdftime import utime
import pygrib
from datetime   import datetime,timedelta
import scipy
import ncepy
#Necessary to generate figs when not running an Xserver (e.g. via PBS)
plt.switch_backend('agg')
import pdb
import subprocess
import weather_modules as wm

######################   USER DEFINED SETTINGS    ############################################
try:
  exp=str(sys.argv[1])
  fhr=str(sys.argv[2]).zfill(2)
  startplot=datetime.strptime(str(sys.argv[3]),'%Y%m%d%H')   # start time
except:
  exp='rw_c008'
  exp='rw_019'
  exp='rw_021'
  exp='rw_022'
  exp='rw_023'
  startplot=datetime.strptime('2015103020','%Y%m%d%H')   # start time
  fhr=10                                               # number of plots to generate plus one for the start time.

outputdir='/home/Donald.E.Lippi/plotting/python/namv4-pygraphics/figs'# output directory
fhr=str(fhr).zfill(2)
fhrstr=fhr
print("exp: "+exp)
print("fhr: "+fhr)
                                                       # 1=model top; 63=model bottom.
                                                       # levs=['all'] will plot all levels.
proj="lcc"                                             # map projection
dom="SC5"                                               # domain (CONUS,NW,NWRFC,NC,NE,SC,SE,SW,MIDATL
                                                       # Great_Lakes,AK,NAK,SAK,SWAK,SEAK,PR,GUAM,HI,
                                                       # POWER,OK,LAKE_VICTORIA,AFRICA,MEDFORD))
varnames=[                                             # uncomment the desired variables below
#          'REFC',\
          'APCP',\
#          'APCP_C',\
#          'MSL',\
#          'SRHL0_1km_max',\
#          'SRHL0_3km_max',\
#          'VUSHR0-6000',\
#          'VVSHR0-6000',\
#          'CAPEsfc',\
#          'CAPE18000-0',\
#          'CAPE9000-0',\
#          'CAPE25500-0',\
##          'MXUPHL2_5km_max',\
#          'MNUPHL2_5km_min',\
#          'MAXREFC_max',\
#          'MAXREF_1km_max',\
         ]
######################   END OF USER DEFINED SETTINGS    ########################################
 #subprocess.call(["wgrib2","-v","-match",varname,nesteddata])      
 #ps=subprocess.Popen(("wgrib2","-v","-match",varname,nesteddata),stdout=subprocess.PIPE)      
 #output=subprocess.check_output(("cut","-f","4","-d",":"),stdin=ps.stdout)
 #varname=subprocess.call(["echo",blah,"|","cut","-f","4","-d",":"])
 #ps=subprocess.Popen(("echo",varname),stdout=subprocess.PIPE)
 #select_name=subprocess.check_output(("cut","-f","2","-d",":"),stdin=ps.stdout)
 #ps.wait()
import re
date=re.sub("[^0-9]","",str(startplot))
ps1=subprocess.Popen(("echo",date),stdout=subprocess.PIPE)
PDY=subprocess.check_output(("cut","-c","1-8"),stdin=ps1.stdout).strip("\n")
ps1.wait()

ps2=subprocess.Popen(("echo",date),stdout=subprocess.PIPE)
YYYY=subprocess.check_output(("cut","-c","1-4"),stdin=ps2.stdout).strip("\n")
ps2.wait()

ps3=subprocess.Popen(("echo",date),stdout=subprocess.PIPE)
MM=subprocess.check_output(("cut","-c","5-6"),stdin=ps3.stdout).strip("\n")
ps3.wait()

ps4=subprocess.Popen(("echo",date),stdout=subprocess.PIPE)
CYC=subprocess.check_output(("cut","-c","9-10"),stdin=ps4.stdout).strip("\n")
ps4.wait()
CYC=int(CYC)
cycstr=str(CYC)
if( exp == 'obs'):
   if(CYC>= 00 and CYC<=05): cyc_dir=06
   if(CYC>= 06 and CYC<=11): cyc_dir=12
   if(CYC>= 12 and CYC<=17): cyc_dir=18
   if(CYC>= 18 and CYC<=23): cyc_dir=00
   datadir="/home/Donald.E.Lippi/imagemagic/data/gefs."+PDY+"/"+str(cyc_dir)+"/ccpa/"
   nesteddata=str(os.path.join(datadir,'ccpa_conus_0.125d_t'+cycstr+'z_03h'))
else:
   datadir="/scratch4/NCEPDEV/stmp4/Donald.E.Lippi/pyplot_work_"+exp+"."+PDY
   if( not os.path.exists(datadir)): os.makedirs(datadir)
   nesteddata=str(os.path.join(datadir,'namrr.t'+cycstr+'z.conusnest.hiresf'+fhrstr+'.tm00.grib2'))     # name of file

print(nesteddata)
exit()

# Download from tape if file doesn't exist locally.
if( not os.path.isfile(nesteddata) and exp != 'obs'):
  HPSSDIR="/NCEPDEV/emc-meso/5year/Donald.E.Lippi/nw"+exp+"/rh"+YYYY+"/"+YYYY+MM+"/"+PDY+"/meso2_noscrub_Donald.E.Lippi_com_namrr_"+exp+"_namrr."+PDY+cycstr+".conusnest.tar" 
  HPSSFILE="namrr.t"+cycstr+"z.conusnest.hiresf"+fhr+".tm00.grib2"
  subprocess.call(["htar","-xvf",HPSSDIR,"./"+HPSSFILE])
  subprocess.call(["mv",HPSSFILE,datadir])

# Create the basemap
# create figure and axes instances
fig = plt.figure(figsize=(11,11))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

# Setup map corners for plotting.  This will give us CONUS
llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,res=ncepy.corners_res('SC',proj=proj)
if(dom == 'SC4'): llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-104.0,28.0,-92.0,35.0
if(dom == 'SC5'): llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat=-104.0,26.0,-92.0,33.0
m = Basemap(llcrnrlon=llcrnrlon,   llcrnrlat=llcrnrlat,
               urcrnrlon=urcrnrlon,  urcrnrlat=urcrnrlat,
               projection=proj, lat_0=35.4,lon_0=-97.6,
               resolution=res)

# Map background stuff to make things look nice
parallels = np.arange(-26.,90,4.)
meridians = np.arange(0.,360.,4.)
#m.drawmapboundary(fill_color='#7777ff')
#m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder = 0)
m.drawcoastlines(linewidth=1.25)
m.drawstates(linewidth=1.25)
m.drawcountries(linewidth=1.25)
m.drawparallels(parallels,labels=[1,0,0,1])
m.drawmeridians(meridians,labels=[1,0,0,1])
m.drawcounties(linewidth=0.2, color='k')

def mkplot(varname):
    print("mkplot - "+str(multiprocessing.current_process()))
    grbs = pygrib.open(nesteddata)


    if(False):
       grb_message_filename="./grb_messages.txt"
       grb_message_file=open(grb_message_filename,"w")
       lines_of_text=[]
       grbs.seek(0)
       for grb_message in grbs:
           lines_of_text.append(str(grb_message)+"\n")
       grb_message_file.writelines(lines_of_text)
       exit()

    global lons,lats
    lats,lons=grbs[1].latlons()
    cyctime=grbs[1].dataTime #Cycle (e.g. 1200 UTC)
    grbtime=str(cyctime).zfill(4)
    date=grbs[1].dataDate #PDY
    idate=repr(date)+grbtime #CDATE
    vdate=ncepy.ndate(idate,int(fhr))
    valpdy=vdate[0:8]
    valcyc=vdate[8:11]

    

#    if fhr==0: cyctime=cytime+'00'
#      #Pad with a zero and convert to a string
#    if cyctime < 1000:
#        grbtime='0'+repr(cyctime)
#    else:
#        grbtime=repr(cyctime)

    date=grbs[1].dataDate    #PDY

    if(varname=='REFC'):            grb=grbs.select(name="Maximum/Composite radar reflectivity")[0].values
    if(varname=='APCP'):            grb=grbs.select(name="Total Precipitation")[0].values
    if(varname=='APCP_C'):            grb=grbs.select(name="Total Precipitation")[0].values-grbs.select(name="Large scale precipitation")[0].values
    if(varname=='SRHL0_1km_max'):   grb=grbs.select(stepType="max",name="Storm relative helicity",typeOfLevel="heightAboveGroundLayer",topLevel=1000,bottomLevel=0)[0].values 
    if(varname=='SRHL0_3km_max'):   grb=grbs.select(stepType="max",name="Storm relative helicity",typeOfLevel="heightAboveGroundLayer",topLevel=3000,bottomLevel=0)[0].values
    #if(varname=='MXUPHL2_5km_max'): grb=grbs.select(name="Updraft Helicity",typeOfLevel="heightAboveGroundLayer",topLevel=5000,bottomLevel=2000)[0].values 
    if(varname=='MXUPHL2_5km_max'): grb=grbs.select(stepType='max',parameterName="199",typeOfLevel="heightAboveGroundLayer",topLevel=5000,bottomLevel=2000)[0].values 
    if(varname=='VUSHR0-6000'):     grb=grbs.select(name="Vertical u-component shear")[0].values 
    if(varname=='VVSHR0-6000'):     grb=grbs.select(name="Vertical v-component shear")[0].values 
    if(varname=='CAPEsfc'):         grb=grbs.select(name="Convective available potential energy",typeOfLevel="surface")[0].values 
    if(varname=='CAPE18000-0'):     grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="18000-0")[0].values 
    if(varname=='CAPE9000-0'):      grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="9000-0")[0].values 
    if(varname=='CAPE25500-0'):     grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="25500-0")[0].values 
    if(varname=='MSL'):
       grb=scipy.ndimage.gaussian_filter(grbs.select(name="Pressure reduced to MSL")[0].values*0.01,2)
       ug  = grbs.select(stepType='instant',name='10 metre U wind component',typeOfLevel='heightAboveGround',level=10)[0].values
       vg  = grbs.select(stepType='instant',name='10 metre V wind component',typeOfLevel='heightAboveGround',level=10)[0].values
       u10=ncepy.ms2kts(ug)
       v10=ncepy.ms2kts(vg)
       #t2m = ncepy.Kelvin2F(grbs.select(stepType='instant',name='2 metre temperature',typeOfLevel='heightAboveGround',level=2)[0].values)
       #t850 = grbs.select(name='Temperature',typeOfLevel='isobaricInhPa',level=850)[0].values
       #td2m=ncepy.Kelvin2F(grbs.select(stepType='instant',name='2 metre dewpoint temperature',typeOfLevel='heightAboveGround',level=2)[0].values)
       #td2mdepression=t2m-td2m
       #td2mdepression=scipy.ndimage.gaussian_filter(td2mdepression,3)
       dbz=grbs.select(name="Maximum/Composite radar reflectivity")[0].values
       #td850=grbs.select(name='Dew point temperature',typeOfLevel='isobaricInhPa',level=850)[0].values
       #td1000=grbs.select(name='Dew point temperature',typeOfLevel='isobaricInhPa',level=1000)[0].values
       #print(np.shape(t850))
       try:
          dx=grbs[1]['DxInMetres']/1000.
          dy=grbs[1]['DyInMetres']/1000.
       except:
          dx=grbs[1]['DiInMetres']/1000.
          dy=grbs[1]['DiInMetres']/1000.
       skip=15 #int(150/dx)+2
       scale=20.
       #datadx,datady=np.gradient(t850,dx,dy)
       #tadvect=-u*datadx -v*datady
       
     
    keep_ax_lst = ax.get_children()[:]
    print(np.max(grb))
  
    dispatcher=plot_Dictionary()
    experiment_dispatcher=exp_Dictionary()
    # Clear off old plottables but keep all the map info
    ncepy.clear_plotables(ax,keep_ax_lst,fig)
    print(varname)
    var_n=grb

    experiment=experiment_dispatcher[exp]
    function=dispatcher[varname]
    var_n,clevs,cm,units,longname=function(var_n)
   
       

    if(varname=='MSL'):
       # mslp
       if(True):
        cs = m.contour(lons,lats,var_n,clevs,colors='black',linewidths=1.75,latlon=True) 
        clab = plt.clabel(cs, inline=1,colors='k',fontsize=10,fmt='%.0f')
        highslowswindow=200
        ncepy.plt_highs_and_lows(m,var_n,lons,lats,mode='reflect',window=highslowswindow)

       # dbz
       if(True):
        clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
        cmap=ncepy.mrms_radarmap()
        cs0 = m.contourf(lons,lats,dbz,clevs,cmap=cmap,latlon=True,extend='max')

       # temperature
       if(False):
        clevs = np.arange(-36.,104.,4)
        cs1 = m.contourf(lons,lats,t2m,clevs,cmap=ncepy.ncl_t2m(),latlon=True,extend='both')
        cs1.set_clim(-36, 100.)
        cbar = m.colorbar(cs1,location='bottom',pad="5%",ticks=clevs,format='%.0f')
        cbar.ax.tick_params(labelsize=8.5)
        cbar.set_label('F')

       #temperature advection 850
       if(False):
        cmap=plt.get_cmap(name='coolwarm')
        cs5 = m.contourf(lons,lats,tadvect,cmap=cmap,latlon=True,extend='both')
        cbar = m.colorbar(cs5,location='bottom',pad="5%")
        cbar.ax.tick_params(labelsize=8.5)
        cbar.set_label('s-1')
     
       #moisture advection 850
       if(False):
        lev=925
        u  = grbs.select(stepType='instant',name='U component of wind',typeOfLevel='isobaricInhPa',level=lev)[0].values
        v  = grbs.select(stepType='instant',name='V component of wind',typeOfLevel='isobaricInhPa',level=lev)[0].values
        tdlev=grbs.select(name='Dew point temperature',typeOfLevel='isobaricInhPa',level=lev)[0].values
        qvap=wm.td_to_mixrat(tdlev,lev*100.)
        datadx,datady=np.gradient(qvap,dx,dy)
        qadvect=-u*datadx -v*datady
        cmap=plt.get_cmap(name='coolwarm')
        cs5 = m.contourf(lons,lats,qadvect,cmap=cmap,latlon=True,extend='both')
        cbar = m.colorbar(cs5,location='bottom',pad="5%")
        cbar.ax.tick_params(labelsize=8.5)
        cbar.set_label('s-1')
     
       # moisture
       if(False):
        clevs = np.arange(1.,3.0,1.)
        cmap=plt.get_cmap(name='summer_r')
        td2mdepression=ma.masked_where(td2mdepression > 2.5,td2mdepression)
        #cs1 =  m.contourf(lons,lats,td2mdepression,clevs,cmap=cmap,latlon=True,linewidths=1.15,extend='min')
        cs1 =   m.contour(lons,lats,td2mdepression,clevs,cmap=cmap,latlon=True,linewidths=1.15,extend='both')
        clab1=plt.clabel(cs1,inline=1,fontsize=8,fmt='%.0f')
        #cbar = m.colorbar(cs1,location='bottom',pad="5%",ticks=clevs,format='%.0f')
        #cbar.ax.tick_params(labelsize=8.5)
        #cbar.set_label('F')

       # instability (cape) 
       if(False):
        cape=grbs.select(name="Convective available potential energy",typeOfLevel="surface")[0].values
        clevs=[500,1000,1500,2000,2500,3000] #CAPE 
        cmap=plt.get_cmap(name='Spectral_r')
        cape=ma.masked_where(cape < 500,cape)
        cs2 = m.contourf(lons,lats,cape,clevs,cmap=cmap,latlon=True,extend='max')
        cs2.set_clim(500, 3000.)
        #clab2=plt.clabel(cs2,inline=1,fontsize=8,fmt='%.0f')
        cbar = m.colorbar(cs2,location='bottom',pad="5%",ticks=clevs,format='%.0f')
        cbar.ax.tick_params(labelsize=8.5)
        cbar.set_label('J/kg')

       # lift

       # shear
       if(False):
        grbu=grbs.select(name="Vertical u-component shear")[0].values
        grbv=grbs.select(name="Vertical v-component shear")[0].values
        shear=np.sqrt(grbu**2 + grbv**2)
        #clevs=np.arange(-20,20.5,2)
        #cmap=plt.get_cmap(name='YlGn_r')
        #cs2 = m.contourf(lons,lats,shear,clevs,cmap=cmap,latlon=True,extend='both')
        #m.barbs(x[::skip,::skip],y[::skip,::skip],grbu[::skip,::skip],grbv[::skip,::skip],length=5,rounding=True,pivot='middle',fill_empty=True,sizes=dict(emptybarb=0.05))
        m.barbs(lons[::skip,::skip],lats[::skip,::skip],grbu[::skip,::skip],grbv[::skip,::skip],latlon=True,length=3.75,sizes={'spacing':0.2},pivot='middle')

       #low level convergence
       if(False):
        conv10m=wm.h_convergence(ug,vg,dx,dy)
        cmap=plt.get_cmap(name='coolwarm')
        lev=round(max(abs(conv10m.max()),abs(conv10m.min())))
        clevs=np.arange(-1*lev*0.15,lev*0.15+0.05,0.05)
        cs4 = m.contourf(lons,lats,conv10m,clevs,cmap=cmap,latlon=True,extend='both')
        cbar = m.colorbar(cs4,location='bottom',pad="5%")
        cbar.ax.tick_params(labelsize=8.5)
        cbar.set_label('s-1')


       #Supercell Composite Parameter - SCP = (muCAPE / 1000 J kg-1) * (ESRH / 50 m2 s-2) * (EBWD / 20 m s-1)
       if(False):
        esrh=grbs.select(stepType="max",name="Storm relative helicity",typeOfLevel="heightAboveGroundLayer",topLevel=3000,bottomLevel=0)[0].values
        ebwd=np.sqrt(grbu**2 + grbv**2)
        SCP=(cape/1000.)*(esrh/50.)*(ebwd/20.)
        clevs=np.arange(1,4.,1.)
        cs3 = m.contourf(lons,lats,SCP,clevs,cmap=cmap,latlon=True,extend='both')

       # sfc winds
       if(True):
        m.barbs(lons[::skip,::skip],lats[::skip,::skip],u10[::skip,::skip],v10[::skip,::skip],latlon=True,length=3.75,sizes={'spacing':0.2},pivot='middle')







    if(varname=='REFC'):    cs = m.contourf(lons,lats,var_n,clevs,cmap=cm,latlon=True,extend='max')
    #norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)
    #cs = m.pcolormesh(lons,lats,var_n,cmap=cm,latlon=True,norm=norm,vmin=clevs[0],vmax=clevs[-1])
    if(varname=='CAPEsfc'):
      grbu=grbs.select(name="Vertical u-component shear")[0].values 
      grbv=grbs.select(name="Vertical v-component shear")[0].values 
      grb=np.sqrt(grbu**2 + grbv**2)
      x,y=m(lons,lats); q=20
      m.barbs(x[::q,::q],y[::q,::q],grbu[::q,::q],grbv[::q,::q],length=5,rounding=True,pivot='middle',fill_empty=True,sizes=dict(emptybarb=0.05))
      grb=grbs.select(name="Convective available potential energy",typeOfLevel="surface")[0].values

    if(not varname=="MSL"):
       norm = c.BoundaryNorm(clevs, cm.N)
       cs = m.contourf(lons,lats,var_n,clevs,cmap=cm,norm=norm,latlon=True,extend='max')
       cbar = m.colorbar(cs,location='bottom',pad="5%",extend="both",ticks=clevs)
       cbar.ax.tick_params(labelsize=8.5)
       cbar.set_label("["+str(units)+"]")
    plt.title("%s CONUSNEST %s \n%s %sZ cycle Fhr %s Valid %s %sZ" % (experiment,longname,date,grbtime,fhr,valpdy,valcyc))
    output=outputdir+"/"+str(date)+str(CYC)
    if not os.path.exists(output): os.makedirs(output)
    plt.savefig(output+'/'+"%s_%s_%s_f%s_CONUSNEST.png" % (experiment,varname,dom,fhr),dpi=125,bbox_inches='tight')



############### useful functions ###########################################
def roundTime(dt=None, roundTo=60):
   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt.replace(tzinfo=None) - dt.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def gemplot(clist):
    gemlist=ncepy.gem_color_list()
    colors=[gemlist[i] for i in clist]
    cm = matplotlib.colors.ListedColormap(colors)
    return cm

############### plot_ functions ###########################################

#-------------- tracer   -------------------------------------------------- 
def plot_sphum(var_n):
    """specific humidity [kg/kg]"""
    longname="specific humidity"; units="g/kg" #units="kg/kg *10-3"
    var_n=var_n*1000
    clevs=np.arange(0.,32.5,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_liq_wat(var_n):
    """liquid water content [g/kg]"""
    exit("Need to fix plot_liq_wat")
    longname="liq_wat"; units="g/kg" #units="kg/kg *10-3"
    clevs=np.arange(0.,32.5,1)
    clevs=[0.00,0.03,0.05,0.25,0.25,1.,3.,4.,5.]
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_o3mr(var_n):
    """ozone mixing ratio [kg/kg]"""
    exit("Need to fix plot_o3mr")
    longname="Ozone mixing ratio"; units="g/kg" #units="kg/kg *10-3"
    var_n=var_n*1000
    clevs=np.arange(0.,32.5,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

#-------------- core    -------------------------------------------------- 
def plot_u(var_n):
    """zonal wind [m/s]"""
    longname="zonal wind"; units="m/s"
    clevs=np.arange(-20,20.5,2)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_v(var_n):
    """meridional wind [m/s]"""
    longname="meridional wind"; units="m/s"
    clevs=np.arange(-20,20.5,2)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_W(var_n):
    """vertical velocity [m/s]"""
    longname="vertical velocity"; units="cm/s"
    if(units=="cm/s"): var_n=var_n*10
    clevs=np.arange(-20,20.5,2)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_DZ(var_n):
    """DZ [??]"""
    exit("Need to fix plot_DZ")
    longname="DZ"; units="??"
    clevs= np.arange(0,10,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_T(var_n):
    """temperature [K]"""
    longname="temperature"; units="F"
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_delp(var_n):
    """delp [??]"""
    exit("Need to fix plot_delp")
    longname="??"; units="??"
    clevs= np.arange(0,10,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_phis(var_n):
    """incremental pressure on each level [Pa]"""
    exit("Need to fix plot_phis")
    longname="incremental pressure on each level"; units="hPa"
    if(units=='hPa'): var_n=var_n*0.01
    clevs=np.arange(0,10,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)


#-------------- nggps2d -------------------------------------------------- 
def plot_ALBDOsfc(var_n):
    """surface albedo (%)"""
    longname="surface albedo"; units="%"
    clevs=np.arange(0.,100.5,5.)
    cm=ncepy.tcamt()
    return(var_n,clevs,cm,units,longname)

def plot_CPRATsfc(var_n):
    """surface convective precipitation rate [kg/m**2/s]"""
    longname="surface convective precipitation rate"; units="in."
    if(units=="in."): 
       var_n = var_n*3*3600/25.4  # inches
       clevs=[0,0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    elif(units=="mm"):
       var_n = var_n*3*3600      # mm
       clevs= [0,0.1,2,5,10,15,20,25,35,50,75,100,125,150,175]  #mm
    clist=[0,23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)


def plot_PRATEsfc(var_n):
    """surface precipitation rate [kg/m**2/s]"""
    longname="surface precipitation rate"; units="in."
    if(units=="in."):
       var_n = var_n*3*3600/25.4  # inches
       clevs=[0,0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    elif(units=="mm"):
       var_n = var_n*3*3600      # mm
       clevs= [0,0.1,2,5,10,15,20,25,35,50,75,100,125,150,175]  #mm
    clist=[0,23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_DLWRFsfc(var_n):
    """surface downward longwave flux [W/m**2]"""
    longname="surface downward longwave flux"; units="W/m**2"
    clevs=np.arange(0,525,25)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_ULWRFsfc(var_n):
    """surface upward longwave flux [W/m**2]"""
    longname="surface upward longwave flux"; units="W/m**2"
    clevs=np.arange(0,525,25)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_DSWRFsfc(var_n):
    """surface downward shortwave flux [W/m**2]"""
    longname="surface downward shortwave flux"; units="W/m**2"
    clevs=np.arange(0,1050,50)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_USWRFsfc(var_n):
    """surface upward shortwave flux [W/m**2]"""
    longname="surface upward shortwave flux"; units="W/m**2"
    clevs=np.arange(0,1050,50)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_DSWRFtoa(var_n):
    """top of atmos downward shortwave flux [W/m**2]"""
    longname="top of atmos downward shortwave flux"; units="W/m**2"
    clevs=np.arange(0,1050,50)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_USWRFtoa(var_n):
    """top of atmos upward shortwave flux [W/m**2]"""
    longname="top of atmos upward shortwave flux"; units="W/m**2"
    clevs=np.arange(0,1050,50)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_ULWRFtoa(var_n):
    """top of atmos upward longwave flux [W/m**2]"""
    longname="top of atmos upward longwave flux"; units="W/m**2"
    clevs=np.arange(0,525,25)
    #clevs=np.arange(0,15000.5,1000)
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname)

def plot_GFLUXsfc(var_n):
    """surface ground heat flux [W/m**2]"""
    longname="surface ground heat flux"; units="W/m**2"
    clevs= [-300.,-200.,-100.,-75.,-50.,-25.0,-10.0,0.,10.0,25.,50.,75.,100.,200.,300.]
    cm=ncepy.ncl_grnd_hflux()
    return(var_n,clevs,cm,units,longname)

def plot_HGTsfc(var_n):
    """surface geopotential height [gpm]"""
    longname="surface geopotential height"; units="gpm"
    clevs=[0,250.,500.,750.,1000.,1500.,2000.,3000.,4000.,5000.,7500.,10000.,15000.,20000.,25000.,30000.,30000.5]
    cm=plt.get_cmap(name='jet')
    return(var_n,clevs,cm,units,longname)

def plot_HPBLsfc(var_n):
    """surface planetary boundary layer height [m]"""
    longname="surface planetary boundary layer height"; units="m"
    clevs=[0,50.,100.,150.,200.,250.,500.,750.,1000.,1500.,2000.,3000.,4000.,5000.,7500.]
    cm=plt.get_cmap(name='jet')
    return(var_n,clevs,cm,units,longname)

def plot_ICECsfc(var_n):
    """surface ice concentration (ice=1; no ice=0) [fraction]"""
    longname="surface ice concentration"; units="(ice=1; no ice=0)"
    clevs=[0,0.5,1]
    clist=[23,27]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SLMSKsfc(var_n):
    """sea-land-ice mask (0-sea, 1-land, 2-ice)"""
    longname="sea-land-ice mask"; units="0-sea, 1-land, 2-ice"
    clevs=[0,1,2,2.01]
    clist=[24,18,27]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 

def plot_LHTFLsfc(var_n):
    """surface latent heat flux [W/m**2]"""
    longname="surface latent heat flux"; units="W/m**2"
    clevs=[-300.,-200.,-100.,-75.,-50.,-25.0,-10.0,-5.,5.,10.0,25.,50.,75.,100.,200.,300]
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_SHTFLsfc(var_n):
    """surface sensible heat flux [W/m**2]"""
    longname="surface sensible heat flux"; units="W/m**2"
    clevs=[-300.,-200.,-100.,-75.,-50.,-25.0,-10.0,0.,10.0,25.,50.,75.,100.,200.,300.]
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_PRESsfc(var_n): # done
    """surface pressure [Pa]"""
    longname="surface pressure"; units="hPa"
    if(units=='hPa'): var_n=var_n*0.01
    var_n=scipy.ndimage.gaussian_filter(var_n, 2) # first pass
    var_n=scipy.ndimage.gaussian_filter(var_n, 2) # second pass
    clevs=np.arange(950.,1050.,4.)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_PWATclm(var_n):
    """atmos column precipitable water [kg/m**2]"""
    longname="atmos column precipitable water"; units="mm"
    clevs= [4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76]
    clist = [0,16,18,30,28,27,25,4,23,3,21,8,5,19,17,31,12,2,7,14]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SOILM(var_n):
    """total column soil moisture content [kg/m**2]"""
    longname="total column soil moisture content"; units="in."
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    if(units=="in."): var_n = var_n/25.4 #inches
    clevs= [4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76]
    clist = [0,16,18,30,28,27,25,4,23,3,21,8,5,19,17,31,12,2,7,14]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SOILW1(var_n):
    """volumetric soil moisture 0-10cm [fraction]"""
    longname="volumetric soil moisture 0-10cm"; units="fraction"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    clevs=[0,.10,.20,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1]
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SOILW2(var_n):
    """volumetric soil moisture 10-40cm [fraction]"""
    longname="volumetric soil moisture 10-40cm"; units="fraction"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    clevs=[0,.10,.20,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1]
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SOILW3(var_n):
    """volumetric soil moisture 40-100cm [fraction]"""
    longname="volumetric soil moisture 40-100cm"; units="fraction"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    clevs=[0,.10,.20,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1]
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SOILW4(var_n):
    """volumetric soil moisture 100-200cm [fraction]"""
    longname="volumetric soil moisture 100-200cm"; units="fraction"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    clevs=[0,.10,.20,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1]
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SPFH2m(var_n):
    """2m specific humidity [kg/kg]"""
    longname="2m specific humidity"; units="g/kg" #units="kg/kg *10-3"
    var_n=var_n*1000
    clevs=np.arange(0.,32.5,1)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_SOILT1(var_n): #done
    """soil temperature 0-10cm [K]"""
    longname="soil temperature 0-10cm"; units="F"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_SOILT2(var_n): #done
    """soil temperature 10-40cm [K]"""
    longname="soil temperature 10-40cm"; units="F"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_SOILT3(var_n): #done
    """soil temperature 40-100cm [K]"""
    longname="soil temperature 40-100cm"; units="F"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_SOILT4(var_n): #done
    """soil temperature 100-200cm [K]"""
    longname="soil temperature 100-200cm"; units="F"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_TMP2m(var_n):
    """2m temperature [K]"""
    longname="2m temperature"; units="F"
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_TMPsfc(var_n):
    """surface temperature [K]"""
    longname="surface temperature"; units="F"
    if(units=="F"): var_n=ncepy.Kelvin2F(var_n) # [F]
    clevs= np.arange(-36.,104.,4)
    cm=ncepy.ncl_t2m()
    return(var_n,clevs,cm,units,longname)

def plot_UGWDsfc(var_n): 
    """surface zonal gravity wave stress [N/m**2]"""
    longname="surface zonal gravity wave stress"; units="N/m**2"
    clevs=[-20,-10,-5,-2.5,-1,-0.05,-0.001,0.001,0.05,1,2.5,5,10,20]
    clevs=[-5,-2.5,-1,-0.05,-0.01,0.01,0.05,1,2.5,5]
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_VGWDsfc(var_n): 
    """surface meridional gravity wave stress [N/m**2]"""
    longname="surface meridional gravity wave stress"; units="N/m**2"
    clevs=[-5,-2.5,-1,-0.05,-0.01,0.01,0.05,1,2.5,5]
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_UFLXsfc(var_n): 
    """surface zonal momentum flux [N/m**2]"""
    longname="surface zonal momentum flux"; units="N/m**2"
    clevs=np.arange(-1,1.05,0.1)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_VFLXsfc(var_n): 
    """surface meridional momentum flux [N/m**2]"""
    longname="surface meridional momentum flux"; units="N/m**2"
    clevs=np.arange(-1,1.05,0.1)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_UGRD10m(var_n): 
    """10 meter u wind [m/s]"""
    longname="10 meter u wind"; units="m/s"
    clevs=np.arange(-20,20.5,2)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_VGRD10m(var_n): 
    """10 meter v wind [m/s]"""
    longname="10 meter v wind"; units="m/s"
    clevs=np.arange(-20,20.5,2)
    cm=plt.get_cmap(name='RdBu_r')
    return(var_n,clevs,cm,units,longname)

def plot_WEASDsfc(var_n): 
    """surface snow water equivalent [kg/m**2]"""
    longname="surface snow water equivalent"; units="in."
    if(units=="in."):
       var_n = var_n/25.4  # inches
       clevs=[0,0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    elif(units=="mm"):
       clevs= [0,0.1,2,5,10,15,20,25,35,50,75,100,125,150,175]  #mm
    clist=[0,23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_SNODsfc(var_n): 
    """surface snow depth [m]"""
    longname="surface snow depth"; units="in."
    if(units=="in."):
       var_n = var_n/0.0254  # inches
       clevs=[0,0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    elif(units=="mm"):
       var_n = var_n/1000.      # mm
       clevs= [0,0.1,2,5,10,15,20,25,35,50,75,100,125,150,175]  #mm
    clist=[0,23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_ZORLsfc(var_n):
    """surface roughness [m]"""
    longname="surface roughness"; units="m"
    clevs=np.arange(0,3.1,0.1)
    cm=plt.get_cmap(name='jet') 
    return(var_n,clevs,cm,units,longname)

def plot_VFRACsfc(var_n):
    """vegetation fraction"""
    longname="vegetation fraction"; units="fraction"
    clevs=np.arange(0.,100.5,5.)
    cm=ncepy.tcamt()
    return(var_n,clevs,cm,units,longname)

def plot_F10Msfc(var_n): 
    """10-meter wind speed divided by lowest model wind speed"""
    longname="10-meter wind speed divided by lowest model wind speed"; units="none"
    clevs=np.arange(0,2.05,.1)
    cm=plt.get_cmap(name='jet') 
    return(var_n,clevs,cm,units,longname)

def plot_VTYPEsfc(var_n): 
    """vegetation type in integer 1-13"""
    longname="vegetation type in integer 1-13"; units="1-13"
    clevs=np.arange(1,15,1)
    cm=plt.get_cmap(name='jet')
    return(var_n,clevs,cm,units,longname)

def plot_STYPEsfc(var_n): 
    """soil type in integer 1-9"""
    longname="soil type"; units="1-9"
    var_n=maskoceans(lons, lats, var_n, inlands=True, resolution=res)
    clevs=np.arange(1,11,1)
    cm=plt.get_cmap(name='jet')
    return(var_n,clevs,cm,units,longname)

def plot_TCDCclm(var_n):
    """atmos column total cloud cover [%]"""
    longname="atmos column total cloud cover"; units="%"
    clevs=[0,10,20,30,40,50,60,70,75,80,85,90,95,100] # percent
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_TCDChcl(var_n):
    """high cloud level total cloud cover [%]"""
    longname="high cloud level total cloud cover"; units="%"
    clevs=[0,10,20,30,40,50,60,70,75,80,85,90,95,100] #percent
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_TCDCmcl(var_n): 
    """mid cloud level total cloud cover [%]"""
    longname="mid cloud level total cloud cover"; units="%"
    clevs=[0,10,20,30,40,50,60,70,75,80,85,90,95,100] #percent
    clist=[0,30,29,27,24,4,23,3,5,19,17,7,2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname)

def plot_TCDClcl(var_n): 
    """low cloud level total cloud cover [%]"""
    longname="low cloud level total cloud cover"; units="%"
    clevs=[0,10,20,30,40,50,60,70,75,80,85,90,95,100] # percent
    clist=[ 0,30,29,27,24, 4,23, 3, 5,19,17, 7, 2]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 

def plot_REFC(var_n):
    """Stoelinga simulated maximum (composite) reflectivity [dbz]"""
    longname="Maximum Composite Reflectivity"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 


def plot_APCP(var_n):
    """Total Precipitation [in.]"""
    longname="Total Precipitation"; units="in."
    if(units=="in."): var_n = var_n/25.5  #mm to in.
    clevs=[0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    clist=[23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 

def plot_APCP_C(var_n):
    """Total Precipitation [in.]"""
    longname="Total Precipitation"; units="in."
    if(units=="in."): var_n = var_n/25.5  #mm to in.
    clevs=[0.01,0.05,0.1,0.25,0.5,0.75,1.,1.5,2.,3.,4.,5.,6.,7.] #inches
    clist=[23,22,21,20,19,10,17,16,15,14,29,28,24,25]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 


def plot_MSL(var_n):
    """Pressure reduced to MSL [hPa]"""
    longname="Mean Sea Level Pressure"; units="hPa"
    clevs = np.arange(900,1100.,4.)
    cm=""
    return(var_n,clevs,cm,units,longname) 

def plot_CAPEsfc(var_n):
    """Surface Based Convective Available Potential Energy [J/kg]"""
    longname="Surface Based Convective Available Potential Energy"; units="J/kg"
    clevs=[100,250,500,1000,1500,2000,2500,3000] #CAPE 
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname) 
    if(varname=='CAPE18000-0'):   grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="18000-0")[0].values 
    if(varname=='CAPE9000-0'):    grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="9000-0")[0].values 
    if(varname=='CAPE25500-0'):   grb=grbs.select(name="Convective available potential energy",typeOfLevel="pressureFromGroundLayer",level="25500-0")[0].values 

def plot_CAPE18000_0(var_n):
    """Convective Available Potential Energy 180mb - 0mb From Ground Layer [J/kg]"""
    longname="Convective Available Potential Energy 180mb - 0mb From Ground Layer"; units="J/kg"
    clevs=[100,250,500,1000,1500,2000,2500,3000] #CAPE 
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname) 

def plot_CAPE9000_0(var_n):
    """Convective Available Potential Energy 90mb - 0mb From Ground Layer [J/kg]"""
    longname="Convective Available Potential Energy 90mb - 0mb From Ground Layer"; units="J/kg"
    clevs=[100,250,500,1000,1500,2000,2500,3000] #CAPE 
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname) 

def plot_CAPE25500_0(var_n):
    """Convective Available Potential Energy 255mb - 0mb From Ground Layer [J/kg]"""
    longname="Convective Available Potential Energy 255mb - 0mb From Ground Layer"; units="J/kg"
    clevs=[100,250,500,1000,1500,2000,2500,3000] #CAPE 
    cm=plt.get_cmap(name='Spectral_r')
    return(var_n,clevs,cm,units,longname) 

def plot_REFD1km(var_n):
    """Stoelinga simulated base (1 km AGL) reflectivity"""
    longname="Stoelinga simulated base (1 km AGL) reflectivity"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 

def plot_REFD4km(var_n):
    """Stoelinga simulated base reflectivity"""
    longname="Stoelinga simulated base reflectivity"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 

def plot_REFDm10C(var_n):
    """Reflectivity at -10C level"""
    longname="Reflectivity at -10C level"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 

def plot_RETOP(var_n):
    """Echo top ( <= 18.5 dBz ) [m]"""
    longname="Echo top ( <= 18.5 dBz )"; units="kft"
    if(units=="kft"): var_n = var_n*304.8 #meters to kft 
    clevs=[0, 4, 7,10,13,16,19,22,25,28,31,34,37,40,43,46,49] # kft
    clist=[ 0,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 

def plot_SRHL0_3km_max(var_n):
    """0-3 km Storm Relative Helicity (max)"""
    longname="0-3 km Storm Relative Helicity (max)"; units="m/s**2"
    clevs=[0,25,50,75,100,150,200,250,300,400,500,600,700,800] # m/s**2 
    clist=[ 0, 4, 25,26,27, 23, 22, 21, 20, 18, 17, 15, 7 ]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 


def plot_SRHL0_1km_max(var_n):
    """0-1 km Storm Relative Helicity (max)"""
    longname="0-1 km Storm Relative Helicity (max)"; units="m/s**2"
    clevs=[0,25,50,75,100,150,200,250,300,400,500,600,700,800] # m/s**2 
    clist=[ 0, 4, 25,26,27, 23, 22, 21, 20, 18, 17, 15, 7 ]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 

def plot_MXUPHL2_5km(var_n):
    """2-5 km Updraft Helicity (max)"""
    longname="2-5 km Updraft Helicity (max)"; units="m/s**2"
    colors = ['lightgray','skyblue','mediumblue','green','orchid','firebrick','orangered','DarkViolet','black']
    clevs = [10.,       25.,      50.,         75.,    100.,  150.,       200.,      250.,        300.]
    cm = c.ListedColormap(colors)
    cm.set_over('white')
    #norm = c.BoundaryNorm(clevs, cm.N)
    return(var_n,clevs,cm,units,longname) 

def plot_MNUPHL2_5km(var_n):
    """2-5 km Updraft Helicity (min)"""
    longname="2-5 km Updraft Helicity (min)"; units="m/s**2"
    clevs=[0,25,50,75,100,150,200,250,300,400,500,600,700,800] # m/s**2 
    clist=[ 0, 4, 25,26,27, 23, 22, 21, 20, 18, 17, 15, 7 ]
    cm=gemplot(clist)
    return(var_n,clevs,cm,units,longname) 


def plot_MAXREFC(var_n):
    """Stoelinga simulated maximum (composite) reflectivity"""
    longname="Stoelinga simulated maximum (composite) reflectivity"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 

def plot_MAXREF_1km(var_n):
    """Stoelinga simulated maximum (composite) reflectivity"""
    longname="Stoelinga simulated maximum (composite) reflectivity"; units="dbz"
    clevs=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75] # dbz 
    cm=ncepy.mrms_radarmap()
    return(var_n,clevs,cm,units,longname) 



############### Dictionary for plot_function calls ###################################
def plot_Dictionary():
    #As fields are added to fv3 output just put those in the following dictionary
    #   according to the syntax used. Then all you have to do is create a function
    #   that defines the clevs, cm, and var_n if it requires unit conversion 
    #   (e.g., plot_PRATEsfc(var_n) )
    """The purpose of this dictionary is so that for each variable name (e.g., "ALBDOsfc") 
       the corresponding function is called (e.g., plot_ALBDOsfc(var_n)) to provide the 
       appropriate variable specific name, units, clevs, clist, and colormap for plotting.
    """
    dispatcher={  
# tracer
        'sphum':plot_sphum,
        'liq_wat':plot_liq_wat,
        'o3mr':plot_o3mr,
# core
        'u':plot_u,
        'v':plot_v,
        'W':plot_W,
        'T':plot_T,
        'DZ':plot_DZ,
        'delp':plot_delp,
        'phis':plot_phis,
# nggps 2d
        'ALBDOsfc':plot_ALBDOsfc,
        'CPRATsfc':plot_CPRATsfc,
        'PRATEsfc':plot_PRATEsfc,
        'DLWRFsfc':plot_DLWRFsfc,
        'ULWRFsfc':plot_ULWRFsfc,
        'DSWRFsfc':plot_DSWRFsfc,
        'USWRFsfc':plot_USWRFsfc,
        'DSWRFtoa':plot_DSWRFtoa,
        'USWRFtoa':plot_USWRFtoa,
        'ULWRFtoa':plot_ULWRFtoa,
        'GFLUXsfc':plot_GFLUXsfc,
        'HGTsfc':plot_HGTsfc,
        'HPBLsfc':plot_HPBLsfc,
        'ICECsfc':plot_ICECsfc,
        'SLMSKsfc':plot_SLMSKsfc,
        'LHTFLsfc':plot_LHTFLsfc,
        'SHTFLsfc':plot_SHTFLsfc,
        'PRESsfc':plot_PRESsfc,
        'PWATclm':plot_PWATclm,
        'SOILM':plot_SOILM,
        'SOILW1':plot_SOILW1,
        'SOILW2':plot_SOILW2,
        'SOILW3':plot_SOILW3,
        'SOILW4':plot_SOILW4,
        'SPFH2m':plot_SPFH2m,
        'SOILT1':plot_SOILT1,
        'SOILT2':plot_SOILT2,
        'SOILT3':plot_SOILT3,
        'SOILT4':plot_SOILT4,
        'TMP2m':plot_TMP2m,
        'TMPsfc':plot_TMPsfc,
        'UGWDsfc':plot_UGWDsfc,
        'VGWDsfc':plot_VGWDsfc,
        'UFLXsfc':plot_UFLXsfc,
        'VFLXsfc':plot_VFLXsfc,
        'UGRD10m':plot_UGRD10m,
        'VGRD10m':plot_VGRD10m,
        'WEASDsfc':plot_WEASDsfc,
        'SNODsfc':plot_SNODsfc,
        'ZORLsfc':plot_ZORLsfc,
        'VFRACsfc':plot_VFRACsfc,
        'F10Msfc':plot_F10Msfc,
        'VTYPEsfc':plot_VTYPEsfc,
        'STYPEsfc':plot_STYPEsfc,
        'TCDCclm':plot_TCDCclm,
        'TCDChcl':plot_TCDChcl,
        'TCDCmcl':plot_TCDCmcl,
        'TCDClcl':plot_TCDClcl,
        'REFC':plot_REFC,
        'APCP':plot_APCP,
        'APCP_C':plot_APCP_C,
        'MSL':plot_MSL,
        'CAPEsfc':plot_CAPEsfc,
        'REFD1km':plot_REFD1km,
        'REFD4km':plot_REFD4km,
        'REFDm10C':plot_REFDm10C,
        'RETOP':plot_RETOP,
        'MXUPHL2_5km_max':plot_MXUPHL2_5km,
        'MNUPHL2_5km_min':plot_MNUPHL2_5km,
        'MAXREFC_max':plot_MAXREFC,
        'MAXREF_1km_max':plot_MAXREF_1km,
               }
    return dispatcher  

def exp_Dictionary():
    experiment_dispatcher={  
        'rw_c008':'control',
        'rw_019':'w_only',
        'rw_021':'w_so_elev5',
        'rw_022':'w_so_elev10',
        'rw_023':'so_elev10',
        'obs':'observations',
                           }
    return experiment_dispatcher



if __name__ == '__main__':
    #pool=multiprocessing.Pool(len(varnames)) # one processor per variable
    #pool=multiprocessing.Pool(8) # 8 processors for all variables. Just a little slower.
    #pool.map(mkplot,varnames) 
    mkplot(varnames[0])
    toc=timer()
    time=toc-tic
    hrs=int(time/3600)
    mins=int(time%3600/60)
    secs=int(time%3600%60)
    print("Total elapsed time: "+str(toc-tic)+" seconds.")
    print("Total elapsed time: "+str(hrs).zfill(2)+":"+str(mins).zfill(2)+":"+str(secs).zfill(2))
