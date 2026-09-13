# ### 1. Importing necessary libraries
import requests
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import detrend
import os
import numpy as np
import math
import sys

basedir=os.getcwd()
#basedir='/home/anewman/GPS/Kivu'
#datadir='/home/anewman/GPS/KivuGNSS'
# processing and base info
procdir=os.path.join(basedir,'../KivuGNSS.processing')
bindir=os.path.join(procdir,'bin')
mapdir=os.path.join(procdir,'mapdata')
sys.path.append(os.path.join(procdir,'external_programs','euler_pole','euler_pole'))
from euler_pole import EulerPole as EP
from euler_pole import cart2sph
# output directories
datadir=os.path.join(basedir,'timeseries')
plotdir=os.path.join(basedir,'plots')
ratedir=os.path.join(basedir,'rates')

yearspernanosec=1/(86400e9*365.25) #convert numeric nanosecond into fractional year

# make necessary output directories if they don't exist
if not os.path.exists(datadir):
    os.mkdir(datadir)
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
if not os.path.exists(ratedir):
    os.mkdir(ratedir)

def wcorr(df):
    """
    df contains dN,E,U and wN,E,U for dataset to determine
    the weighed covarience matrix
    """
    Q=np.array(df[['dN','dE','dU']])
    W=np.array(df[['wN','wE','wU']])
    QW=Q*W
    C=QW.T.dot(QW)/W.T.dot(W)  # covariance
    CorrNE=C[0][1]/(C[0][0]*C[1][1])**0.5
    CorrNU=C[0][2]/(C[0][0]*C[2][2])**0.5
    CorrEU=C[1][2]/(C[1][1]*C[2][2])**0.5
    return CorrNE,CorrNU,CorrEU

def xyz2blh(x, y, z):
    """Convert XYZ coordinates to BLH,
    return tuple(latitude, longitude, height).
    
    from https://github.com/purpleskyfall/XYZ2BLH/blob/master/xyz2blh.py
    """
    A = 6378137.0
    B = 6356752.314245
    e = math.sqrt(1 - (B**2)/(A**2))
    # calculate longitude, in radians
    longitude = math.atan2(y, x)

    # calculate latitude, in radians
    xy_hypot = math.hypot(x, y)

    lat0 = 0
    latitude = math.atan(z / xy_hypot)

    while abs(latitude - lat0) > 1E-9:
        lat0 = latitude
        N = A / math.sqrt(1 - e**2 * math.sin(lat0)**2)
        latitude = math.atan((z + e**2 * N * math.sin(lat0)) / xy_hypot)

    # calculate height, in meters
    N = A / math.sqrt(1 - e**2 * math.sin(latitude)**2)
    if abs(latitude) < math.pi / 4:
        R, phi = math.hypot(xy_hypot, z), math.atan(z / xy_hypot)
        height = R * math.cos(phi) / math.cos(latitude) - N
    else:
        height = z / math.sin(latitude) - N * (1 - e**2)

    # convert angle unit to degrees
    longitude = math.degrees(longitude)
    latitude = math.degrees(latitude)

    return latitude, longitude, height

def getProcessedGNSS(Stats,analysisCenter='ovsi',refFrame='igs14',processing='Uncleanded',refCoord='from_analysis_center',datadir='timeseries'):
    '''
    Pull data from one of the automatic GNSS processing centers ('unr', or 'cwu') and save into relevant files.
    At present all flags above after analysisCenter only pertain to data pulled from 'cwu'.  All
    data pulled from 'unr' will assume defaults for their 'tenv' file format.
    
    Format 'tenv' is headerless.  Info on this is here: http://geodesy.unr.edu/gps_timeseries/README_tenv.txt

    Recently updated to use authentification tokens for UNAVCO: 9/27/2023
    Now using .pos files format version 1.1.1 (hopefully these remain more consistent.  Certainly much more data-rich)
    '''

    import requests
    from pathlib import Path
    from earthscope_sdk.auth.device_code_flow import DeviceCodeFlowSimple
    from earthscope_sdk.auth.auth_flow import NoTokensError
    import os
    # get, or update token if necessary
    #token_path='/Users/an77/tokens/'  # note that as of now, this will not work when you create your own file
    token_path= '/Users/lhylu/Library/Application Support/earthscope-cli/sso_tokens.json'  # note that as of now, this will not work when you create your own file
    print(token_path)

    device_flow = DeviceCodeFlowSimple(Path(token_path))
    try:
        # get access token from local path
        device_flow.get_access_token_refresh_if_necessary()
    except NoTokensError:
        # if no token was found locally, do the device code flow
        device_flow.do_flow()
    token = device_flow.access_token
    #file open(token_path,'r')
    #token=file.read()

    nStats=len(Stats)
    count=0
    # # if not os.path.exists(datadir):
    #     os.mkdir(datadir)
    
    for stat in Stats:
        count += 1
        print('Requesting data for: '+stat+'  ['+str(count)+'/'+str(nStats)+']') # usable output

        if analysisCenter == 'cwu':
            #URL='https://web-services.unavco.org/gps/data/position/' + stat +    '/v3?analysisCenter=' + analysisCenter +    '&referenceFrame=' + refFrame +    '&report=long&dataPostProcessing=' + processing +    '&refCoordOptio='+ refCoord
            #URL='https://data.unavco.org/archive/gnss/products/position/' + stat + '/' + stat + '.' + analysisCenter +  '.' + refFrame + '.csv'
            URL='https://data.unavco.org/archive/gnss/products/position/' + stat + '/' + stat + '.' + analysisCenter +  '.' + refFrame + '.pos'
            filename=os.path.join(datadir,stat + '.' + analysisCenter +  '.' + refFrame + '.pos')
        #     req = requests.get(URL, headers={"authorization": f"Bearer {token}"})
        elif analysisCenter == 'unr':
            URL='http://geodesy.unr.edu/gps_timeseries/tenv/IGS14/' + stat +    '.tenv'
            filename=os.path.join(datadir,stat+'.tenv')
            req = requests.get(URL) # get data
        elif analysisCenter == 'ovsi':
            #datadir = os.path.join(basedir, 'CR_GNSS')  # Update to use the timeseries folder
            #filename = os.path.join(datadir, stat + '.ovs.final_igb14.pos')
            src_path = os.path.join(basedir, 'CR_GNSS', stat + '.ovs.final_igb14.pos')
            dst_path = os.path.join('timeseries', stat + '.ovs.final_igb14.pos')
            shutil.copy2(src_path, dst_path)
            print(f"Copied {src_path} → {dst_path}")
        
        else:
            print("ERROR:   analysic Center '"+ analysisCenter+"' not recognized.  Exiting.")
            exit(-1)

        # if req.status_code == requests.codes.ok:
        #     # save the file
        #     # print("  File URL: "+ URL)
        #     url_content = req.content
        #     file = open(filename, 'wb')
        #     file.write(url_content)
        #     file.close()
        # else:
        #     #problem occured
        #     print(f"failure: {req.status_code}, {req.reason}")


def GNSSTimeSeries2Pandas(stat,analysisCenter='cwu',datadir='timeseries',refFrame='igs14'):
    '''
    Convert data pulled from an automatic GNSS processing center ('unr', or 'cwu') into a standard pandas format.
    
    The 'cwu' csv file has poor headers for Pandas, and the 'tenv' format is headerless.  
    Info on the 'tenv' format is here: http://geodesy.unr.edu/gps_timeseries/README_tenv.txt
    '''
    if not os.path.exists(datadir):
        print("Error:  Cannot find '"+datadir+"'.  Exiting.")
        exit(-1)
    if analysisCenter == 'cwu':
        filename=os.path.join(datadir,stat+'.'+analysisCenter+'.'+refFrame+'.pos')
        # find length of header
        with open(filename,'r') as file:
            for num, line in enumerate(file,1):
                if "End Field Description" in line:
                    hlength=num
                    break
        df = pd.read_csv(filename,delim_whitespace=True,header=hlength)
        df = df.rename(columns = {str('*YYYYMMDD'):str('YYYYMMDD')})
        df['Date']=pd.to_datetime(df['YYYYMMDD'], format='%Y%m%d')
        
        Comps=('N','E', 'U')
        for comp in Comps:
            df = df.rename(columns = {str('d'+comp):str(comp+'pos')}) 
            df = df.rename(columns = {str('S'+comp.lower()):str(comp+'err')})
        Cors=('NE','EU', 'NU')
        for cor in Cors:
            df = df.rename(columns = {str('R'+cor.lower()):str(cor+'cor')})   # lowercase in input
    elif analysisCenter == 'unr':
        filename=os.path.join(datadir,stat+'.tenv')
        df = pd.read_csv(filename, header=None, delim_whitespace=True,
           names = ['YYYYMMDD','HHMMSS','JJJJJ.JJJJ','X','Y','Z','Sx','Sy','Sz',
                   'Rxy','Rxz','Ryz','NLat','Elong','Height','dN', 'dE', 'dU', 'Sn', 'Se', 'Su', 'Rne', 'Rnu', 'Reu', 'Soln' ]
           )
        df['Date'] = pd.to_datetime(df['YYYYMMDD'], format='%Y%m%d')
        print('Date')


    elif analysisCenter == 'ovsi':
        print("Testing ovsi pandas")
        datadir = '/Users/lhylu/KivuGNSS.processing/CR_GNSS'
        filename=os.path.join(datadir,stat+'.ovs.final_igb14.pos')
        df = pd.read_csv(filename, header=None, delim_whitespace=True,
           names = ['STAT','Date','DecimalYear','MJD','GPSWeek','GPSWeekDay','Epos','Npos','Upos',
                   'AntHeight','Eerr','Nerr','Uerr','NEcor','EUcor','NUcor']
           )
        df['Date']=pd.to_datetime(df.Date, yearfirst=True)
    else:
        print("ERROR:   analysic Center '"+ analysisCenter+"' not recognized.  Exiting.")
        exit(-1)
    df['NumDate']=pd.to_numeric(df.Date)*yearspernanosec # converts to date in NANOseconds
    return df 

def getUNRlocs(stat): 
    '''
    directly pulls info from their station information page
    '''
    URL='http://geodesy.unr.edu/NGLStationPages/stations/'+stat+'.sta'
    req = requests.get(URL)
    #statinfo = req.content
    # read all lines using readline()
    lat=lon=height=-9999.9  # defaults
    for row in req.iter_lines():
       # check if string present on a current line 
                      
       row=str(row)
       if row.find('Latitude') != -1:
          lat=float(row.split()[2])
       if row.find('Longitude') != -1:
          lon=float(row.split()[2])
       if row.find('Height') != -1:
          height=float(row.split()[2])
    return lat, lon, height

def getCWUlocs(stat,analysisCenter='cwu',refFrame='igs14'): 
    '''
    directly pulls info from their station information page
    '''
    if analysisCenter == 'unr':
        #creating pandas readable csv from the url
        file=os.path.join(datadir,stat+'.'+analysisCenter+'.'+refFrame+'.pos')
        lat=lon=height=-9999.9  # defaults
        # get xyz coords from header (reread file)
        #  worked before version 1.2 of file structure
        # csv_file=open(file,'r')
        # content = csv_file.readlines()
        # csv_file.close()
        # xyzline=content[6]
        # words=xyzline.split()
        # xyz=[float(words[5]),float(words[7]),float(words[9])] # xyz coordinates
        # lat,lon,height=xyz2blh(xyz[0],xyz[1],xyz[2])
    
        with open(file,'r') as csv_file:
            for line in csv_file:
                if line.startswith("NEU Reference position"):
                    words=line.split()
                    lat,lon,height=float(words[4]),float(words[5]),float(words[6])
                    break
    elif analysisCenter == 'ovsi':
        datadir = '/Users/lhylu/KivuGNSS.processing/CR_GNSS'
        file = os.path.join(datadir, stat + '.ovs.final_igb14.pos')
        df = pd.read_csv(file, header=None, sep='\s+', comment='*',
           names = ['YYYYMMDD','HHMMSS','JJJJJ.JJJJ','X','Y','Z','Sx','Sy','Sz',
                   'Rxy','Rxz','Ryz','NLat','Elong','Height','dN', 'dE', 'dU', 'Sn', 'Se', 'Su', 'Rne', 'Rnu', 'Reu', 'Soln' ]
           )
        df['Date'] = pd.to_datetime(df['YYYYMMDD'], format='%Y%m%d')
        df['lat'] = df['X']
        df['lon'] = df['Y']
        df['height'] = df['Z']
        print(lat, lon)


    return lat,lon,height

def ITRF14EulerPole(Plate):
    '''
    Will output euler pole as defined by ITRF-14.  Will report error if it doesn't reognize the plate.
    Described plates are:
          ANTA, ARAB, AUST, EURA, INDI, NAZC, NOAM, NUBI, PCFC, SOAM, SOMA
    ''' 
    ITRFplatesfile = os.path.join(mapdir,'ITRF14','Table-S1.txt')
    ITRFplates = pd.read_csv(ITRFplatesfile, delim_whitespace=True, comment='#')
    try:
        ITRFplate =ITRFplates[ITRFplates["plate"]==Plate].values
        Plat,Plon,Prot = cart2sph(ITRFplate[0,1],
                                 ITRFplate[0,2],
                                 ITRFplate[0,3]) # cart to geog
    except:
        print("local Euler Pole Rates not described") 
    Pole=EP(Plat,Plon,Prot)
    return Pole

def SariaEulerPole(lPlate):
    '''
    Gets local Euler Pole relative to the Nubian Plate as described within 
    Saria et al., JGR, 2014.
    '''
    lplatesfile = os.path.join(mapdir,'Saria_etal_JGR_2014_table.txt')
    try: 
        lplates = pd.read_csv(lplatesfile, delim_whitespace=True, comment='#')
        lplate =lplates[lplates["Plate"]==lPlate].values
        lPlat,lPlon,lProt =(float(lplate[0,1]), float(lplate[0,2]), float(lplate[0,3]))
        lPole=EP(lPlat,lPlon,lProt)
    except:
        print("local Euler Pole Rates not described") 
    return lPole

def parseTimeSeries(df,analysisCenter='cwu'):
    Comps=('N','E', 'U')
    for comp in Comps:
        ycolComp=str(comp+'pos')
        df[ycolComp]=df[ycolComp]-df[ycolComp].mean()
    # establish rapid and final datasets 
    if analysisCenter == 'unr':
        dfFinal=df  # only final are used here
        dfRapid=pd.DataFrame()
    elif analysisCenter == 'cwu': 
        dfF1=df.loc[df['Soln']=='final']          #final solutions
        dfF2=df.loc[df['Soln']=='suppl']          #final solutions
        dfF3=df.loc[df['Soln']=='suppf']          #final solutions
        dfFinal=pd.concat([dfF1,dfF2,dfF3]).sort_index()
        dfRapid=df.loc[df['Soln']=='rapid'].reset_index()   #rapid solutions
    return dfFinal, dfRapid

def linearRates(Stats,ratesfile,platename,Euler,analysisCenter,yscale=1000):
    """
    Create velocity determinations from the data that is pulled from either the 'cwu' or 'unr' analysis centers.
    Requires a Stats list of stations, the name of the ratesfile, the Euler pole to remove from the data, 
    as well as a scale for plotting [default of yscale=1000 converts m to mm].
    In addition to writing the ratesfile, process will return a dictionary of fit values useful for futher plotting, etc.
    """
    ratesFile=open(os.path.join(ratedir,ratesfile), 'w')  # new velocity file
    ratesFile.write('# local Plate motion is defined as %s Plate \n' % (platename))
    ratesFile.write('STAT        Lat      Long   Height      Nvel     Evel     Nvel-loc Evel-loc Uvel      Nerr     Eerr     Uerr    NEcor   NUcor   EUcor         Sdate       Edate\n') 
    ratesFile.write('#           °        °      m          mm/yr    mm/yr    mm/yr    mm/yr    mm/yr     mm/yr    mm/yr    mm/yr                             YEAR-MO-DY  YEAR-MO-DY\n') 
    ratesFile.write('#--------------------------------------------------------------------------------------------------------------------------------------------------------------\n') 
    fitDict = {}  # dictionary of fit values for individual stations and components.  this info is returned for plotting later
    for stat in Stats:
        #creating pandas readable csv from the url
        if analysisCenter == 'unr':
            lat, lon, height = getUNRlocs(stat)  
        elif analysisCenter == 'cwu':
            lat, lon, height = getCWUlocs(stat)  
        elif analysisCenter == 'ovsi':
            lat, lon, height = getCWUlocs(stat)  
        else:
            print("Could not get coordinates from "+analysisCenter+" for "+stat+", setting to defaults.")
        df=GNSSTimeSeries2Pandas(stat,analysisCenter=analysisCenter) 
        dfFinal,dfRapid=parseTimeSeries(df,analysisCenter=analysisCenter)
    
        xFND=dfFinal['NumDate']
        sp=0
        slopes=np.zeros(3) # store slopes and errors
        errs=np.zeros(3)
        dfFNEU=pd.DataFrame()
        Comps=('N','E', 'U')
        for comp in Comps:
            ycolComp=str(comp+'pos')
            wcolComp=str(comp+'err')
            sdate=dfFinal.Date[0].strftime("%Y-%m-%d")
            edate=dfFinal.Date[len(dfFinal)-1].strftime("%Y-%m-%d")
            yF=dfFinal[ycolComp]*yscale
            wghtF=1/dfFinal[wcolComp]/yscale
            
            # performing linear fits to data
            fit,var=np.polyfit(xFND,yF,1,w=wghtF, full=False, cov=True)
            fitDict[stat+'-'+comp+'-fit'] = fit
            fitDict[stat+'-'+comp+'-var'] = var
            fity=np.polyval(fit,xFND)
            errs[sp]=np.sqrt(var[0][0])
            slopes[sp]=fit[0]
            dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
            dfFNEU['w'+comp]=(wghtF)   #  weights 
            sp+=1
        CorrNE,CorrNU,CorrEU=wcorr(dfFNEU)
        pAZ,pRate = Euler.velocity(lat,lon) 
        pNvel = pRate * np.cos(np.radians(pAZ)) 
        pEvel = pRate * np.sin(np.radians(pAZ))
        #print("Euler Pole: %s (%5.1f,%6.1f): N= %5.1f, E= %5.1f (Mag = %5.1f, %6.1f°)"
        #     %(stat,lat,lon,pNvel,pEvel,pRate,pAZ))
        ratesFile.write('%s   %8.4f %9.4f %8.1f  %8.1f %8.1f %8.1f %8.1f %8.1f  %8.1f %8.1f %8.1f  %7.4f %7.4f %7.4f    %s  %s\n' % 
                        (stat, lat, lon, height, 
                        slopes[0], slopes[1],  # in NNR
                        slopes[0]-pNvel, slopes[1]-pEvel,    # in SO plate ref.
                        slopes[2],    # vert
                        errs[0],  errs[1], errs[2],   # errors
                        CorrNE, CorrNU, CorrEU, 
                        sdate, edate)) 
    ratesFile.close()
    return fitDict
       
def timeSeriesPlots(Stats,fitDict,analysisCenter,yscale=1000):
    """
    Create displacement timeseries from the data that is pulled from either the 'cwu' or 'unr' analysis centers.
    Requirs a Stats list of stations, the fit dictionary for those stations (created in linearRates() function), 
    as well as a scale for plotting [default of yscale=1000 converts m to mm].
    """
    # Plotting design options
    # errors
    csize=2; elw=0.8; ecol='k' # head size, width, color
    # markers 
    mec='k'; mew=0.8; mfmt='o'; msz='4'  #edge color, width, shape, size
    # final results
    falpha=1; fcol='blue'  # opacity, color
    # rapid results
    ralpha=0.5; rcol='red'  # opacity, color
    # grid
    galpha=0.5; gvwidth=0.5; ghwidth=1; gcol='gray'
    for stat in Stats:
        df=GNSSTimeSeries2Pandas(stat,analysisCenter=analysisCenter) 
        # split data into final and rapid, removing mean along the way
        dfFinal,dfRapid=parseTimeSeries(df,analysisCenter=analysisCenter)
        
        #creating Station (3 component) plots 
        f, ax = plt.subplots(3, 1, figsize=(14,8),sharey=False, sharex=True) 
        f.tight_layout(h_pad=0)
        
        xF=dfFinal['Date']
        xFND=dfFinal['NumDate']
        if len(dfRapid) > 0 : 
            xR=dfRapid['Date']
        sp=0
        slopes=np.zeros(3) # store slopes and errors
        errs=np.zeros(3)
        dfFNEU=pd.DataFrame()
    
        Comps=('N','E', 'U')
        for comp in Comps:
            ycolComp=str(comp+'pos')
            wcolComp=str(comp+'err')
            yF=dfFinal[ycolComp]*yscale
            yeF=dfFinal[wcolComp]*yscale
            wghtF=1/dfFinal[wcolComp]/yscale
            if len(dfRapid) > 0 : 
                yR=dfRapid[ycolComp]*yscale
                yeR=dfRapid[wcolComp]*yscale
            
            # performing linear fits to data
            fit= fitDict[stat+'-'+comp+'-fit']
            var= fitDict[stat+'-'+comp+'-var']
            errs[sp]=np.sqrt(var[0][0])
            fity=np.polyval(fit,xFND)
            slopes[sp]=fit[0]
            dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
            dfFNEU['w'+comp]=(wghtF)   #  weights 

            #plot
            ax[sp].plot(xF, fity)
            ax[sp].errorbar(x=xF, y=yF, yerr=yeF, 
                fmt=mfmt, ms=msz, capsize=csize, label='Final', mfc=fcol, mec=mec, mew=mew, 
                ecolor=ecol, elinewidth=elw, alpha=falpha)
            if len(dfRapid) > 0 : 
                ax[sp].errorbar(x=xR, y=yR, yerr=yeR, 
                    fmt=mfmt, ms=msz, capsize=csize, label='Rapid', mfc=rcol, mec=mec, mew=mew, 
                    ecolor=ecol, elinewidth=elw, alpha=ralpha)
            #plot labels and legend
            ax[sp].grid(axis='x', linestyle='-', color=gcol, linewidth=gvwidth, alpha=galpha)
            ax[sp].axhline(0, linestyle='-', color=gcol,linewidth=ghwidth, alpha=galpha) 
            if sp == 0 :
                ax[sp].legend(loc='upper left',fancybox=True, shadow=True)
                ax[sp].set_ylabel('North [mm]')
            if sp == 1 :
                ax[sp].set_ylabel('East [mm]')
            elif sp == 2 :
                ax[sp].set_ylabel('Vertical [mm]')
                ax[sp].set_xlabel('Date')
                [xmin,xmax,ymin,ymax]=plt.axis()
                ax[sp].text(xmin+(xmax-xmin)*.005, ymin+(ymax-ymin)*.01, 'Daily positions processed @ '+analysisCenter.upper(), 
                    horizontalalignment='left', 
                    verticalalignment='bottom')
                # put title atop as a last thing (includes values)
                ax[0].set_title('Kivu Rift Project: '+stat+'  Rates: N=%.1f±%.1f, E=%.1f±%.1f, V=%.1f±%.1f [mm/yr]' % (slopes[0], errs[0],slopes[1], errs[1],slopes[2], errs[2]))
            sp+=1
        f.savefig(os.path.join(plotdir,stat+'_TS.png'), dpi=150, facecolor='white', bbox_inches='tight', pad_inches=0.5)#!/usr/bin/env python
        plt.close('all')


def UKF(Stats,fitDict,analysisCenter,yscale=1000,UKFerr=[3,0.01,5]):
    """
    Similar to timeSeriesPlots, will create a plot of the time series data, but will also show the continuously
    updated Unscented Kalman Filter using the prior results.
    Requires a Stats list of stations, the fit dictionary for those stations (created in linearRates() function), 
    the analysiscenter (cwu, or unr),
    as well as a scale for plotting [default of yscale=1000 converts m to mm].

    This version, as implemented, will predict the 3D position (N,E,U), given the invididual data, their standard
    error, and covarience. 

    UKFerr are  [std measurement error, std process error, initial uncertainty in mm] 


    adaptation for GPS data using instructions from  below
    https://filterpy.readthedocs.io/en/latest/kalman/UnscentedKalmanFilter.html
    """
    from filterpy.kalman import MerweScaledSigmaPoints, UnscentedKalmanFilter 
    from filterpy.common import Q_discrete_white_noise
    
    tscale=365.25 # converting results from per day to per/year for velocities
    
    basedir=os.getcwd()
    # output directories
    plotdir=os.path.join(basedir,'plots')
    
    dt = 1 # 1 day 
    dim_z = 3 # x,y,z
    dim_x = 2 * dim_z   # state variables, pos + vel for dim_z components
    
    z_std = UKFerr[0]  # std of measurement error
    q_std = UKFerr[1]  # std of process error
    startP = UKFerr[2] # initial uncertainty  (in mm)
    
    # plotting params
    # errors
    csize=2; elw=0.8; ecol='k' # head size, width, color
    # markers 
    mec='k'; mew=0.8; mfmt='o'; msz='4'  #edge color, width, shape, size
    # final results
    falpha=1; fcol='blue'  # opacity, color
    # rapid results
    ralpha=0.5; rcol='red'  # opacity, color
    # grid
    galpha=0.5; gvwidth=0.5; ghwidth=1; gcol='gray' 
     
    
    def fx(x,dt=dt):
        # build state transition function for constant velocity (change later)
        F = np.array([
            [1,dt,0,0,0,0],  # velocity x
            [0, 1,0,0,0,0],  # position x
            [0,0,1,dt,0,0],  # velocity y
            [0,0,0, 1,0,0],  # position y
            [0,0,0,0,1,dt],  # velocity z
            [0,0,0,0,0, 1]   # position z
            ], dtype=float)
        return np.dot(F,x)
    def hx(x):
        return np.array([x[0],x[2],x[4]])  # 3 component
    # create sigma points to use in the filter. This is standard for Gaussian processes
    points = MerweScaledSigmaPoints(dim_x, alpha=.1, beta=2., kappa=0)
    
    for stat in Stats:
        df=GNSSTimeSeries2Pandas(stat,analysisCenter=analysisCenter) 
        # split data into final and rapid, removing mean along the way
        # parse and remove mean
        Comps=('N','E', 'U')
        for comp in Comps:
            df[str(comp+"pos")] *= yscale  # convert components to mm
            df[str(comp+"err")] *= yscale

        dfFinal,dfRapid=parseTimeSeries(df,analysisCenter=analysisCenter)
        
        kf = UnscentedKalmanFilter(dim_x=dim_x, dim_z=dim_z, dt=dt, fx=fx, hx=hx, points=points)
        kf.P *= startP
        kf.x = np.array([df["Npos"][0],0,df["Epos"][1],0,df[ "Upos"][1],0], dtype=float)
        kf.R = np.diag([z_std**2, z_std**2, z_std**2]) # measurement noise matrix
        kf.Q = Q_discrete_white_noise(dim=dim_z, dt=dt, var=q_std**2, block_size=2) # process noise matrix
        #kf.inv = np.linalg.pinv  # force a specific inversion tool
        Estdf=pd.DataFrame()
        for index, row in df.iterrows():
            z = np.array([row["Npos"],row[ "Epos"],row[ "Upos"]])
            varN=row["Nerr"]**2
            varE=row["Eerr"]**2
            varU=row["Uerr"]**2
            #R = np.diag([varN,varE,varU]) # ignoring data covariance
            covNE=row['NEcor']*row["Nerr"]*row["Eerr"]
            covNU=row['NUcor']*row["Nerr"]*row["Uerr"]
            covEU=row['EUcor']*row["Eerr"]*row["Uerr"]
            # include full 3D covariance matrix
            R = np.array([
                [varN,covNE,covNU],  
                [covNE,varE,covEU],  
                [covNU,covEU,varU]  
                 ], dtype=float)
            
            kf.predict()  # predicts next data point
            kf.update(z,R=R)  # updates parameters with true values at time predicted
            pdrow = pd.DataFrame({
                'Nest':[kf.x[0]], 'Eest':[kf.x[2]], 'Uest':[kf.x[4]],    # estimated positions
                'Nvest':[kf.x[1]*tscale],'Evest':[kf.x[3]*tscale],'Uvest':[kf.x[5]*tscale],  # estimated velocities
                'LogLikelihood':[kf.log_likelihood],
                'Mahalanobis':[kf.mahalanobis],
                })
            Estdf=pd.concat([Estdf,pdrow])
            #print(kf.x,kf.log_likelihood)
        #print(kf.log_likelihood,np.array([kf.x[1],kf.x[3],kf.x[5]])*tscale)

        #creating Station (3 component) plots 
        f, ax = plt.subplots(3, 1, figsize=(14,8),sharey=False, sharex=True) 
        f.tight_layout(h_pad=0)
        
        x=df['Date']
        xF=dfFinal['Date']
        xFND=dfFinal['NumDate']
        if len(dfRapid) > 0 : 
          xR=dfRapid['Date']
        sp=0
        slopes=np.zeros(3) # store slopes and errors
        errs=np.zeros(3)
        dfFNEU=pd.DataFrame()
        
        for comp in Comps:
            ycolComp=str(comp+'pos')
            wcolComp=str(comp+'err')
            ecolComp=str(comp+'est')
            yF=dfFinal[ycolComp]
            if len(dfRapid) > 0 : 
                yR=dfRapid[ycolComp]
                yeR=dfRapid[wcolComp]
            yeF=dfFinal[wcolComp]
            wghtF=1/dfFinal[wcolComp]
            
            # performing linear fits to data
            fit= fitDict[stat+'-'+comp+'-fit']
            var= fitDict[stat+'-'+comp+'-var']
            errs[sp]=np.sqrt(var[0][0])
            fity=np.polyval(fit,xFND)
            slopes[sp]=fit[0]
            dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
            dfFNEU['w'+comp]=(wghtF)   #  weights 
            
            #plot
            ax[sp].plot(xF, fity, label='linear fit')
            ax[sp].plot(x, Estdf[ecolComp],'r-', label='constant-Vel UKF', alpha=0.5)
            ax[sp].errorbar(x=xF, y=yF, yerr=yeF, 
                fmt=mfmt, ms=msz, capsize=csize, label='Final', mfc=fcol, mec=mec, mew=mew, 
                ecolor=ecol, elinewidth=elw, alpha=falpha)
            if len(dfRapid) > 0 : 
                ax[sp].errorbar(x=xR, y=yR, yerr=yeR, 
                    fmt=mfmt, ms=msz, capsize=csize, label='Rapid', mfc=rcol, mec=mec, mew=mew, 
                    ecolor=ecol, elinewidth=elw, alpha=ralpha)
            #ax[sp].text(right,bottom,, ha='bottom',va='right')
            
            #plot labels and legend
            ax[sp].grid(axis='x', linestyle='-', color=gcol, linewidth=gvwidth, alpha=galpha)
            ax[sp].axhline(0, linestyle='-', color=gcol,linewidth=ghwidth, alpha=galpha) 
            if sp == 0 :
                ax[sp].legend(loc='upper left',fancybox=True, shadow=True)
                ax[sp].set_ylabel('North [mm]')
            if sp == 1 :
                ax[sp].set_ylabel('East [mm]')
            elif sp == 2 :
                ax[sp].set_ylabel('Vertical [mm]')
                ax[sp].set_xlabel('Date')
                [xmin,xmax,ymin,ymax]=plt.axis()
                ax[sp].text(xmin+(xmax-xmin)*.005, ymin+(ymax-ymin)*.01, 'Daily positions processed @ '+analysisCenter.upper(), 
                    horizontalalignment='left', 
                    verticalalignment='bottom')
                ax[0].set_title('Kivu Rift Project: '+stat+'  Rates: N=%.1f±%.1f, E=%.1f±%.1f, V=%.1f±%.1f [mm/yr]' % (slopes[0], errs[0],slopes[1], errs[1],slopes[2], errs[2]))
            sp+=1
        f.savefig(os.path.join(plotdir,stat+'_UKF.png'), dpi=150, facecolor='white', bbox_inches='tight', pad_inches=0.5)
        plt.close('all')


def commonMode(Stats,commonratesfile,analysisCenter,meanList=None, yscale=1000,**kwargs):
    """
    Create a local velocity determination removing the common-mode signal from each daily solution.
    Data is pulled from either the 'cwu' or 'unr' analysis centers.
    Requires a Stats list of stations, the name of the commonratesfile, 
    as well as a scale for plotting [default of yscale=1000 converts m to mm].
    In addition to writing the ratesfile, process will return a dictionary of fit values useful for futher plotting, etc.
    """
    ratesFile=open(os.path.join(ratedir,commonratesfile), 'w')  # new velocity file
    ratesFile.write('# local Plate motion is defined by data meanList: \n')
    stringMeanList=" ".join(str(x) for x in meanList)
    ratesFile.write('# '+stringMeanList+' \n')
    ratesFile.write('STAT        Lat      Long   Height      Nvel     Evel    Uvel      Nerr     Eerr     Uerr    NEcor   NUcor   EUcor         Sdate       Edate\n') 
    ratesFile.write('#           °        °      m          mm/yr    mm/yr   mm/yr     mm/yr    mm/yr    mm/yr                             YEAR-MO-DY  YEAR-MO-DY\n') 
    ratesFile.write('#-------------------------------------------------------------------------------------------------------------------------------------------------\n') 
    fitDict = {}  # dictionary of fit values for individual stations and components.  this info is returned for plotting later

    # build composite catalog ordered by date
    dfPrimary=pd.DataFrame()
    for stat in Stats:
        df=GNSSTimeSeries2Pandas(stat,analysisCenter=analysisCenter) 
        df['stat']=stat
        df.set_index('Date', drop=False, inplace=True)
        #df.drop(['Datetime', ' X',' Y',' Z',' Std Dev X',' Std Dev Y',' Std Dev Z',' Corr XY',' Corr XZ',' Corr YZ', 
        #' N latitude', ' E longitude', ' Height'], axis=1, inplace=True)
        df=df[['stat','Date','NumDate','Npos','Epos','Upos','Nerr','Eerr','Uerr']] # unneccessary, but will save memory?
        dfPrimary=pd.concat([dfPrimary,df])
    idxold=None
    dfPrimary.sort_index(inplace=True)
    for idx,row in dfPrimary.sort_index().iterrows():
        if idx != idxold:  # only unique dates
            if meanList:
                Nmean=dfPrimary.loc[dfPrimary['stat'].isin(meanList)].loc[idx,'Npos'].mean()
                Emean=dfPrimary.loc[dfPrimary['stat'].isin(meanList)].loc[idx,'Epos'].mean()
                Umean=dfPrimary.loc[dfPrimary['stat'].isin(meanList)].loc[idx,'Upos'].mean()
            else:
                Nmean=dfPrimary['Npos'][idx].mean()
                Emean=dfPrimary['Epos'][idx].mean()
                Umean=dfPrimary['Upos'][idx].mean()
            # Add means to the primary database        
            dfPrimary.loc[idx,'Nmean']=Nmean
            dfPrimary.loc[idx,'Emean']=Emean
            dfPrimary.loc[idx,'Umean']=Umean
        idxold=idx

    dfPrimary['Nposl']=dfPrimary['Npos']-dfPrimary['Nmean']   
    dfPrimary['Eposl']=dfPrimary['Epos']-dfPrimary['Emean']   
    dfPrimary['Uposl']=dfPrimary['Upos']-dfPrimary['Umean']   
    for stat in Stats:
        #creating pandas readable csv from the url
        if analysisCenter == 'unr':
            lat, lon, height = getUNRlocs(stat)  
        elif analysisCenter == 'cwu':
            lat, lon, height = getCWUlocs(stat)  
        else:
            print("Could not get coordinates from "+analysisCenter+" ofor "+stat+", setting to defaults.")
        df=dfPrimary.loc[dfPrimary['stat'] == stat]
        sdate=df.Date[0].strftime("%Y-%m-%d")
        edate=df.Date[len(df)-1].strftime("%Y-%m-%d")
        xFND=df['NumDate']
        sp=0
        slopes=np.zeros(3) # store slopes and errors
        errs=np.zeros(3)
        dfFNEU=pd.DataFrame()
        Comps=('N','E', 'U')
        for comp in Comps:
            ycolComp=str(comp+'posl')
            ycolComp0=str(comp+'pos')
            wcolComp=str(comp+'err')
            mcolComp=str(comp+'mean')
            EcolComp=str(comp+'errl')
            yF=df[ycolComp]*yscale
            #wghtF=(df[wcolComp]**-2)/yscale
            # below is an error estimate using the highly correlated difference between the daily positions and the mean
            meanErr=detrend(df[mcolComp]).std()  # a singular value
            compErr=df[wcolComp]  # a series 
            crosscorr = np.corrcoef(df[ycolComp0],df[mcolComp])[0][1]
            #  sigma_x-y = sqrt(sigma_x^2 + sigma_y^2 - 2 corr *  sigma_x * sigma_y)
            Err=(compErr**2 + meanErr**2 -2* crosscorr*compErr*meanErr)**0.5
            #print("meanErr=",meanErr)
            #print("compErr=",compErr)
            #print("crosscorr",crosscorr)
            wghtF=(Err**-2)/yscale # (1/sigma^2)

            # performing linear fits to data
            fit,var=np.polyfit(xFND,yF,1,w=wghtF, full=False, cov=True)
            fitDict[stat+'-'+comp+'-fit'] = fit
            fitDict[stat+'-'+comp+'-var'] = var
            fity=np.polyval(fit,xFND)
            errs[sp]=np.sqrt(var[0][0])
            slopes[sp]=fit[0]
            dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
            dfFNEU['w'+comp]=(wghtF)   #  weights 
            #df= pd.concat([df,pd.Series(Err).rename(EcolComp)], axis=1)
            #print(df[EcolComp].head())
            dfPrimary.loc[dfPrimary['stat']==stat,EcolComp]=pd.Series(Err)
            sp+=1
        CorrNE,CorrNU,CorrEU=wcorr(dfFNEU)
        ratesFile.write('%s   %8.4f %9.4f %8.1f  %8.2f %8.2f %8.2f  %8.2f %8.2f %8.2f  %7.4f %7.4f %7.4f    %s  %s\n' % 
                        (stat, lat, lon, height, 
                        slopes[0], slopes[1],  slopes[2],    #  N, E, U
                        errs[0],  errs[1], errs[2],   # errors
                        CorrNE, CorrNU, CorrEU, 
                        sdate, edate)) 
    ratesFile.close()
    TSfile = os.path.join(datadir,'Primary_ts.csv')
    dfPrimary.to_csv(TSfile)
    return dfPrimary, fitDict

def timeSeriesPlotsCommon(dfPrimary,fitDict,analysisCenter,meanList,yscale=1000):
    """
    Create displacement timeseries from the data that is pulled from either the 'cwu' or 'unr' analysis centers.
    Requirs a Stats list of stations, the fit dictionary for those stations (created in linearRates() function), 
    as well as a scale for plotting [default of yscale=1000 converts m to mm].
    """
    # Plotting design options
    # errors
    csize=2; elw=0.8; ecol='k' # head size, width, color
    # markers 
    mec='k'; mew=0.8; mfmt='o'; msz='4'  #edge color, width, shape, size
    # final results
    falpha=1; fcol='blue'  # opacity, color
    # rapid results
    ralpha=0.5; rcol='red'  # opacity, color
    # grid
    galpha=0.5; gvwidth=0.5; ghwidth=1; gcol='gray'
    for stat in dfPrimary['stat'].unique(): # Stats:
        df=dfPrimary.loc[dfPrimary['stat'] == stat]
        
        #creating Station (3 component) plots 
        f, ax = plt.subplots(3, 1, figsize=(14,8),sharey=False, sharex=True) 
        f.tight_layout(h_pad=0)
        
        xF=df['Date']
        xFND=df['NumDate']
        sp=0
        slopes=np.zeros(3) # store slopes and errors
        errs=np.zeros(3)
        dfFNEU=pd.DataFrame()
    
        Comps=('N','E', 'U')
        for comp in Comps:
            ycolComp=str(comp+'posl')
            wcolComp=str(comp+'errl')
            yF=df[ycolComp]*yscale
            yeF=df[wcolComp]*yscale
            wghtF=df[wcolComp]**-2/yscale
            
            # performing linear fits to data
            fit= fitDict[stat+'-'+comp+'-fit']
            var= fitDict[stat+'-'+comp+'-var']
            errs[sp]=np.sqrt(var[0][0])
            fity=np.polyval(fit,xFND)
            slopes[sp]=fit[0]
            dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
            dfFNEU['w'+comp]=(wghtF)   #  weights 

            #plot
            ax[sp].plot(xF, fity)
            ax[sp].errorbar(x=xF, y=yF, yerr=yeF, 
                fmt=mfmt, ms=msz, capsize=csize, label='Final', mfc=fcol, mec=mec, mew=mew, 
                ecolor=ecol, elinewidth=elw, alpha=falpha)
            #plot labels and legend
            ax[sp].grid(axis='x', linestyle='-', color=gcol, linewidth=gvwidth, alpha=galpha)
            ax[sp].axhline(0, linestyle='-', color=gcol,linewidth=ghwidth, alpha=galpha) 
            if sp == 0 :
                ax[sp].legend(loc='upper left',fancybox=True, shadow=True)
                ax[sp].set_ylabel('North [mm]')
            if sp == 1 :
                ax[sp].set_ylabel('East [mm]')
            elif sp == 2 :
                ax[sp].set_ylabel('Vertical [mm]')
                ax[sp].set_xlabel('Date')
                [xmin,xmax,ymin,ymax]=plt.axis()
                ax[sp].text(xmin+(xmax-xmin)*.005, ymin+(ymax-ymin)*.01, 'Daily positions processed @ '+analysisCenter.upper()+". Removing mean of stations:"+" ".join(str(x) for x in meanList), 
                    horizontalalignment='left', 
                    verticalalignment='bottom')
                # put title atop as a last thing (includes values)
                ax[0].set_title('Kivu Rift Project: '+stat+'  Rates: N=%.1f±%.1f, E=%.1f±%.1f, V=%.1f±%.1f [mm/yr]' % (slopes[0], errs[0],slopes[1], errs[1],slopes[2], errs[2]))
            sp+=1
        f.savefig(os.path.join(plotdir,stat+'_TS_Common.png'), dpi=150, facecolor='white', bbox_inches='tight', pad_inches=0.5)#!/usr/bin/env python
        plt.close('all')

