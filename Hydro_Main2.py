'''
Created on Nov 20, 2014

@author: bogdan
'''
import math
import numpy
import Water_Level
import Temperature
import Hydrodynamic



# This class provides the functionality we want. You only need to look at
# this if you want to know how this works. It only needs to be defined
# once, no need to muck around with its internals.
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args:  # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False


#------------------------------------------------------------------------------------------------------
# Globals
#------------------------------------------------------------------------------------------------------
paths = {1:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/WL/csv_press_corr',
       2:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/ADCP-TorHarb',
       3:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/ADCP-TorHarb/PC-ADP/processed',
       4:'/home/bogdan/Documents/UofT/PhD/Data_Files/2010/Toberymory_tides'}

#-----------------------------------------------------------------------------------------------------

def WL_FFT_analysis_all():
    print "WL fft analysis"
    filenames = {'Emb A':'10279443_corr.csv',
                 'Emb B':'1115681_corr.csv',
                 'Emb C':'10238147_corr.csv',
                 'Cell 1':'10279696_corr.csv',
                 'Cell 2':'10279693_corr.csv',
                 'Cell 3':'10279699_corr.csv',
                 'Out Harb':'10279444_corr.csv'}


    num_segments = 10
    wla = Water_Level.WaterLevelAnalysis(paths[1], filenames, num_segments)
    for key in sorted(wla.getDict().iterkeys()):
        fname = wla.getDict()[key]

        [dates, depths] = wla.read_press_corr_file(paths[1], fname)

        # plot the original Lake oscillation input
        xlabel = 'Time [days]'
        ylabel = 'Z(t) [m]'
        ts_legend = [key + ' - water levels [m]']
        wla.plotTimeSeries("Lake levels", xlabel, ylabel, dates, depths, ts_legend)
        # end plot

        tunits = 'day'
        window = 'hanning'
        log = False
        filter = None
        [y, Time, fftx, NumUniquePts, mx, f, power, x05, x95] = \
            wla.doFFTSpectralAnalysis(dates, depths, tunits = tunits, window = window, filter = filter, log = log)

        data = []
        data.append([mx])
        data.append([key])
        data.append([x05])
        data.append([x95])
        data.append([f])

        y_label = '|Z(t)| [m]'
        title = 'Single-Sided Amplitude Spectrum vs freq'
        funits = 'cph'
        logarithmic = False
        grid = True
        plottitle = False
        ymax = None  # 0.01
        wla.plotSingleSideAplitudeSpectrumFreq(data, funits = funits, y_label = y_label, title = title, log = logarithmic, \
                                            fontsize = 20, tunits = tunits, plottitle = plottitle, grid = grid, \
                                            ymax = ymax)

        wla.plotWaveletScalogram(dates, depths, tunits, title = title)

def WL_FFT_pairs():
    num_segments = 10
    filenames = {'Emb A':'10279443_corr.csv',
                 'Emb B':'1115681_corr.csv',
                 'Emb C':'10238147_corr.csv',
                 'Cell 1':'10279696_corr.csv',
                 'Cell 2':'10279693_corr.csv',
                 'Cell 3':'10279699_corr.csv',
                 'Out Harb':'10279444_corr.csv',
                 'Inn Harb':'13320-01-MAY-30-NOV-2013_out.csv'}
    filenames = {'Out Harb':'10279444_corr.csv','Inn Harb':'13320-01-MAY-30-NOV-2013_short.csv'}
    
    for key, value in filenames.iteritems():
        fnames = [filenames['Out Harb'], value]
        names = ['Out Harb', key]
        wla = Water_Level.WaterLevelAnalysis(paths[1], fnames, num_segments)
        wla.doDualSpectralAnalysis(paths[1], fnames, names, b_wavelets = False, window = "hanning", \
                                   num_segments = num_segments, tunits = 'day', \
                                   funits = "cph", filter = None, log = False, doy = True, grid = False)

def Vel_FFT_pairs(date, plotFFT = True, skipRDI = False):

    # Process RDI-Teledyne
    num_segments = 10
    filenames = {'OH':'600mhz-DPL_002.000',
                 'EmbC':'1200mhz-EMBC_004.000'}

    data = []
    if not skipRDI :
        for key, value in filenames.iteritems():
            hyd = Hydrodynamic.Hydrodynamic(value, key, paths[2], num_segments, date, ctlname = None)
 
            hyd.readRawBinADCP()
            bins = [0, 4]
            # type can be only ampl here
            type = 'ampl'
            # type = 'power'
            if plotFFT:
                hyd.plot_FFT(key, bins, tunits = "day", funits = "cph", log = False, grid = False, type = type, withci = True)
            else:
                hyd.select_data_dates()
                data.append(hyd.getData())

    # process PC ADP
    ctlfile = 'MODE006.ctl'
    fnames = ['MODE006.ve', 'MODE006.vn', 'MODE006.vu']
    num_segments = 10
    name = 'Cell 3'
    hydpcadp = Hydrodynamic.Hydrodynamic(name, name, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
    hydpcadp.readPcAdpVel()
    bins = [0, 3]
    if plotFFT:
        hydpcadp.plot_FFT(name, bins, tunits = "day", funits = "cph", log = False, grid = False, type = type, withci = True, sel_dates = False)
    else:
        data.append(hydpcadp.getData())
    return data


def get_Dz_Dt(adptype, path, num_segments, date, dt):
     # get the water level variations dz
    print "WL analysis"
    filenames = {'OH':'10279444_corr.csv', 'EmbC':'10238147_corr.csv', 'Cell3':'10279699_corr.csv'}
    num_segments = 10
    wla = Water_Level.WaterLevelAnalysis(path, filenames, num_segments, date)
   
    [dates, depths] = wla.read_press_corr_file(path, filenames[adptype])

    print "dtime2-0=%f, dtime2-1=%f, size=%d" % (dates[0], dates[len(dates) - 1], len(dates))
    rtime, rdepths,  rdzdt, dzdt = wla.delta_z_resample(dates, depths, dt)
    #rtime2, dzdt = wla.delta_z_mov_average(dates, depths, dt)
    print "rtime2-0=%f, rtime2-1=%f, size=%d" % (rtime[0], rtime[len(rtime) - 1], len(rtime))
    return wla, rtime, rdepths, rdzdt, dates, depths, dzdt
    
def get_Velocities(adptype, date, num_segments): 
    if adptype == 'EmbC' or adptype == 'OH':
        print "Process RDI-Teledyne"
        filenames = {'EmbC':'1200mhz-EMBC_004.000', 'OH':'600mhz-DPL_002.000'}
        # hyd = Hydrodynamic.Hydrodynamic(filenames['OH'], paths[2], num_segments, date)
        hyd = Hydrodynamic.Hydrodynamic(filenames[adptype],adptype,  paths[2], num_segments, date)
        hyd.readRawBinADCP()
        hyd.select_data_dates()
    else:
        print "Process PC ADP" 
        ctlfile = 'MODE006.ctl'
        fnames = ['MODE006.ve', 'MODE006.vn', 'MODE006.vu']
        
        hyd = Hydrodynamic.Hydrodynamic(adptype, adptype, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
        hyd.readPcAdpVel()
            # these ones have dates already selected from reading. no select_data_dates is necessary
    #endif
    return hyd

def Dz_Dt_Du(date, bin, adptype = 'Cell3'):

    #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          #+ direction is TO outer harbour and to Lake Ontario    
    num_segments = 10
        
    factor = 1.  # 1 hour in days
    factor = 6.  # 10 minutes  in days

    dt = 1. / 24. / factor
    labels = [' velocity [m/s]', 'dz/dt [m/h]']

    hyd = get_Velocities(adptype, date, num_segments)
        
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)
    print "htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0]))
    # we are interested in u - along the longitudinal axis of the bay
    # sample average u values
    rtime1, udt = Hydrodynamic.Hydrodynamic.resample(hyd.time, up, dt, bin)
    # rtime1, udt = Hydrodynamic.Hydrodynamic.resample(hyd.time, hyd.results_u + 1j * hyd.results_v, dt, bin)

    #print "rtime1-0=%f, rtime1-1=%f, size=%d" % (rtime1[0], rtime1[len(rtime1) - 1], len(rtime1))
    wla, rtime2, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, date, dt)
    print "plot dz/dt - u'"
    # multiply by factor to get m/hour to match the axes labels
    wla.plot_dzdt_up_line(udt[:-3], rdzdt[:-3] * factor, labels = labels)
    wla.plot_cross_spectogram_u_dz(rtime1, udt, rtime1, numpy.array(rdzdt), scaleunit = 'hour', da = [6,300]) 
   



# Analysis of Boat Passage in FFNMP
#--------------------------------------------------------------------------------------------
def wct_lake_bp():
    date = ['10/07/25 00:00:00', '10/07/30 00:00:00']
    dt = 1. / 24. / 30.  # 10 minutes hour in days

    # get the water level variations dz
    print "WL analysis tobermory"
    filenames = {'CIH':'LL4.csv', 'IBP':'LL1.csv', 'OBP':'LL2.csv', 'HIL':'LL3.csv'}
    num_segments = 4

    wla = Water_Level.WaterLevelAnalysis(paths[4], filenames, num_segments, date)
    fname = filenames['HIL']
    [ldates, ldepths] = wla.read_press_corr_file(paths[4], fname)
    fname = filenames['CIH']
    # fname = filenames['IBP']
    # fname = filenames['OBP']
    [edates, edepths] = wla.read_press_corr_file(paths[4], fname)

    rtime1, rdepths1, rdz1, dz1 = wla.delta_z_resample(ldates, ldepths, dt)
    print "rtime1=%f, rtime1=%f, size=%d" % (rtime1[0], rtime1[len(rtime1) - 1], len(rtime1))
    rtime2, rdepths2, rdz2, dz2 = wla.delta_z_resample(edates, edepths, dt)
    print "rtime2=%f, rtime2=%f, size=%d" % (rtime2[0], rtime2[len(rtime2) - 1], len(rtime2))

    # print "plot dz/dt - u'"
    # wla.plot_dzdt_up_line(udt[:-3], dzdt[:-3])
    da = [8, 200]
    wla.plot_cross_spectogram_u_dz(rtime1, rdz1, rtime1, rdz2, scaleunit = 'hour', da = da)
#--------------------------------------------------------------------------------------------


def Vel_hodographs(date, dt, modd):
    data = Vel_FFT_pairs(date, plotFFT = False, skipRDI = False)
    tunit = 'min'
    vunit = 'm/s'
    lunit = 'm'
    bins = [[0, 2, 3], [0, 2, 4], [0, 2, 3]]
    if modd != None:
        modd = modd
    
    i = 0
    for d in data:
        [name, time, results_u, results_v] = d
        Hydrodynamic.Hydrodynamic.plotPVD(d, dt, bins[i], tunit, vunit, lunit, fontsize = 20, title = True, modd = modd)
        i += 1

def Vel_windrose(date):
    def calculateVector(results_u, results_v):
        rad = 4.*math.atan(1.0) / 180.  # degress to radians
        wspd = numpy.sqrt(results_u ** 2 + results_v ** 2)
        wdir = numpy.arctan2(results_v , results_u) / rad  # in degrees
        wdir[ wdir < 0 ] = wdir[ wdir < 0 ] + 360
        return wdir, wspd

    data = Vel_FFT_pairs(date, plotFFT = False, skipRDI = False)

    tunit = 'min'
    vunit = 'm/s'
    lunit = 'm'
    bins = [[0, 2, 3], [0, 2, 4], [0, 2, 3]]
    i = 0
    for d in data:
        [name, time, results_u, results_v] = d
        print "==============================="
        print "Location: %s" % name
        for b in bins[i]:
            print "    bin:%d" % b
            curdir, curspd = calculateVector(results_u[b], results_v[b])
            Hydrodynamic.Hydrodynamic.draw_windrose(curdir, curspd, 'bar', loc = (1, 0.05), fontsize = 20, unit = r'[$\mathsf{m\cdot s^{-1}}$]')

        i += 1

    
    
def Vel_profiles(date, adptype = 'Cell3', datetimes = None, firstbin=1, interval =1, save= False, showDZ=False):
    '''
     Based on the rotation transformation with counter clockwise 
      (+) is going to the bay 
      (-) in coming out of the bay into Outer harbour 
    
    '''
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    num_segments = 10
    
    factor = 1.  # 1 hour in days
    factor = 6.  # 10 minutes  in days

    dt = 1. / 24. / factor
    labels = [' velocity [m/s]', 'dz/dt [m/h]']

    hyd = get_Velocities(adptype, date, num_segments)

    if adptype == 'OH':
        Yrange= [0.0, 7.5]
        Xrange= [-0.70, 0.70]
    elif adptype == "EmbC":
        Yrange= [0.0, 5.5]
        Xrange= [-0.70, 0.70]
    elif adptype == "Cell3":
        Xrange= [-0.7, 0.9]
        Yrange= [0.0, 5.5]
  
    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)
    print "htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0]))
  
    #loop through the bins and plot
    profiles = []
    for date in datetimes:
        profile = []
        j = hyd.getDateIndex(date, hyd.time[0])
        for i in range(0, hyd.adcp.goodbins):  #len(up)):
            #print "date:%s i=%d, j=%d" % (date, i,j)
            profile.append(up[i][j])
        #end for
        profiles.append(profile)
        
    if showDZ:
        date = [datetimes[0], datetimes[-1]]
        dt = 1. / 24. # weeneed 1 h interval
        wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, date, dt)
        #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
        rdzdt.append(rdzdt[-1])    
        
    hyd.printVelProfiles(adptype, profiles, datetimes, firstbin, interval, Xrange, Yrange, save, numpy.array(dzdt))    

def plot_FFT_V_T_WL(adptype, date, scale = 'log', drawslope = False, resample = False):
    # 1 Get velocities
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    num_segments = 10
    
    bin = 1
    
    factor = 6.  # 10 minutes  in days
    factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    hyd = get_Velocities(adptype, date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)
    
    if resample:
        rutime, rup = hyd.resample(hyd.time, up, dt, bin)
    else:
        rup = up[bin]
        rutime= hyd.time[bin]    
    
    # 2 get Dz
    daterange = [date[0], date[-1]]
   
    wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, daterange, dt)
    #add one more to dzdt since is intepolated between nodes and profiles are on nodes.
    
    if not resample:
        rdepths = depths
        rtime = dates
    
    #3) Get Temperatures
    FORMAT="%y/%m/%d %H:%M:%S"
    
    if adptype == "OH":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-OuterHarbour/csv_processed/AboveBottom/1_Station 21"
        fname = "Bot_St21.csv"
    elif adptype == "Cell3":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap"
        fname = "Cell3.csv"
    start_num = hyd.get_date_num(daterange[0], FORMAT)
    end_num= hyd.get_date_num(daterange[1], FORMAT)
    dateTime, temp, results = hyd.get_data_from_file(fname, timeinterv = [start_num, end_num], rpath = harbour_path)
    
    
    #4) Plot 
    hyd.plot_FFT_V_T_WL(rutime, rup, dateTime, temp, rdepths, rtime, scale = 'log', drawslope = drawslope )
    
    
    
def subpl_wl_dz_vel(date, adptype, resampled = False, hourgrid = False, doy = False, img=False, cbrange=[-1,1]):
    '''
    Creates a 3 layer subplot with:
        1) water levels
        2) Dz/Dt
        3) Velocities
    '''
    
    # 1 Get velocities
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    maxdepth = {'OH':6.5, 'EmbC':5.0,'Cell3':4.5}
    
    num_segments = 10
    
    factor = 6.  # 10 minutes  in days
    factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    hyd = get_Velocities(adptype, date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)

    rup = []
    if img:
        if adptype == "OH":
            up = up[:-3]
        #=======================================================================
        # if resampled:
        #     for i in range(0,len(up)):
        #         rutime,rupi = hyd.resample(hyd.time, up, dt, i)
        #         rup.append(rupi[:-10])
        #=======================================================================
        else:
            rup = up
    else:
        if resampled:
            rutime, rup0 = hyd.resample(hyd.time, up, dt, bins[adptype][0])
            rutime, rup2 = hyd.resample(hyd.time, up, dt, bins[adptype][1])
            rutime, rup3 = hyd.resample(hyd.time, up, dt, bins[adptype][2])
            rup.append(rup0[:-10])
            rup.append(rup2[:-10])
            rup.append(rup3[:-10])
        else:
            rup.append(up[bins[adptype][0]])
            rup.append(up[bins[adptype][1]])
            rup.append(up[bins[adptype][2]])
        
    # 2 get Dz
    daterange = [date[0], date[-1]]
    wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, daterange, dt)
    #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
    
    legend = []
    legend.append("Water level") 
    legend.append("$\Delta Z/ \Delta t$")
    bin_legend = []
    for i in range(0,len(bins[adptype])):
        text = ("bin %d") % bins[adptype][i] 
        bin_legend.append(text)
    legend.append(bin_legend)  
    
    dataarr = []
    if resampled:
        dataarr.append(rdepths[:-10])
        dataarr.append(rdzdt[:-10])
        dataarr.append(rup)
        date = rtime[:-10]
    else:    
        dataarr.append(depths[:-1])
        dataarr.append(dzdt[:-1])
        dataarr.append(rup)
        date = dates[:-1]
        
    labels = ['Z [$m$]', 'dz/dt [$m$]','velocity [$m s^{-1}$]']

    hyd.display_subplots(date, hyd.time[0], dataarr, dnames = labels, yday = doy, tick = None, legend = legend, hourgrid = hourgrid, img=img, cbrange=cbrange, maxdepth = maxdepth[adptype])


if __name__ == '__main__':

    date = ['13/06/25 00:00:00', '13/08/16 00:00:00']  # this interval is needed to cover completely the intersection between velocity data and water level data
    # zoom in
    # date = ['13/07/18 00:00:00', '13/07/27 00:00:00']

    v = 'tobermory'  # for doing the XCT analysis on FFNMP data
    v = 'hodographs'
    #v = 'windrose_vel'
    v = 'dz_dt'
    v = 'subpl_wl_dz_vel'
    #v = 'vel-profiles'
    #v = 'wl_fft_all'
    #v = 'vel_fft_pairs'
    #v = 'wl_fft_pairs'
    #v = 'plot_fft_v_T_wl'
    # map the inputs to the function blocks
    for case in switch(v):
        if case('wl_fft_all'):
            WL_FFT_analysis_all()
            break
        if case('wl_fft_pairs'):
            WL_FFT_pairs()
            break
        if case('vel_fft_pairs'):
            Vel_FFT_pairs(date)
            break
        if case ('dz_dt'):
            bin = 1
            #adptype = 'OH', 'EmbC' or 'Cell3'
            adptype = 'Cell3'
            adptype = 'EmbC'
            Dz_Dt_Du(date, bin, adptype)
            break
        if case ('tobermory'):
            wct_lake_bp()
            break
        if case ('hodographs'):
            date = ['13/07/19 00:00:00', '13/07/20 00:00:00']  # for 24h progressive vector diagram - hodograph
            date = ['13/07/15 00:00:00', '13/07/16 00:00:00']
            
            # for long progressive vector diagram - hodograph
            date = ['13/07/15 00:00:00', '13/07/27 00:00:00']  
            factor = 1.  # 1 hour in days
            factor = 6.  # 10 minutes  in days
            dt = 1. / 24. / factor
            #modd = None
            modd = 24 # every nth value of the dataset gets an arrow
            Vel_hodographs(date, dt, modd)
            break
        if case ('windrose_vel'):
            Vel_windrose(date)
            break
        if  case ('vel-profiles'):
            #location can be 'Cell3' 'EmbC' 'OH'
            location = 'Cell3' 
            location = 'EmbC'
            #location = 'OH' 
            save = True
            showWL=True
            #datetimes = ['13/07/25 00:00:00','13/07/26 00:00:00','13/07/27 00:00:00','13/07/28 00:00:00','13/07/29 00:00:00']
            datetimes_storm = ['13/07/19 00:00:00','13/07/19 01:00:00','13/07/19 02:00:00','13/07/19 03:00:00','13/07/19 04:00:00','13/07/19 05:00:00',
                         '13/07/19 06:00:00','13/07/19 07:00:00','13/07/19 08:00:00','13/07/19 09:00:00','13/07/19 10:00:00','13/07/19 11:00:00',
                         '13/07/19 12:00:00','13/07/19 13:00:00','13/07/19 14:00:00','13/07/19 15:00:00','13/07/19 16:00:00','13/07/19 17:00:00',
                         '13/07/19 18:00:00','13/07/19 19:00:00','13/07/19 20:00:00','13/07/19 21:00:00','13/07/19 22:00:00','13/07/19 23:00:00',
                         
                         '13/07/20 00:00:00','13/07/20 01:00:00','13/07/20 02:00:00','13/07/20 03:00:00','13/07/20 04:00:00','13/07/20 05:00:00',
                         '13/07/20 06:00:00','13/07/20 07:00:00','13/07/20 08:00:00','13/07/20 09:00:00','13/07/20 10:00:00','13/07/20 11:00:00',
                         '13/07/20 12:00:00','13/07/20 13:00:00','13/07/20 14:00:00','13/07/20 15:00:00','13/07/20 16:00:00','13/07/20 17:00:00',
                         '13/07/20 18:00:00','13/07/20 20:00:00','13/07/20 20:00:00','13/07/20 21:00:00','13/07/20 22:00:00','13/07/20 23:00:00',
                         
                         '13/07/21 00:00:00','13/07/21 01:00:00','13/07/21 02:00:00','13/07/21 03:00:00','13/07/21 04:00:00','13/07/21 05:00:00',
                         '13/07/21 06:00:00','13/07/21 07:00:00','13/07/21 08:00:00','13/07/21 09:00:00','13/07/21 10:00:00','13/07/21 11:00:00',
                         '13/07/21 12:00:00','13/07/21 13:00:00','13/07/21 14:00:00','13/07/21 15:00:00','13/07/21 16:00:00','13/07/21 17:00:00',
                         '13/07/21 18:00:00','13/07/21 21:00:00','13/07/21 20:00:00','13/07/21 21:00:00','13/07/21 22:00:00','13/07/21 23:00:00',
                         
                         '13/07/22 00:00:00','13/07/22 01:00:00','13/07/22 02:00:00','13/07/22 03:00:00','13/07/22 04:00:00','13/07/22 05:00:00',
                         '13/07/22 06:00:00','13/07/22 07:00:00','13/07/22 08:00:00','13/07/22 09:00:00','13/07/22 10:00:00','13/07/22 11:00:00',
                         '13/07/22 12:00:00','13/07/22 13:00:00','13/07/22 14:00:00','13/07/22 15:00:00','13/07/22 16:00:00','13/07/22 17:00:00',
                         '13/07/22 18:00:00','13/07/22 22:00:00','13/07/22 20:00:00','13/07/22 21:00:00','13/07/22 22:00:00','13/07/22 23:00:00',
                         
                         '13/07/23 00:00:00','13/07/23 01:00:00','13/07/23 02:00:00','13/07/23 03:00:00','13/07/23 04:00:00','13/07/23 05:00:00',
                         '13/07/23 06:00:00','13/07/23 07:00:00','13/07/23 08:00:00','13/07/23 09:00:00','13/07/23 10:00:00','13/07/23 11:00:00',
                         '13/07/23 12:00:00','13/07/23 13:00:00','13/07/23 14:00:00','13/07/23 15:00:00','13/07/23 16:00:00','13/07/23 17:00:00',
                         '13/07/23 18:00:00','13/07/23 23:00:00','13/07/23 20:00:00','13/07/23 21:00:00','13/07/23 23:00:00','13/07/23 23:00:00',
                         ]
            
            datetimes_quiet = ['13/08/01 00:00:00','13/08/01 01:00:00','13/08/01 02:00:00','13/08/01 03:00:00','13/08/01 04:00:00','13/08/01 05:00:00',
                         '13/08/01 06:00:00','13/08/01 07:00:00','13/08/01 08:00:00','13/08/01 09:00:00','13/08/01 10:00:00','13/08/01 11:00:00',
                         '13/08/01 12:00:00','13/08/01 13:00:00','13/08/01 14:00:00','13/08/01 15:00:00','13/08/01 16:00:00','13/08/01 17:00:00',
                         '13/08/01 18:00:00','13/08/01 19:00:00','13/08/01 20:00:00','13/08/01 21:00:00','13/08/01 22:00:00','13/08/01 23:00:00',
                         
                         '13/08/02 00:00:00','13/08/02 01:00:00','13/08/02 02:00:00','13/08/02 03:00:00','13/08/02 04:00:00','13/08/02 05:00:00',
                         '13/08/02 06:00:00','13/08/02 07:00:00','13/08/02 08:00:00','13/08/02 09:00:00','13/08/02 10:00:00','13/08/02 11:00:00',
                         '13/08/02 12:00:00','13/08/02 13:00:00','13/08/02 14:00:00','13/08/02 15:00:00','13/08/02 16:00:00','13/08/02 17:00:00',
                         '13/08/02 18:00:00','13/08/02 20:00:00','13/08/02 20:00:00','13/08/02 21:00:00','13/08/02 22:00:00','13/08/02 23:00:00',
                         
                         '13/08/03 00:00:00','13/08/03 01:00:00','13/08/03 02:00:00','13/08/03 03:00:00','13/08/03 04:00:00','13/08/03 05:00:00',
                         '13/08/03 06:00:00','13/08/03 07:00:00','13/08/03 08:00:00','13/08/03 09:00:00','13/08/03 10:00:00','13/08/03 11:00:00',
                         '13/08/03 12:00:00','13/08/03 13:00:00','13/08/03 14:00:00','13/08/03 15:00:00','13/08/03 16:00:00','13/08/03 17:00:00',
                         '13/08/03 18:00:00','13/08/03 21:00:00','13/08/03 20:00:00','13/08/03 21:00:00','13/08/03 22:00:00','13/08/03 23:00:00',
                         
                         '13/08/04 00:00:00','13/08/04 01:00:00','13/08/04 02:00:00','13/08/04 03:00:00','13/08/04 04:00:00','13/08/04 05:00:00',
                         '13/08/04 06:00:00','13/08/04 07:00:00','13/08/04 08:00:00','13/08/04 09:00:00','13/08/04 10:00:00','13/08/04 11:00:00',
                         '13/08/04 12:00:00','13/08/04 13:00:00','13/08/04 14:00:00','13/08/04 15:00:00','13/08/04 16:00:00','13/08/04 17:00:00',
                         '13/08/04 18:00:00','13/08/04 22:00:00','13/08/04 20:00:00','13/08/04 21:00:00','13/08/04 22:00:00','13/08/04 23:00:00',
                         
                         '13/08/05 00:00:00','13/08/05 01:00:00','13/08/05 02:00:00','13/08/05 03:00:00','13/08/05 04:00:00','13/08/05 05:00:00',
                         '13/08/05 06:00:00','13/08/05 07:00:00','13/08/05 08:00:00','13/08/05 09:00:00','13/08/05 10:00:00','13/08/05 11:00:00',
                         '13/08/05 12:00:00','13/08/05 13:00:00','13/08/05 14:00:00','13/08/05 15:00:00','13/08/05 16:00:00','13/08/05 17:00:00',
                         '13/08/05 18:00:00','13/08/05 23:00:00','13/08/05 20:00:00','13/08/05 21:00:00','13/08/05 23:00:00','13/08/05 23:00:00',
                         ]
            
            firstbins = {"Cell3":0.3, "EmbC":1.0, "OH":1.0}
            intervals= {"Cell3":1.0, "EmbC":1.0, "OH":1.0}
            Vel_profiles(date, adptype = location, datetimes = datetimes_storm, firstbin = firstbins[location], \
                         interval = intervals[location], save = save, showDZ=showWL)
            break
        if case ('subpl_wl_dz_vel'):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            date = ['13/07/18 00:00:00', '13/07/30 00:00:00']    # shows enough of storm and quiet
            #date = ['13/07/24 00:00:00', '13/07/24 23:00:00']    # one day shows hourly alternation in DZ
            date = ['13/07/19 23:00:00', '13/07/21 13:00:00']    # stormy period
            resampled = False
            hourgrid = True
            DOY = True
            img=True
            location = "Cell3"
            subpl_wl_dz_vel(date, location, resampled, hourgrid, doy = DOY, img=img, cbrange=[None, None]) # cbrange=[-0.5,0.5])
            break
        if case ('plot_fft_v_T_wl'):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            location = "OH"
            #location ="Cell3"
            drawslope = False
            plot_FFT_V_T_WL(location, date, scale = 'log', drawslope = drawslope)
            break
        if case():  # default, could also just omit condition or 'if True'
            print ("something else!")
            # No need to break here, it'll stop anyway

    print "Done!"
