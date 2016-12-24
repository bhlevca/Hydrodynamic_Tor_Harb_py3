'''
Created on Nov 20, 2014

@author: bogdan
'''
import math
import numpy
import csv
import Water_Level
import Temperature
import Hydrodynamic
import os
import locale
import matplotlib.dates as dates
import matplotlib.mathtext
from datetime import datetime
import utools
import resample


windows = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
window_6hour = "window_6hour"  # 30 * 6 for a 2 minute sampling
window_hour = "window_hour"  # 30
window_day = "window_day"  # 30 * 24
window_half_day = "window_half_day"  # 30 * 12
window_3days = "window_3days"  # 3 * 30 * 24
window_7days = "window_7days"  # 7 * 30 * 24


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
       4:'/home/bogdan/Documents/UofT/PhD/Data_Files/2010/Toberymory_tides',
       5:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Volume-Calculations',
       6:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Carleton-Nov2013/csv_processed',
       7:'/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap',
       8:'/home/bogdan/Documents/UofT/MITACS-TRCA/3DModel/Output_data',
       9:'/home/bogdan/Documents/UofT/MITACS-TRCA/3DModel/Output_data/Model-discharge'}


def Temp_FFT_analysis_all():
    print("WL fft analysis")
    filenames = {'Stn41':'Station41November.csv',
                 }
    filenames = {'Bot_CB':                     'Bot_CB.csv',                                       
                'TC4':                    'Bot_TC4.csv',                                      
                'Cell3':                                  'Cell3.csv',                                                    
                'EmbA':            'Emb_A_Back_Water_Lagoon_1.csv',                              
                'EmbC':                'Emb_C_West_Embayment.csv',                                  
                'SpaDock':                    'SpaDock.csv',                                      
                'Stn1':                             'Stn1.csv',                               
                'Stn27':                            'Stn27.csv',                              
                'Stn43':                            'Stn43.csv',                              
                'Stn61':                      'Stn61.csv',                                        
                'TC1':                    'Bot_TC1.csv',                                      
                'CB':                         'CB.csv',                                           
                'CIFerry':                    'Center_Island_Ferry.csv',                                      
                'EmbB':                 'Emb_B_Inside_Wetland.csv',                                   
                'JarDock':                             'JarDock.csv',                                               
                'Stn13':                      'Stn13.csv',                                        
                'Stn20':                            'Stn20.csv',                              
                #'Stn40':                            'Stn40.csv',                              
                'Stn44':                            'Stn44.csv',                              
                'Stn8':                      'Stn8.csv',                                        
                'TC2':                    'Bot_TC2.csv',                                      
                'Cell1':                    'Cell_1A.csv',                                      
                'FerryBT':            'East_of_Ferry_Boat_Terminal.csv',                              
                 #'EmbC':                                 'EmbC.csv',                                                   
                 #'MidDock':                             'MidDock.csv',                                               
                'Stn14':         'Stn14.csv',                                        
                'Stn21':         'Stn21.csv',                              
                'Stn41':         'Stn41.csv',                              
                'Stn49':         'Stn49.csv',                              
                'TIBBM':         'Toronto_Islands_Blockhouse_Bay_Mouth.csv',                 
                'TC3':           'Bot_TC3.csv',                                      
                'Cell2':         'Cell2.csv',                                        
                 #'EGap':        'EGap.csv',                                        
                'EmbCOB':        'Emb_C_Outside_Berm.csv',                                     
                'SimcoeDock':    'Simcoe_Slip.csv',                                           
                'Stn19':         'Stn19.csv',                                        
                'Stn22':         'Stn22.csv',                              
                #'Stn42':        'Stn42.csv',                              
                'Stn50':         'Stn50.csv',                              
                'WGap':          'WGap.csv',                                                     
                 
                 }
    

    #tpath = paths[6]
    tpath = paths[7]
    num_segments = 10
    ta = Temperature.Temperature(tpath, filenames, num_segments)
    for key in sorted(ta.getDict().keys()):
        fname = ta.getDict()[key]
        print("File: %s" % fname)
        [dates, depths] = ta.read_temp_file(tpath, fname)
        
        #interpoate if time samolungf is > 10 min
        dt = (dates[2]-dates[1])*86400
        if  dt > 600:  #10 mio
            newdates = numpy.linspace(dates[0], dates[-1], int(len(dates)*dt/600))
            depths = numpy.interp(newdates, dates, depths)
            dates = newdates

        # plot the original Lake oscillation input
        #xlabel = 'Time (days)'
        #ylabel = 'Z[t] (m)'
        #ts_legend = [key + ' - water levels [m]']
        #ta.plotTimeSeries("Lake levels", xlabel, ylabel, dates, depths, ts_legend)
        # end plot

        tunits = 'day'
        window = 'hanning'
        log = 'loglog'
        filter = None
        [y, Time, fftx, NumUniquePts, mx, f, power, x05, x95] = \
            ta.doFFTSpectralAnalysis(dates, depths, tunits = tunits, window = window, filter = filter, log = log)

        data = []
        data.append([mx])
        data.append([key])
        data.append([x05])
        data.append([x95])
        data.append([f])

        y_label = 'Temperature [$^\circ$C]'
        title = 'Single-Sided Amplitude Spectrum vs freq'
        funits = 'cph'
        logarithmic = 'loglog'
        grid = False
        plottitle = False
        ymax = None  # 0.01
        ta.plotSingleSideAplitudeSpectrumFreq(data, funits = funits, y_label = y_label, title = title, log = logarithmic, \
                                            fontsize = 20, tunits = tunits, plottitle = plottitle, grid = grid, \
                                            ymax = ymax)

        #ta.plotWaveletScalogram(dates, depths, tunits, title = title)

def WL_FFT_analysis_all():
    print("WL fft analysis")
    filenames = {'Emb A':'10279443_corr.csv',
                 'Emb B':'1115681_corr.csv',
                 'Emb C':'10238147_corr.csv',
                 'Cell 1':'10279696_corr.csv',
                 'Cell 2':'10279693_corr.csv',
                 'Cell 3':'10279699_corr.csv',
                 'Out Harb':'10279444_corr.csv'}

    wav = False 
    num_segments = 10
    wla = Water_Level.WaterLevelAnalysis(paths[1], filenames, num_segments)
    for key in sorted(wla.getDict().keys()):
        fname = wla.getDict()[key]

        [dates, depths] = wla.read_press_corr_file(paths[1], fname)

        # plot the original Lake oscillation input
        xlabel = 'Time [days]'
        # ylabel = 'Z(t) [m]'
        ylabel = '|Z(f)| [m]'
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

        title = 'Single-Sided Amplitude Spectrum vs freq'
        funits = 'cph'
        logarithmic = 'loglog'
        grid = False
        plottitle = False
        ymax = None  # 0.01
        wla.plotSingleSideAplitudeSpectrumFreq(data, funits = funits, y_label = ylabel, title = title, log = logarithmic, \
                                            fontsize = 18, tunits = tunits, plottitle = plottitle, grid = grid, \
                                            ymax = ymax)

        if wav:
            wla.plotWaveletScalogram(dates, depths, tunits, title = title)

def WL_FFT_pairs(dateinterval=None):
    num_segments = 10
    filenames = {'Emb A':'10279443_corr.csv',
                 'Emb B':'1115681_corr.csv',
                 'Emb C':'10238147_corr.csv',
                 'Cell 1':'10279696_corr.csv',
                 'Cell 2':'10279693_corr.csv',
                 'Cell 3':'10279699_corr.csv',
                 'Out Harb':'10279444_corr.csv',
                 'IH data':'13320-01-AUG-13-AUG-2013.csv',
                 'TC4':'01-13-Aug-WL-TC4_out.csv',
                 'Lake Ontario':'10279444_corr.csv',
                 # ---------------------------
                 'IH model':'01-13-Aug-Water_level_Spadinapai.csv',
                 'IH model T=5h':'Model-HWL-Spadina-T=5h_out.csv',
                 'OH model T=5h':'Model-HWL-TC4-T=5h_out.csv',
                 'IH model T=1h':'Model-HWL-Spadina-T=1h_out.csv',
                 'OH model T=1h':'Model-HWL-TC4-T=1h_out.csv',
                 'OH model':'Model-WL-OHSt18-5-7Aug_out.csv',
                 'IH model':'Model-WL-Spadina-5-7Aug_out.csv',
                 'OH model':'Model-WL-TC4-5-7Aug_out.csv',
                 'Cell 3':'Cell3-cross_section_discharge-7-9Aug_out.csv',
                 'Emb C':'EmbC-cross_section_discharge-7-9Aug_out.csv',
                 'E Gap':'EG-cross_section_discharge-7-9Aug_out.csv',
                 'IH model':'Model-HWL-Spadina-T=0.33h_out.csv',
                 'OH model':'Model-HWL-TC4-T=0.33h_out.csv',
                 'IH model':'Model-HWL-Spadina-T=1.7h_out.csv',
                 'OH model':'Model-HWL-TC4-T=1.7h_out.csv',
                 'IH model':'Model-HWL-Spadina-T=3.2h_out.csv',
                 'OH model': 'Model-HWL-TC4-T=3.2h_out.csv'}
                 
    filenames = {'Inn Harb':'13320-01-MAY-30-NOV-2013_short.csv', 'Out Harb':'10279444_corr.csv'}
    #filenames = {'Emb A':'10279443_corr.csv', 'Lake Ontario':'10279444_corr.csv'}
    #filenames = {'TC4':'01-13-Aug-Water_level_TC4_out.csv','Spadina':'01-13-Aug-Water_level_Spadina_out.csv'}
    #Original output with input from Spadina real data
    #filenames = {'IH model':'01-13-Aug-WL-Spadina_out.csv','IH data':'13320-01-AUG-13-AUG-2013.csv'}
    #filenames = {'IH data':'13320-01-AUG-13-AUG-2013.csv', 'IH model':'01-13-Aug-WL-Spadina_out.csv'}
    #Output with input from Spadina real data divided by 3 +shift 36 min
    #filenames = {'IH data':'13320-01-AUG-13-AUG-2013.csv', 'IH model':'01-13-Aug-Water_level_Spadina_out.csv'}
    #filenames = {'IH data':'13320-01-AUG-13-AUG-2013.csv', 'IH model':'Model-Spadina-Water level-2Aug2013_out.csv'}
    #Output with input from Spadina real data divided by 3 and increased Manning friction (double to 0.046)  +shift 36 min
    #filenames = {'IH data':'13320-01-AUG-13-AUG-2013.csv', 'IH model':'Model-Spadina-Water level-2Aug2013_out-fric0.046.csv'}
    #Output with input from Spadina real data divided by 7 and increased Manning friction (double to 0.035)  +shift 36 min
    #filenames = {'IH data':'13320-01-AUG-13-AUG-2013.csv', 'IH model':'Model-Spadina-Water level-2Aug2013_out_X7.csv'}
    #filenames = {'IH model T=5h':'Model-HWL-Spadina-T=5h_out.csv', 'OH model T=5h':'Model-HWL-TC4-T=5h_out.csv'}
    #filenames = {'IH model T=1h':'Model-HWL-Spadina-T=1h_out.csv', 'OH model T=1h':'Model-HWL-TC4-T=1h_out.csv'}
    #filenames = {'IH model':'Model-WL-Spadina-5-7Aug_out.csv', 'IH data':'Spadina-AUG-2013.csv'}
    #filenames = {'OH model':'Model-WL-TC4-5-7Aug_out.csv', 'OH data':'St18_corr.csv'}
    #filenames = {'IH model T=3.2h':'Model-HWL-Spadina-T=3.2h_out.csv', 'OH model T=3.2h': 'Model-HWL-TC4-T=3.2h_out.csv'}
    #filenames = {'IH model T=1.7h':'Model-HWL-Spadina-T=1.7h_out.csv', 'OH model T=1.7h': 'Model-HWL-TC4-T=1.7h_out.csv'}
    #filenames = {'IH model T=0.33h':'Model-HWL-Spadina-T=0.33h_out.csv', 'OH model T=0.33h': 'Model-HWL-TC4-T=0.33h_out.csv'}
    #filenames = {'IH model T=1.0h':'Model-HWL-Spadina-T=1.0h_out.csv', 'OH model T=1.0h': 'Model-HWL-TC4-T=1.0h_out.csv'}
    
    #filenames = {'Emb C':'EmbC-cross_section_discharge-7-9Aug_out.csv','E Gap':'EG-cross_section_discharge-7-9Aug_out.csv'}
    path_no=1
    #path_no=9
    #path_no=8
    for key, value in filenames.items():
        
        names = ['Out Harb', key]
        fnames = [filenames['Inn Harb'], filenames['Out Harb']]
        #names = ['Emb A', key]
        #fnames = [filenames['Emb A'], filenames['Lake Ontario']]
        #fnames = [filenames['IH model'], value]
        #fnames = [filenames['OH model'], filenames['OH data']]
        #fnames = [filenames['Emb C'], filenames['E Gap']]
        #fnames = [filenames['IH model T=5h'], filenames['OH model T=5h']]
        #fnames = [filenames['IH model T=3.2h'], filenames['OH model T=3.2h']]
        #fnames = [filenames['IH model T=1.7h'], filenames['OH model T=1.7h']]
        #fnames = [filenames['IH model T=0.33h'], filenames['OH model T=0.33h']]
        #fnames = [filenames['IH model T=1.0h'], filenames['OH model T=1.0h']]
        #names = ['Out Harb', key]
        #names = ['TC4', key]
        #names = ['IH model T=5h', 'OH model T=5h']
        #names = ['Emb C', 'E GAp']
        #names = ['OH model', 'OH data']
        #names = ['IH model T=3.2h', 'OH model T=3.2h']
        #names = ['IH model T=1.7h', 'OH model T=1.7h']
        #names = ['IH model T=0.33h', 'OH model T=0.33h']
        #names = ['IH model T=1.0h', 'OH model T=1.0h']

        log = 'loglog'
        #log = False
        wla = Water_Level.WaterLevelAnalysis(paths[path_no], fnames, num_segments)
        wla.doDualSpectralAnalysis(paths[path_no], fnames, names, b_wavelets=False, window="hanning",
                                   num_segments=num_segments, tunits = 'day',
                                   funits="cph", filter=None, log=log, doy=True, grid=False,
                                   dateinterval=dateinterval, d3d=False)

def Vel_FFT_pairs(date, plotFFT = True, skipRDI = False):

    # Process RDI-Teledyne
    num_segments = 3
    filenames = {'OH':'600mhz-DPL_002.000',
                 'EmbC':'1200mhz-EMBC_004.000'}

    log = 'loglog'
     # type can be only ampl here
    type = 'ampl'
    data = []
    if not skipRDI :
        for key, value in filenames.items():
            hyd = Hydrodynamic.Hydrodynamic(value, key, paths[2], num_segments, date, ctlname = None)
 
            hyd.readRawBinADCP()
            bins = [0, 4]
           
            
            # type = 'power'
            if plotFFT:
                hyd.plot_FFT(key, bins, tunits = "day", funits = "cph", log = log, grid = False, type = type, withci = True)
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
        hydpcadp.plot_FFT(name, bins, tunits = "day", funits = "cph", log = log, grid = False, type = type, withci = True, sel_dates = False)
    else:
        data.append(hydpcadp.getData())
    return data


def Temp_FFT(date, plotFFT = True, skipRDI = False):

    # Process RDI-Teledyne
    num_segments = 10
    filenames = {'OH':'600mhz-DPL_002.000',
                 'EmbC':'1200mhz-EMBC_004.000'}

    log = 'loglog'
     # type can be only ampl here
    type = 'ampl'
    data = []
    if not skipRDI :
        for key, value in filenames.items():
            hyd = Hydrodynamic.Hydrodynamic(value, key, paths[2], num_segments, date, ctlname = None)
 
            hyd.readRawBinADCP()
            # type = 'power'
            if plotFFT:
                hyd.plot_Temp_FFT(key, tunits="day", funits="cph", log=log, grid=False, type=type, withci=True)
            else:
                hyd.select_data_dates()
                data.append(hyd.getData())

    # process PC ADP
    ctlfile = 'MODE006.ctl'
    fnames = ['MODE006.hdr']
    num_segments = 10
    name = 'Cell 3'
    hydpcadp = Hydrodynamic.Hydrodynamic(name, name, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
    time, temp = hydpcadp.readPcAdpTemp()
    if plotFFT:
        hydpcadp.plot_Temp_FFT(name, tunits = "day", funits = "cph", log = log, grid = False, type = type, withci = True, sel_dates = False)
    else:
        data.append(hydpcadp.getData())
    return data



def get_Dz_Dt(adptype, path, num_segments, date, dt):
     # get the water level variations dz
    print("WL analysis")
    filenames = {'OH':'10279444_corr.csv', 'EmbC':'10238147_corr.csv', 'Cell3':'10279699_corr.csv', \
                 'Cell1':'10279696_corr.csv', 'Cell2':'10279693_corr.csv', \
                 'EmbA':'10279443_corr.csv', 'EmbB':'1115681_corr.csv'}
    num_segments = 10
    wla = Water_Level.WaterLevelAnalysis(path, filenames, num_segments, date)
   
    [dates, depths] = wla.read_press_corr_file(path, filenames[adptype])

    print("dtime2-0=%f, dtime2-1=%f, size=%d" % (dates[0], dates[len(dates) - 1], len(dates)))
    rtime, rdepths,  rdzdt, dzdt = wla.delta_z_resample(dates, depths, dt)
    #rtime2, dzdt = wla.delta_z_mov_average(dates, depths, dt)
    print("rtime2-0=%f, rtime2-1=%f, size=%d" % (rtime[0], rtime[len(rtime) - 1], len(rtime)))
    return wla, rtime, rdepths, rdzdt, dates, depths, dzdt
    
def get_Velocities(adptype, date, num_segments, delft3d=False, d3d_int=3): 
    if adptype == 'EmbC' or adptype == 'OH':
        print("Process RDI-Teledyne")
        filenames = {'EmbC':'1200mhz-EMBC_004.000', 'OH':'600mhz-DPL_002.000'}
        # hyd = Hydrodynamic.Hydrodynamic(filenames['OH'], paths[2], num_segments, date)
        hyd = Hydrodynamic.Hydrodynamic(filenames[adptype],adptype,  paths[2], num_segments, date)
        hyd.readRawBinADCP()
        hyd.select_data_dates(delft3d=delft3d, d3d_int=d3d_int)
    else:
        print("Process PC ADP") 
        ctlfile = 'MODE006.ctl'
        fnames = ['MODE006.ve', 'MODE006.vn', 'MODE006.vu']
        
        hyd = Hydrodynamic.Hydrodynamic(adptype, adptype, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
        hyd.readPcAdpVel()
        hyd.select_data_dates(delft3d=delft3d, d3d_int=d3d_int)
            # these ones have dates already selected from reading. no select_data_dates is necessary
    #endif
    return hyd

def Dz_Dt_Du(date, bin, adptype = 'Cell3'):

    #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise   and vp the along dir
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
    print("htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0])))
    # we are interested in u - along the longitudinal axis of the bay
    # sample average u values
    rtime1, udt = Hydrodynamic.Hydrodynamic.resample(hyd.time, vp, dt, bin)
    # rtime1, udt = Hydrodynamic.Hydrodynamic.resample(hyd.time, hyd.results_u + 1j * hyd.results_v, dt, bin)

    #print "rtime1-0=%f, rtime1-1=%f, size=%d" % (rtime1[0], rtime1[len(rtime1) - 1], len(rtime1))
    wla, rtime2, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, date, dt)
    print("plot dz/dt - u'")
    # multiply by factor to get m/hour to match the axes labels
    wla.plot_dzdt_up_line(udt[:-3], rdzdt[:-3] * factor, labels = labels)
    wla.plot_cross_spectogram_u_dz(rtime1, udt, rtime1, numpy.array(rdzdt), scaleunit = 'hour', da = [6,300]) 
   

def wct_V_T(date, bin, adptype = 'Cell3'):
    #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  and vp the along dir
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
    print("htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0])))
    # we are interested in u - along the longitudinal axis of the bay
    # sample average u values
    rtime1, udt = Hydrodynamic.Hydrodynamic.resample(hyd.time, vp, dt, bin)
    #print "rtime1-0=%f, rtime1-1=%f, size=%d" % (rtime1[0], rtime1[len(rtime1) - 1], len(rtime1))
    wla, rtime2, rdepths, rdzdt, times, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, date, dt)
    
    dt = datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
    start_num = dates.date2num(dt)
    dt = datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
    end_num = dates.date2num(dt)
    
    #3) Get Temperatures
    if adptype == "OH":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-OuterHarbour/csv_processed/AboveBottom/1_Station_21"
        fname = "Bot_St21.csv"
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-OuterHarbour/csv_processed/AboveBottom/2_Temperature Chain 3"
        fname = "TC3_1m.csv"
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-OuterHarbour/csv_processed/AboveBottom/0_Temperature Chain 4"
        fname = "TC4_1m.csv"
        
        dateTime, temp, results = hyd.get_data_from_file(fname, timeinterv = [start_num, end_num], rpath = harbour_path)
    elif adptype == "Cell3":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap"
        fname = "Cell3.csv"
        ctlfile = 'MODE006.ctl'
        fnames = ['MODE006.hdr']
        name = 'Cell 3'
        num_segments = 10
        hydpcadp = Hydrodynamic.Hydrodynamic(name, name, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
        dateTime, temp = hydpcadp.readPcAdpTemp()
    elif adptype == "EmbC":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap"
        fname = "EmbC.csv"
        fname = "Emb_C_West_Embayment.csv"
        dateTime, temp, results = hyd.get_data_from_file(fname, timeinterv = [start_num, end_num], rpath = harbour_path)
    

    wla.plot_cross_spectogram_u_dz(rtime1, udt, dateTime, numpy.array(temp), scaleunit = 'hour', da = [6,300])
    

def wct_lake_bp():
    # Analysis of Boat Passage in FFNMP
    date = ['10/07/25 00:00:00', '10/07/30 00:00:00']
    dt = 1. / 24. / 30.  # 10 minutes hour in days

    # get the water level variations dz
    print("WL analysis tobermory")
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
    print("rtime1=%f, rtime1=%f, size=%d" % (rtime1[0], rtime1[len(rtime1) - 1], len(rtime1)))
    rtime2, rdepths2, rdz2, dz2 = wla.delta_z_resample(edates, edepths, dt)
    print("rtime2=%f, rtime2=%f, size=%d" % (rtime2[0], rtime2[len(rtime2) - 1], len(rtime2)))

    # print "plot dz/dt - u'"
    # wla.plot_dzdt_up_line(udt[:-3], dzdt[:-3])
    da = [8, 200]
    wla.plot_cross_spectogram_u_dz(rtime1, rdz1, rtime1, rdz2, scaleunit = 'hour', da = da)

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

def Vel_windrose(date, skipRDI = False):
    def calculateVector(results_u, results_v):
        rad = 4. * math.atan(1.0) / 180.  # degress to radians
        wspd = numpy.sqrt(results_u ** 2 + results_v ** 2)
        wdir = numpy.arctan2(results_v , results_u) / rad  # in degrees
        wdir[ wdir < 0 ] = wdir[ wdir < 0 ] + 360
        return wdir, wspd

    data = Vel_FFT_pairs(date, plotFFT = False, skipRDI = skipRDI)

    tunit = 'min'
    vunit = 'm/s'
    lunit = 'm'
    bins = [[0, 2, 4], [0, 2, 5], [0, 2, 4]]
    #bins = [[3], [4], [3]]
    i = 0
    for d in data:
        [name, time, results_u, results_v] = d
        print("===============================")
        print("Location: %s" % name)
        for b in bins[i]:
            print("    bin:%d" % b)
            curdir, curspd = calculateVector(results_u[b], results_v[b])
            Hydrodynamic.Hydrodynamic.draw_windrose(curdir, curspd, 'bar', loc = (1, 0.05), fontsize = 20,\
                                                    unit = '[m s$^{-1}}$]', angle = 67.5) 

        i += 1

    
    
def Vel_profiles(date, adptype = 'Cell3', datetimes = None, firstbin=1, interval =1, save= False, showDZ=False):
    '''
     Based on the rotation transformation with counter clockwise 
      (+) is going to the bay 
      (-) in coming out of the bay into Outer harbour 
    
    '''
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise and vp the along dir 
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
    print("htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0])))
  
    #loop through the bins and plot
    profiles = []
    for date in datetimes:
        profile = []
        j = hyd.getDateIndex(date, hyd.time[0])
        for i in range(0, hyd.goodbins):  #len(up)):
            #print "date:%s i=%d, j=%d" % (date, i,j)
            profile.append(vp[i][j])
        #end for
        profiles.append(profile)
        
    if showDZ:
        date = [datetimes[0], datetimes[-1]]
        dt = 1. / 24. # weeneed 1 h interval
        wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, date, dt)
        #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
        rdzdt.append(rdzdt[-1])    
        
    hyd.printVelProfiles(adptype, profiles, datetimes, firstbin, interval, Xrange, Yrange, save, numpy.array(dzdt))    


def calc_eddy_visc(hyd, adptype, dateint, num_segments, dz, angles_from_N):
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    
    revert = True # True 0 bin bottom of graph because it is reversing limits in the routine below 
    
    if adptype != "Cell3":
        [time1,  tresults_u, tresults_v, tresults_z] = hyd.select_data_dates(dateint = dateint, velprofile = False)
        up, vp = hyd.rotation_transform(tet, clockwise=False, u=tresults_u, v=tresults_v)
        diff = hyd.goodbins - len(vp)
        if diff ==0 :
            timearr = time1[:]
            velarr = vp[:]
            welarr = tresults_z[:]
        else:
            timearr = time1[:diff]
            velarr = vp[:diff]
            welarr = tresults_z[:diff]
        return hyd.eddy_viscosity(velarr, welarr, dz)
    else:
        hyd2 = get_Velocities(adptype, dateint, num_segments)
        #reproject on the direction of thde flow
        up, vp = hyd2.rotation_transform(tet, clockwise=False)
        diff = hyd2.goodbins - len(vp)        
        if diff ==0 :
            timearr = hyd2.time[:]
            velarr = vp[:]
            welarr = hyd2.results_z[:]
        else:
            timearr = hyd2.time[:diff]
            velarr = vp[:diff]
            welarr = hyd2.results_z[:diff]
            
        return hyd2.eddy_viscosity(velarr, welarr, dz)
        
    

def Avg_Vel_profiles(date, adptype = 'Cell3', firstbin=1, interval =1,\
                     eddyvisc = False, ADCPprof = False, grad_valong = False):
    '''
     Based on the rotation transformation with counter clockwise 
      (+) is going to the bay 
      (-) in coming out of the bay into Outer harbour 
    
    '''
    matplotlib.mathtext.SHRINK_FACTOR = 0.9


    # angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  and vp the along dir
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    num_segments = 10
    
    factor = 1.  # 1 hour in days
    factor = 6.  # 10 minutes  in days

    dt = 1. / 24. / factor
    labels = [' velocity [m/s]']

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
    
    #if counterclockwise vp is the along the channel veloccity and + (PLUS) is toward the lake 
    #if clockwise up is the along the channel veloccity and + (PLUS) is toward the embayment 
    print("htime-0=%f, htime0-1=%f, size=%d" % (hyd.time[0][0], hyd.time[0][len(hyd.time[0]) - 1], len(hyd.time[0])))
  
    upw_str = '13/07/20 12:00:00'
    dt = datetime.strptime(upw_str, "%y/%m/%d %H:%M:%S")
    upwell = dates.date2num(dt)
    
    dw_str = '13/07/26 12:00:00'
    dt = datetime.strptime(dw_str, "%y/%m/%d %H:%M:%S")
    downwell = dates.date2num(dt)
    
    stab_str = '13/08/05 12:00:00'
    dt = datetime.strptime(stab_str, "%y/%m/%d %H:%M:%S")
    stable = dates.date2num(dt)
    profiledates = [stable, upwell, downwell]
    
    
    #profiledates = None
    revert = True # True 0 bin bottom of graph because it is reversing limits in the routine below 
    
    diff = hyd.goodbins - len(vp)
    if diff ==0 :
        timearr = hyd.time[:]
        velarr = vp[:]
    else:
        timearr = hyd.time[:diff]
        velarr = vp[:diff]
    utools.display_data.display_avg_vertical_temperature_profiles_err_bar_range([timearr], [velarr], startdeptharr = [0], \
                                                                                profiledates = profiledates, doy = True, \
                                                                                revert = revert, \
                                                                                legendloc = 4, xlabel = ' Velocity [$m s^{-1}$]',\
                                                                                errbar= True, rangebar= True,\
                                                                                debug = False)  
    #calculate the eddy viscosity
    
    if eddyvisc:
        upw_end_str = '13/07/21 12:00:00'
        dw_end_str = '13/07/27 12:00:00'
        stab_end_str = '13/08/06 12:00:00'
        datesarr = [[stab_str, stab_end_str], [upw_str, upw_end_str],[dw_str, dw_end_str]]
        
        dz =1
        Kv = []
        for d in datesarr:
            Kv.append(calc_eddy_visc(hyd, adptype, d, num_segments, dz, angles_from_N))
         
        startdepth = 0  
        legend= ['DOY 218', 'DOY 201', 'DOY 207'] 
        utools.display_data.display_vertical_velocity_profiles(numpy.array(Kv)*0.5, startdepth, revert = True, legendloc = 4, legend = legend,\
                                                               grid = False, xlabel = r'$K_\mathsf{\nu}$ [$\mathsf{m^2s^{-1}}$]', title = None)
        
    
    if grad_valong:
        # 3) Mixed water, air ,img data
        # ToDO: Add short and long radiation
        print("Grad Valong Start display mixed subplots  ")
        grad = []
        V_T  = numpy.transpose(vp) 
        U_T  = numpy.transpose(up)
               
        for i in range(0, len(hyd.time[0])):
            dvdz = []
            for j in range(0, len(hyd.time)-1):
                W1 = math.sqrt(V_T[i][j+1]** 2 + U_T[i][j+1]** 2)
                W2 = math.sqrt(V_T[i][j] ** 2 + U_T[i][j]** 2)
                dvdz.append((W1-W2)/dz/2)
            #edn for 
            dvdz.insert(0, (0+dvdz[0])/2.)    #insert fist item to match the velocity array size
            grad.append(dvdz)
            #grad.append(numpy.gradient(, dz))
        #end for
        
        
        gradarr = numpy.array(grad)
        grad_T = numpy.transpose(gradarr)

        if adptype == "Cell3":
            dd = 1 # for Cell 3
        elif adptype == "EmbC":
            dd=2
        elif adptype == "OH":
            dd = 3   
        dateTimes3 = [hyd.time[:-dd], hyd.time[:-dd], hyd.time[:-dd]]
        ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]"] 
        
        imgs = [vp[:-dd], up[:-dd], grad_T[:-dd]]
        #imgs = [hyd.results_u[:-dd], hyd.results_v[:-dd]] #Harbour first
        
        if adptype == 'OH':
            t11 = ['0', '2', '4.5', '7']
            t12 = [7, 4.5, 2, 0]
            
            t21 = ['0', '2', '4.5', '7']
            t22 = [7, 4.5, 2, 0]
            
            t31 = ['0', '2', '4.5', '7']
            t32 = [7, 4.5, 2, 0]
        
            #maxdepth = [9, 27] # Harbour first
            maxdepth = [7, 7, 7]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0, 0]
            mindepths = [0, 0, 0]
            maxtemp = [0.2, 0.2, 0.05]
            mintemps = [-0.2, -0.2, -0.05]
        elif adptype == "EmbC" or adptype == "Cell3":
            t11 = ['0', '1.5', '3.0', '4.5']
            t12 = [4.5, 3, 1.5, 0]
            
            t21 = ['0', '1.5', '3.0', '4.5']
            t22 = [4.5, 3, 1.5, 0]
            
            t31 = ['0', '1.5', '3.0', '4.5']
            t32 = [4.5, 3, 1.5, 0]
        
            #maxdepth = [9, 27] # Harbour first
            maxdepth = [4.5, 4.5, 4.5]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0, 0]
            mindepths = [0, 0, 0]
            maxtemp = [0.2, 0.2, 0.05]
            mintemps = [-0.2, -0.2, -0.05]
            
       
        
        tick = [[t11, t12], [t21, t22], [t31, t32]]
        
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
            
        
        utools.display_data.display_mixed_subplot(dateTimes1=[], data=[], varnames=[], ylabels1=[],
                                           dateTimes2=[], groups=[], groupnames=[], ylabels2=[],
                                           dateTimes3=dateTimes3, imgs=imgs, ylabels3=ylabels3, ticks=tick,
                                           maxdepths=maxdepth, mindepths=mindepths, mintemps=mintemps,
                                           firstlogs=firstlogdepth, maxtemps=maxtemp,
                                           fnames=None, revert=True, custom=None, maxdepth=None, tick=None,
                                           firstlog=None, yday=True, title=False, grid=False, limits=limits,
                                           sharex=True, fontsize=18, group_first=False, interp=2,
                                           cblabel=[r'$\mathsf{U_{alg}}$ [$\mathsf{m\ s^{-1}}$]',
                                                    r'$\mathsf{V_{crs}}$ [$\mathsf{m\ s^{-1}}$]',
                                                    r'$\mathsf{dV\ dz^{-1}}$ [$\mathsf{s^{-1}}$]'])

    if ADCPprof:
        # 3) Mixed water, air ,img data
        # ToDO: Add short and long radiation
        print("Start display mixed subplots ")
        dd = 1
        dateTimes3 = [hyd.time[:-dd], hyd.time[:-dd], hyd.time[:-dd]]
        ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]"] 
        #imgs = [TH_resultsArr, resultsArr] #Harbour first
        imgs = [ hyd.results_u[:-dd], hyd.results_v[:-dd], hyd.results_z[:-dd]]
        
        
        t11 = ['0', '2', '4.5', '7']
        t12 = [7, 4.5, 2, 0]
        
        t21 = ['0', '2', '4.5', '7']
        t22 = [7, 4.5, 2, 0]
        
        t31 = ['0', '2', '4.5', '7']
        t32 = [7, 4.5, 2, 0]
        
        tick = [[t11, t12], [t21, t22], [t31, t32]]
        
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
            
        #maxdepth = [9, 27] # Harbour first
        maxdepth = [7, 7, 7]
        #firstlogdepth = [0, 3] Harbour first
        firstlogdepth = [0, 0, 0]
        
        maxtemp = [0.4, 0.4, 0.02]
        mintemps = [-0.4, -0.4, -0.02]
        mindepths = [0, 0, 0]
        
        utools.display_data.display_mixed_subplot(dateTimes1 = [], data = [], varnames = [], ylabels1 = [], \
                                           dateTimes2 = [], groups = [], groupnames = [], ylabels2 = [], \
                                           dateTimes3 = dateTimes3, imgs = imgs, ylabels3 = ylabels3, ticks = tick, maxdepths = maxdepth, \
                                            mindepths = mindepths, mintemps = mintemps, firstlogs = firstlogdepth, maxtemps = maxtemp, \
                              fnames = None, revert = True, custom = None, maxdepth = None, tick = None, firstlog = None, yday = True, \
                              title = False, grid = False, limits = limits, sharex = True, fontsize = 18, group_first = False, interp = 2,\
                              cblabel =  ["$V_E$ [$m s^{-1}$]", "$V_N$ [$m s^{-1}$]","$V_U$ [$m s^{-1}$]"])


def plot_FFT_V_T_WL(adptype, date, scale = 'log', drawslope = False, resample = False):
    # 1 Get velocities
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    num_segments = 10
    
    bin = 0
    
    factor = 6.  # 10 minutes  in days
    factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    hyd = get_Velocities(adptype, date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)
    
    if resample:
        rvtime, rvp = hyd.resample(hyd.time, vp, dt, bin)
    else:
        rvp = vp[bin]
        rvtime= hyd.time[bin]    
    
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
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-OuterHarbour/csv_processed/AboveBottom/1_Station_21"
        fname = "Bot_St21.csv"
        start_num = hyd.get_date_num(daterange[0], FORMAT)
        end_num= hyd.get_date_num(daterange[1], FORMAT)
        dateTime, temp, results = hyd.get_data_from_file(fname, timeinterv = [start_num, end_num], rpath = harbour_path)
    elif adptype == "Cell3":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap"
        fname = "Cell3.csv"
        ctlfile = 'MODE006.ctl'
        fnames = ['MODE006.hdr']
        name = 'Cell 3'
        num_segments = 10
        hydpcadp = Hydrodynamic.Hydrodynamic(name, name, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
        dateTime, temp = hydpcadp.readPcAdpTemp()
    elif adptype == "EmbC":
        harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/ClimateMap"
        fname = "EmbC.csv"
        fname = "Emb_C_West_Embayment.csv"
        # or 
        #harbour_path = "/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/AllHarbour/csv_processed"
        #fname = "1020950.csv"
        
        start_num = hyd.get_date_num(daterange[0], FORMAT)
        end_num= hyd.get_date_num(daterange[1], FORMAT)
        dateTime, temp, results = hyd.get_data_from_file(fname, timeinterv = [start_num, end_num], rpath = harbour_path)
    
    
    #4) Plot 
    hyd.plot_FFT_V_T_WL(rvtime, rvp, dateTime, temp, rdepths, rtime, scale = 'loglog', drawslope = drawslope )
    
    
    
def subpl_wl_dz_vel(date, adptype, resampled=False, hourgrid=False, doy=False, img=False, cbrange=[-1,1],
                    wl=True, plot_interval=1./6):
    '''
    Creates a 3 layer subplot with:
        1) water levels
        2) Dz/Dt'
        3) Velocities
    '''
    matplotlib.mathtext.SHRINK_FACTOR = 0.9
    # 1 Get velocities
    # angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}   # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
                                                          # + direction is TO outer harbour and to Lake Ontario
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    maxdepth = {'OH':6.5, 'EmbC':5.0,'Cell3':4.5}
    
    num_segments = 10
    factor = 1. / plot_interval  #  10 minutes  in days if plot_interval=1/6
    dt = 1. / 24. / factor

    hyd = get_Velocities(adptype, date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)

    rvp = []
    if img:
        # imgs do not resample
        if adptype == "OH":
            vp = vp[:-3].tolist()
        elif adptype == "EmbC":
            vp = vp[:-2].tolist()
        elif adptype == "Cell3": 
            vp = vp[:-1].tolist()
        else:
            vp = vp.tolist()
                
        if resampled:
            for i in range(0,len(vp)):
                rvtime,rvpi = hyd.resample(hyd.time, vp, dt, i)
                rvp.append(rvpi[:-10])
            rvtime = rvtime[:-10]
        else:
            rvp=vp
    else:
        if resampled:
            rvtime, rvp0 = hyd.resample(hyd.time, vp, dt, bins[adptype][0])
            rvtime, rvp2 = hyd.resample(hyd.time, vp, dt, bins[adptype][1])
            rvtime, rvp3 = hyd.resample(hyd.time, vp, dt, bins[adptype][2])
            rvp.append(rvp0[:-10])
            rvp.append(rvp2[:-10])
            rvp.append(rvp3[:-10])
        else:
            rvp.append(vp[bins[adptype][0]])
            rvp.append(vp[bins[adptype][1]])
            rvp.append(vp[bins[adptype][2]])
    
    # 2 get Dz
    daterange = [date[0], date[-1]]
    wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, daterange, dt)
    # add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
    
    legend = []
    if wl:
        legend.append("Water level") 
    legend.append("$\Delta Z/ \Delta t$")
    bin_legend = []
    for i in range(0,len(bins[adptype])):
        text = ("bin %d") % bins[adptype][i] 
        bin_legend.append(text)
    legend.append(bin_legend)  
    

    # let's do some calculations since the dz/dt will be siplayed in m/h and the information is in m/day
    # and dt0 between 2 consecutive read values in 3 minutes in thi case but shoudl be made generic since
    # we have in dt in days. Therefore we need to divide dz/dt by 24 h

    dataarr = []
    if resampled:
        if wl:
            dataarr.append(rdepths[:-10])
        dataarr.append(numpy.array(rdzdt[:-10]) / 24)
        dataarr.append(rvp)
        date = rtime[:-10]
        dateImg = rvtime
    else:
        if wl:
            dataarr.append(depths[:-1])

        dataarr.append(numpy.array(dzdt[:-1]) / 24)
        dataarr.append(rvp)
        date = dates[:-1]
        if type(hyd.time[0]) is list:
            dateImg = hyd.time[0]
        else:
            dateImg = hyd.time[0].tolist()
    # labels = ['Z [$m$]', 'dz/dt [$m$]','velocity [$m s^{-1}$]']

    if wl:
        labels = [r'$\mathsf{Z\ [m]}$', r'$\mathsf{dZ\ dt^{-1}\ [m\ h^{-1}}$]','Depth [m]']
    else:
        labels = [r'$\mathsf{dZ\ dt^{-1} [m\ h^{-1}}$]',"Depth [m]"]
    
    hyd.display_subplots(date, dateImg, dataarr, dnames=labels, yday=doy, tick=None, legend=legend,
                         hourgrid=hourgrid, img=img, cbrange=cbrange, maxdepth=maxdepth[adptype])

def calc_velocities(adptype, date, avg=True, interv = 1, resampled =False, delft3d=False):
    '''
    Calculates the Depth averaged velocity at the vertical of the ADCP
    @param interv: = 1 #hours
    @param avg: False/True. False: calculate the velocity per vertical slice every interval; 
                            True : calculate one average velocity per whole section every interval
    '''
    
    if adptype not in ['OH', 'EmbC','Cell3']:
        return
    # 1 Get velocities
     #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  and vp the along dir
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    CrossArea = {'OH':7500, 'EmbC':300,'Cell3':186}
    
    num_segments = 10
    factor  = 1./interv #factor = 1. id inter =  1 hour ; factor = 6.  # 10 minutes  in days
    
    dt = 1. / 24. / factor
    hyd = get_Velocities(adptype, date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)



    rvp = []
    for i in range(0, hyd.goodbins):
        if resampled:
            rvtime, rvp0 = hyd.resample(hyd.time, vp, dt, i)
            rvp.append(rvp0[:-10])
        else:
            rvp.append(vp[i])
    
    tottime=0
     
    if avg:
        pos_vel= []
        neg_vel= []
        
        if resampled:
            time = rvtime
            vel = numpy.array(rvp)
            dt0 = (rvtime[2] - rvtime[1])
            lenght = len(time[:-10])
        else:
            time = hyd.time
            vel= numpy.array(vp)
            dt0 = (time[0][2] - time[0][1])
            lenght = len(time[0])
        
        vel_T = numpy.transpose(vel)
        for j in range(0, lenght):
            #Mean Depth Averaged Velocity
            mvel = numpy.mean(vel_T[j])
            tottime +=dt0
            if mvel > 0: 
                pos_vel.append(mvel)
            else: 
                neg_vel.append(mvel)    
            
    volplus =  CrossArea[adptype]*numpy.sum(numpy.array(pos_vel)*dt0)*84600   
    volminus =  CrossArea[adptype]*numpy.sum(numpy.array(neg_vel)*dt0)*84600

    print("Loc=%s Vol+=%f m^3  Vol-=%f m^3" % (adptype, volplus, volminus))
    qplus=volplus/tottime 
    qminus=volminus/tottime
    print("Loc=%s Q+=%f m^3/day Q-=%f m^3/day" % (adptype, qplus, qminus))
    fname = paths[5] + "/"+ adptype+"volumes.csv"
    data = [adptype, date, qplus, qminus]
    hyd.writeVeltoCSV(fname, data, append = True, delft3d=delft3d)
    return qplus, qminus


def get_delft3d_velocities(date, adptype = 'Cell3', interv=3):

    #angles_from_N = {'OH':37, 'EmbC':127,'Cell3':137}    # for clockwise
                                                          # + sign is into Cell3   and to Cherry beach
    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise   and vp the along dir
                                                          #+ direction is TO outer harbour and to Lake Ontario    
    num_segments = 10
        
    factor = 1.  # 1 hour in days
    factor = 6.  # 10 minutes  in days
    d3d_int_min = interv
    dt = 1. / 24. / factor
    labels = [' velocity [m/s]', 'dz/dt [m/h]']

    hyd = get_Velocities(adptype, date, num_segments, delft3d=True, d3d_int=d3d_int_min)
    
   
def calc_flush_time(adptype, date, level, avg=True, interv = 1, resampled =False, umaxcorr= True):
    filemapping = {"Cell3": 'Cell3_Elev_level_maxdepth_volume_area.csv', 
                   "EmbC":'EmbC_Elev_level_maxdepth_volume_area.csv', 
                   "OH":'OuterHarbourElev_level_maxdepth_volume_area.csv' }
                   #"OH":'InnerOuterHarbourElev_level_maxdepth_volume_area.csv' }    
    qplus, qminus = calc_velocities(adptype, date, avg=avg, interv = interv, resampled=resampled)
    if umaxcorr:
        # From On the Distribution of Velocity in a V-shaped channel M. A Mohammadi, Civil Engineering Vol 16. No 1 pp 78-86
        qplus_corr = qplus/1.1934
        qminus_corr = qminus/1.1934
    else:
        qplus_corr = qplus
        qminus_corr = qminus
        
    fname = paths[5] + "/"+ filemapping[adptype]
    maxdepth, volume, area = Hydrodynamic.Hydrodynamic.ReadVolfromCSV(fname, level)
    FDT = volume/((qplus_corr-qminus_corr)/2)
    print("%s  Flushing time scale = %f [days]" % (adptype, FDT))
    fname = paths[5] + "/"+ adptype+"_FDT.csv"
    Hydrodynamic.Hydrodynamic.writeVeltoCSV(fname, [adptype, date, "Flushing TS [days]", FDT], append = True)


def calc_velocities_by_dh(adptype, date, avg=True, interv = 1, resampled =False, filemapping = None, trebitz=False):
    maxdepth = {'OH':6.5, 'EmbC':5.0,'Cell3':4.5}
    pairs = {'Cell1':["Cell1", "Cell2"], 
             'Cell2':["Cell2", "Cell3"],
             'Cell3':["Cell3", "EmbC"],
             'EmbC':["EmbC", "OH"],
             'OH':["OH", "LO"],
             'EmbA':["EmbA", "OH"],
             'EmbB':["EmbB", "OH"]}
    
    num_segments = 10
    
    factor = 6.  # 10 minutes  in days
    #factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    # 2 get Dz
    daterange = [date[0], date[-1]]
    wla1, rtime1, rdepths1, rdzdt1, dates1, depths1, dzdt1 = get_Dz_Dt(pairs[adptype][0], paths[1], num_segments, daterange, dt)
    #wla2, rtime2, rdepths2, rdzdt2, dates2, depths2, dzdt2 = get_Dz_Dt(pairs[adptype][1], paths[1], num_segments, daterange, dt)
    
    fname = paths[5] + "/"+ filemapping[adptype]
    maxdepth, volume, area = Hydrodynamic.Hydrodynamic.ReadVolfromCSV(fname, level)
    
    volplus=0
    volminus=0
    tottime=0
    Z=0  #trebitz initial sum
    for i in range(0, len(dates1)-1):
        dz = depths1[i+1]-depths1[i]
        tottime +=dates1[1]-dates1[0]  
        if trebitz:
            Z+=0.5*abs(dz)
            V=area*Z
            qplus=V/tottime 
            qminus=V/tottime
        else:
            V=area*dz
              
            if V> 0:
                volplus+=V
            else:
                volminus+=V
        
            qplus=volplus/tottime 
            qminus=volminus/tottime
    
    print("Loc=%s Q+=%f m^3/day Q-=%f m^3/day" % (adptype, qplus, qminus))
    return qplus, qminus

    
    #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.

def calc_flush_time_by_dh(adptype, date, level, avg=True, interv = 1, resampled =False, umaxcorr= True, trebitz=False):
    filemapping = {"Cell1": 'Cell1_Elev_level_maxdepth_volume_area.csv',
                   "Cell2": 'Cell2_Elev_level_maxdepth_volume_area.csv',
                   "Cell3": 'Cell3_Elev_level_maxdepth_volume_area.csv', 
                   "EmbC" : 'EmbC_Elev_level_maxdepth_volume_area.csv', 
                   "OH"   : 'OuterHarbourElev_level_maxdepth_volume_area.csv',
                   "EmbA" : 'EmbA_Elev_level_maxdepth_volume_area.csv', 
                   "EmbB" : 'EmbB_Elev_level_maxdepth_volume_area.csv', }
                   #"OH":'InnerOuterHarbourElev_level_maxdepth_volume_area.csv' }    
    qplus, qminus = calc_velocities_by_dh(adptype, date, avg=avg, interv = interv, resampled=resampled, \
                                          filemapping=filemapping,trebitz=trebitz)
        
    fname = paths[5] + "/"+ filemapping[adptype]
    maxdepth, volume, area = Hydrodynamic.Hydrodynamic.ReadVolfromCSV(fname, level)
    print("maxdepth = %f [m] volume = %f [m^3] area = %f [m^2]" % (maxdepth, volume, area))
    if trebitz:
        FDT = volume/qplus 
    else:
        FDT = volume/((qplus-qminus)/2)
    print("%s  Flusing time scale = %f [days]" % (adptype, FDT))
    fname = paths[5] + "/"+ adptype+"_DZ_FDT.csv"
    Hydrodynamic.Hydrodynamic.writeVeltoCSV(fname, [adptype, date, "Flushing TS [days]", FDT], append = True)

def convert_wl_to_delft3d_tim(path,fn,step_min, date,start_WL, timesince2001):
    num_segments = 4
    wla = Water_Level.WaterLevelAnalysis(path, [fn], num_segments, date)
    wla.convert_wl_to_delft3d_tim(path, fn, step_min,start_WL, timesince2001)

def calculate_min_since_20010101000000(path,fn, step_min, date):
    wla = Water_Level.WaterLevelAnalysis(path, [fn], 1, date)
    return wla.calculate_min_since_20010101000000(step_min)

def plot_ctd():
    path="/home/bogdan/Documents/UofT/PhD/Data_Files/2014/CTD/Individual"
    files = ["CTD01.csv","CTD02.csv","CTD03.csv","CTD04.csv","CTD05.csv","CTD06.csv","CTD07.csv","CTD15.csv",]
    Temperature.Temperature.plot_ctds(path, files)


    

def read_lake_and_harbour_data(str_date, date, water_path, harbour_path, just_lake=False):
    locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')

    print("Start read_lake_and_harboud_data")
    start_num = date[0]
    end_num = date[1]
     # 1) read all lake data
    base, dirs, files = next(iter(os.walk(water_path)))
    sorted_files = sorted(files, key = lambda x: x.split('.')[0])

    dateTimeArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
    tempArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
    resultsArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
    k = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
    i = 0
    for fname in sorted_files:
        dateTime, temp, results = utools.readTempHoboFiles.get_data_from_file(fname, window_hour, windows[1], timeinterv = date, rpath = water_path)
        maxidx = 300000
        dateTimeArr[i] = numpy.append(dateTimeArr[i], dateTime[:maxidx])
        resultsArr[i] = numpy.append(resultsArr[i], results[:maxidx])
        tempArr[i] = numpy.append(tempArr[i], temp[:maxidx])
        k[i] = numpy.append(k[i], i)
        i += 1
    # end for
    if not just_lake:
        # 1') read all harbour data (EG + Jarvis Dock
        base, dirs, files = next(iter(os.walk(harbour_path)))
        sorted_files = sorted(files, key = lambda x: x.split('.')[0])
    
        TH_dateTimeArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
        TH_tempArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
        TH_resultsArr = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
        TH_k = numpy.zeros(len(sorted_files), dtype = numpy.ndarray)
        i = 0
        for fname in sorted_files:
            dateTime, temp, results = utools.readTempHoboFiles.get_data_from_file(fname, window_hour, windows[1],
                                                                                  timeinterv=date, rpath=harbour_path)
            maxidx = 300000
            TH_dateTimeArr[i] = numpy.append(TH_dateTimeArr[i], dateTime[:maxidx])
            TH_resultsArr[i] = numpy.append(TH_resultsArr[i], results[:maxidx])
            TH_tempArr[i] = numpy.append(TH_tempArr[i], temp[:maxidx])
            TH_k[i] = numpy.append(TH_k[i], i)
            i += 1
        # end for
        return [dateTimeArr, resultsArr, tempArr, TH_dateTimeArr, TH_resultsArr, TH_tempArr]
    return [dateTimeArr, resultsArr, tempArr, [], [], []]


###############################################################################################################
def convert_tempchain_data_for_lakeStability(harbour_path=None, writepath=None, resmpl=1./24, timeinterv=None):
    # 1') read all harbour data (EG + Jarvis Dock
    base, dirs, files = next(iter(os.walk(harbour_path)))
    sorted_files = sorted(files, key=lambda x: x.split('.')[0])

    TH_dateTimeArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
    TH_tempArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
    TH_resultsArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
    TH_k = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
    i = 0
    for fname in sorted_files:
        dateTime, temp, results = utools.readTempHoboFiles.get_data_from_file(fname, window_hour, windows[1],
                                                                              timeinterv=timeinterv, rpath=harbour_path)
        maxidx = 300000
        TH_dateTimeArr[i] = dateTime[:] #numpy.append(TH_dateTimeArr[i], dateTime[:maxidx])
        TH_resultsArr[i] = results[:]  #numpy.append(TH_resultsArr[i], results[:maxidx])
        TH_tempArr[i] =  temp[:] #numpy.append(TH_tempArr[i], temp[:maxidx])
        TH_k[i] = i #numpy.append(TH_k[i], i)
        i += 1

    # resample if necessary
    if resmpl is not None:
        rDateTimeArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
        rTempArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
        rResultsArr = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)
        rk = numpy.zeros(len(sorted_files), dtype=numpy.ndarray)

        for k in range(0, len(sorted_files)):
            rdate, rdata = resample.resample(TH_dateTimeArr[k], TH_tempArr[k], resmpl)
            rDateTimeArr[k] = rdate[:]
            rResultsArr[k] = rdata[:]
            rTempArr[k] = rdata[:]
            rk[k] = k
    else:
        rDateTimeArr = TH_dateTimeArr
        rTempArr = TH_tempArr
        rResultsArr = TH_resultsArr
        rk = TH_k

    writepath = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/LakeStability'
    ofile = open(writepath + '/water_april-nov-2013.csv', "w")
    writer = csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

    # write the header first
    header = ['Date']
    for i in range(0,len(TH_dateTimeArr)):
        st="%d " % i
        header.append(st)

    writer.writerow(header)

    waterTempMatrix = numpy.zeros((len(rDateTimeArr) + 1,
                                   len(rDateTimeArr[0])), dtype=numpy.float )
    # write first row with dates
    waterTempMatrix[0] = rDateTimeArr[0]

    # write tenperature one depth per row
    for j in  range(0, len(rTempArr)):
        waterTempMatrix[j+1] = rTempArr[j]



    # transpose
    waterTempMatrix_T = numpy.transpose(waterTempMatrix)

    for i in range(0, len(waterTempMatrix_T)):
        writer.writerow(waterTempMatrix_T[i])

    ofile.close()

    return [rDateTimeArr, rResultsArr, rTempArr]


def subplot_lake_harbour(str_date, date, adptype, just_lake=False):
    resampled =False
    locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')

    water_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-LakeOntario/csv_processed'
    water_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2012/MOE-Apr-May_2012-Thermistor_chain/csv_processed'
    harbour_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/AllHarbour/csv_processed/EGap-JarvisDock'
    water_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2015/loboviz2-temp/csv_processed'
    print("Start subplot_lake_harbour ")
    start_num = date[0]
    end_num = date[1]

    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  and using vp along 
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    maxdepth = {'OH':6.5, 'EmbC':5.0,'Cell3':4.5}
    
    num_segments = 10
    
    factor = 6.  # 10 minutes  in days
    factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    if not just_lake:
        hyd = get_Velocities(adptype, str_date, num_segments)
    
        #reproject on the direction of thde flow
        Theta = angles_from_N[adptype] #35.1287  # degrees
        tet = 2 * math.pi * Theta / 360
        up, vp = hyd.rotation_transform(tet, clockwise=False)
    
        rvp = []
        #imgs do not resample
        if adptype == "OH":
            vp = vp[:-2].tolist()
            hyd.time = hyd.time[:-2]
        elif adptype == "EmbC":
            vp = vp[:-1].tolist()
            hyd.time = hyd.time[:-1]
        elif adptype == "Cell3": 
            #up = up[:-1].tolist()
            vp = vp.tolist()
        else:
            vp = vp.tolist()
             
        #reverse the list for plotting   
        vp = vp[::-1]
                
        if resampled:
            for i in range(0,len(up)):
                rvtime,rvpi = hyd.resample(hyd.time, vp, dt, i)
                rvp.append(rvpi[:-10])
            rvtime = rvtime[:-10]
        else:
            rvp=vp
            
        # 2 get Dz
        daterange = [str_date[0], str_date[-1]]
        wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, daterange, dt)
        #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
        
        legend = []
        legend.append("Water level") 
    #end if
    
    # 1) read all lake data
    # 1') read all harbour data (EG + Jarvis Dock

    [dateTimeArr, resultsArr, tempArr, TH_dateTimeArr, TH_resultsArr, TH_tempArr] = \
        read_lake_and_harbour_data(str_date, date, water_path, harbour_path, just_lake)

    if not just_lake: 
        # 2'') read temperature diff at TC4
        diffpath = water_path + "/../../TC-OuterHarbour/colormap"
        base, dirs, files = next(iter(os.walk(diffpath)))
    
        sorted_files = sorted(files, key = lambda x: x.split('.')[0])
    
        cmdateTimeArr = []
        cmtempArr = []
        cmresultsArr = []
        cmk = []
        j = i = 0
        nsorted_files = sorted_files[:]
        for fname in sorted_files:
            print("Filename = %s" % fname)
            dateTime, temp, results = utools.readTempHoboFiles.get_data_from_file(fname, window_6hour, windows[1], timeinterv = date, rpath = diffpath)
            cmdateTimeArr.append(dateTime)
            cmresultsArr.append(results)
            cmtempArr.append(temp)
            cmk.append(j)
            i += 1
            j += 1
        # end for
    else:    
        # plot the temperature not the smoothed ones - ONLY for test
        datetype = 'dayofyear'
        # display_data.display_temperatures(dateTimeArr, tempArr, k, fnames = sorted_files, custom = "Temperature Toronto Waterfront Zones $^oC$", \
        #                                  datetype = datetype)
            
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
            
        maxdepth = 27
        firstlogdepth = 3
        maxtemp = 26
        mintemps = 0
        mindepths = 0
        t31 = ['0', '6', '13', '20', '27']
        t32 = [27, 20, 13.5, 6.5, 0]
        tick = [t31, t32]

        # trim all the arrays to the min size
        ms = len(dateTimeArr[0])
        for i in range(0, len(dateTimeArr)):
            ms = min(ms, len(dateTimeArr[i]))

        for i in range(0, len(dateTimeArr)):
            dateTimeArr[i] = dateTimeArr[i][:ms]
            resultsArr[i] = resultsArr[i][:ms]
            tempArr[i] = tempArr[i][:ms]
        
        utools.display_data.display_img_temperatures(dateTimeArr, tempArr, resultsArr, 1, tick, maxdepth, firstlogdepth, maxtemp, \
                                                     fontsize = 18, datetype = datetype, thermocline = True, interp = None, \
                                                     ycustom = None, cblabel = "T [$^\circ$C]", draw_hline = False, hline_freq = 2)


    
    if not just_lake: 
        # Time domain analysis
        lowcut = 1.0 / (24 * 10) / 3600
        highcut = 1.0 / (24 * 3) / 3600
        tunits = 'day'
    
        if tunits == 'day':
            factor = 86400
        elif tunits == 'hour':
            factor = 3600
        else:
            factor = 1
    
        
        # 3) Mixed water, air ,img data
        custom = numpy.array([r'Wind spd [$\mathsf{km\ h^{-1}}$]', 'Air [T($^\circ$C]', 'Water T[$^\circ$C]', ])
        # ToDO: Add short and long radiation
        print("Start display mixed subplots ")
        
        dateTimes3 = [hyd.time, TH_dateTimeArr, dateTimeArr]
        ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]"] 
        #imgs = [TH_resultsArr, resultsArr] #Harbour first
        imgs = [vp,TH_resultsArr, resultsArr]
        
        
        t11 = ['0', '2', '4.5', '7']
        t12 = [7, 4.5, 2, 0]
        
        t21 = ['0', '3', '6', '10']
        t22 = [10, 6, 3, 0]
        
        t31 = ['0', '6', '13', '20', '27']
        t32 = [27, 20, 13.5, 6.5, 0]
        
        
        tick = [[t11, t12], [t21, t22], [t31, t32]]
        
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
            
        #maxdepth = [9, 27] # Harbour first
        maxdepth = [7, 9, 27]
        #firstlogdepth = [0, 3] Harbour first
        firstlogdepth = [0, 0, 3]
        
        maxtemp = [0.5, 26, 26]
        mintemps = [-0.5, 0, 0]
        mindepths = [0, 0, 0]
        
        utools.display_data.display_mixed_subplot(dateTimes1 = [rtime[:-1]], data = [rdepths[:-1]], varnames = ["Water level"], ylabels1 = ["Depth [m]"], \
                                           dateTimes2 = [], groups = [], groupnames = [], ylabels2 = [], \
                                           dateTimes3 = dateTimes3, imgs = imgs, ylabels3 = ylabels3, ticks = tick, maxdepths = maxdepth, \
                                            mindepths = mindepths, mintemps = mintemps, firstlogs = firstlogdepth, maxtemps = maxtemp, \
                              fnames = None, revert = False, custom = None, maxdepth = None, tick = None, firstlog = None, yday = True, \
                              title = False, grid = False, limits = limits, sharex = True, fontsize = 18, group_first = False, interp = 2)
    #end if



def subplot_Dz_ADCP_T_harbour(str_date, date, adptype, one_V = False):
    '''

    :param str_date:
    :param date:
    :param adptype:
    :param one_V:
    :return:
    '''

    resampled =False
    locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')
    matplotlib.mathtext.SHRINK_FACTOR = 0.9

    water_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/TC-LakeOntario/csv_processed'
    harbour_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/AllHarbour/csv_processed/EGap-JarvisDock'

    print("Start sublpot_Dz_ADCP_harbour ")
    start_num = date[0]
    end_num = date[1]

    

    angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  and using vp along 
                                                          #+ direction is TO outer harbour and to Lake Ontario 
    bins = {'OH':[0,2,4], 'EmbC':[0,2,3],'Cell3':[0,2,3]}
    maxdepth = {'OH':6.5, 'EmbC':5.0,'Cell3':4.5}
    
    num_segments = 10
    
    factor = 6.  # 10 minutes  in days
    factor = 1.  # 1 hour in days
    dt = 1. / 24. / factor

    hyd = get_Velocities(adptype, str_date, num_segments)

    #reproject on the direction of thde flow
    Theta = angles_from_N[adptype] #35.1287  # degrees
    tet = 2 * math.pi * Theta / 360
    up, vp = hyd.rotation_transform(tet, clockwise=False)

    rvp = []
    #imgs do not resample
    if adptype == "OH":
        vp = vp[:-2].tolist()
        hyd.time = hyd.time[:-2]
    elif adptype == "EmbC":
        vp = vp[:-1].tolist()
        hyd.time = hyd.time[:-1]
    elif adptype == "Cell3": 
        #up = up[:-1].tolist()
        vp = vp.tolist()
    else:
        vp = vp.tolist()
         
    #Do not reverse the list for plotting   
    #vp = vp[::-1]
            
    if resampled:
        for i in range(0,len(up)):
            rvtime,rvpi = hyd.resample(hyd.time, vp, dt, i)
            rvp.append(rvpi[:-10])
        rvtime = rvtime[:-10]
    else:
        rvp=vp

    # 3) ADCP data 
    # ToDO: Add short and long radiation
    print("Grad Vel calc ")
    grad = []
    dz =1
    V_T  = numpy.transpose(vp) 
    U_T  = numpy.transpose(up)
           
    for i in range(0, len(hyd.time[0])):
        dvdz = []
        for j in range(0, len(hyd.time)-1):
            W1 = math.sqrt(V_T[i][j+1]** 2 + U_T[i][j+1]** 2)
            W2 = math.sqrt(V_T[i][j] ** 2 + U_T[i][j]** 2)
            dvdz.append((W1-W2)/dz/2)
        #edn for 
        dvdz.insert(0, (0+dvdz[0])/2.)    #insert fist item to match the velocity array size
        grad.append(dvdz)
        #grad.append(numpy.gradient(, dz))
    #end for
    
    gradarr = numpy.array(grad)
    grad_T = numpy.transpose(gradarr)

    if adptype == "Cell3":
        dd = 1 # for Cell 3
    elif adptype == "EmbC":
        dd= 1  #2
    elif adptype == "OH":
        dd = 1 #3
        
    # 2 get Dz
    daterange = [str_date[0], str_date[-1]]
    wla, rtime, rdepths, rdzdt, dates, depths, dzdt = get_Dz_Dt(adptype, paths[1], num_segments, daterange, dt)
    #add one more to dzdt since is intepolated between nodes and prfiles are on nodes.
    limits1 = [[-1.2, 1.2]]
    legend = []
    legend.append("Water level") 
    
    # 1) read all lake data
    # 1') read all harbour data (EG + Jarvis Dock)

    [dateTimeArr, resultsArr, tempArr, TH_dateTimeArr, TH_resultsArr, TH_tempArr] = \
        read_lake_and_harbour_data(str_date, date, water_path, harbour_path)

    # 2'') read temperature diff at TC4
    diffpath = water_path + "/../../TC-OuterHarbour/colormap"
    base, dirs, files = next(iter(os.walk(diffpath)))

    sorted_files = sorted(files, key = lambda x: x.split('.')[0])

    cmdateTimeArr = []
    cmtempArr = []
    cmresultsArr = []
    cmk = []
    j = i = 0
    nsorted_files = sorted_files[:]
    for fname in sorted_files:
        print("Filename = %s" % fname)
        dateTime, temp, results = utools.readTempHoboFiles.get_data_from_file(fname, window_6hour, windows[1], timeinterv = date, rpath = diffpath)
        cmdateTimeArr.append(dateTime)
        cmresultsArr.append(results)
        cmtempArr.append(temp)
        cmk.append(j)
        i += 1
        j += 1
    # end for

    # plot the temperature not the smoothed ones
    ndateTimeArr = numpy.array(cmdateTimeArr)
    ntempArr = numpy.array(cmtempArr)
    nresultsArr = numpy.array(cmresultsArr)
    nk = numpy.array(cmk)
    diffarr = numpy.array(ntempArr[1]) - numpy.array(ntempArr[0])

    # plot the temperature not the smoothed ones - ONLY for test
    # datetype = 'dayofyear'
    # display_data.display_temperatures(dateTimeArr, tempArr, k, fnames = sorted_files, custom = "Temperature Toronto Waterfront Zones $^oC$", \
    #                                  datetype = datetype)
    # display_data.display_img_temperatures(dateTimeArr, tempArr, resultsArr, k, tick, maxdepth, firstlogdepth, maxtemp, fontsize = 18, datetype = datetype)

    # Time domain analysis
    lowcut = 1.0 / (24 * 10) / 3600
    highcut = 1.0 / (24 * 3) / 3600
    tunits = 'day'

    if tunits == 'day':
        factor = 86400
    elif tunits == 'hour':
        factor = 3600
    else:
        factor = 1

    
    # 3) Mixed water, air ,img data
    custom = numpy.array(['Water level', ])
    
        
    
    if one_V:
        
        if adptype == "Cell3":
            dateTimes3 = [hyd.time, TH_dateTimeArr]
            imgs = [vp, TH_resultsArr[::-1]]
        else:
            dateTimes3 = [hyd.time[:-dd], TH_dateTimeArr]
            imgs = [vp[:-dd], TH_resultsArr[::-1]]
            
        ylabels3 = ["Depth [m]", "Depth [m]"] 
        
        if adptype == 'OH':
            t11 = ['0', '2', '4.5', '7']
            t12 = [7, 4.5, 2, 0]

            #maxdepth = [9, 27] # Harbour first
            maxdepth = [7, 10]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0]
            mindepths = [0, 0]
            maxtemp = [0.2,  24]
            mintemps = [-0.2, 0]
            
        elif adptype == "EmbC" or adptype == "Cell3":
            t11 = ['0', '1.5', '3.0', '4.5']
            t12 = [4.5, 3, 1.5, 0]
           
            #maxdepth = [9, 27] # Harbour first
            maxdepth = [4.5, 10]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0]
            mindepths = [0, 0]
            maxtemp = [0.2, 24]
            mintemps = [-0.2, 0]
            
       
        t41 = ['0', '3', '6', '10']
        t42 = [10, 6, 3, 0]
        
        tick = [[t11, t12], [t41, t42]  ]
        
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
        clabel =  ["$V_{alg}$ [$m s^{-1}$]", "Temp [$\circ$ C]"]   
    else:
        if adptype == "Cell3":
            dateTimes3 = [hyd.time, hyd.time, hyd.time, TH_dateTimeArr]
            imgs = [vp, up, grad_T, TH_resultsArr[::-1]]
        else:
            dateTimes3 = [hyd.time[:-dd], hyd.time[:-dd], hyd.time[:-dd], TH_dateTimeArr]
            imgs = [vp[:-dd], up[:-dd], grad_T[:-dd], TH_resultsArr[::-1]]
            
        ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]"] 
        if adptype == 'OH':
            t11 = ['0', '2', '4.5', '7']
            t12 = [7, 4.5, 2, 0]
            
            t21 = ['0', '2', '4.5', '7']
            t22 = [7, 4.5, 2, 0]
            
            t31 = ['0', '2', '4.5', '7']
            t32 = [7, 4.5, 2, 0]
        
            #maxdepth = [9, 27] # Harbour first
            maxdepth = [7, 7, 7, 10]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0, 0, 0]
            mindepths = [0, 0, 0, 0]
            maxtemp = [0.2, 0.2, 0.05, 24]
            mintemps = [-0.2, -0.2, -0.05, 0]
            
        elif adptype == "EmbC" or adptype == "Cell3":
            t11 = ['0', '1.5', '3.0', '4.5']
            t12 = [4.5, 3, 1.5, 0]
            
            t21 = ['0', '1.5', '3.0', '4.5']
            t22 = [4.5, 3, 1.5, 0]
            
            t31 = ['0', '1.5', '3.0', '4.5']
            t32 = [4.5, 3, 1.5, 0]
        
            #maxdepth = [9, 27] # Harbour first
            maxdepth = [4.5, 4.5, 4.5, 10]
            #firstlogdepth = [0, 3] Harbour first
            firstlogdepth = [0, 0, 0, 0]
            mindepths = [0, 0, 0, 0, 0]
            maxtemp = [0.2, 0.2, 0.05, 24]
            mintemps = [-0.2, -0.2, -0.05, 0]
            
       
        t41 = ['0', '3', '6', '10']
        t42 = [10, 6, 3, 0]
        
        tick = [[t11, t12], [t21, t22], [t31, t32], [t41, t42]  ]
        
        #limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers 
        limits = None
        clabel = [r'$\mathsf{U_{alg}}$ [$\mathsf{m\ s^{-1}}$]',
                   r'$\mathsf{V_{crs}}$ [$\mathsf{m\ s^{-1}}$]',
                   r'$\mathsf{dV\ dz^{-1}}$ [$\mathsf{s^{-1}}$]',
                   r'$\mathsf{T_{water}}$ [$^\circ$C]']
    
    
    utools.display_data.display_mixed_subplot(dateTimes1 = [rtime[:-1]], data = [rdzdt[:-1]], varnames = [],
                                              ylabels1 = [r"dZ dt$^{-1}$ [$\mathsf{m\ h^{-1}}$]"], limits1 = limits1,
                                              dateTimes2 = [], groups = [], groupnames = [], ylabels2 = [],
                                              dateTimes3 = dateTimes3, imgs = imgs, ylabels3 = ylabels3,
                                              ticks = tick, maxdepths = maxdepth,
                                              mindepths=mindepths, mintemps=mintemps, firstlogs=firstlogdepth,
                                              maxtemps = maxtemp,
                                              fnames = None, revert = True, custom = None, maxdepth = None,
                                              tick=None, firstlog=None, yday=True,
                                              title = False, grid=False, limits=limits, sharex=True,
                                              fontsize=18, group_first=False, interp=2, cblabel=clabel)
    
if __name__ == '__main__':

    # interval where velocity data is available ( temerature data is available from 1 June-- 30 October)
    date = ['13/06/25 00:00:00', '13/08/16 00:00:00']  # this interval is needed to cover completely the intersection between velocity data and water level data
    # zoom in
    # date = ['13/07/18 00:00:00', '13/07/27 00:00:00']
    # Upwelling - mostly
    #date = ['13/07/09 00:00:00', '13/07/29 00:00:00']
     # No upwelling - mostly
    #date = ['13/07/29 00:00:00', '13/08/12 00:00:00']
    
    #upw  - one day
    #date = ['13/07/18 13:27:00', '13/07/18 13:50:00']
    #relaxation one day
    #date = ['13/07/21 00:00:00', '13/07/24 00:00:00']

    #v = 'tobermory'  # for doing the XCT analysis on FFNMP data
    #v = 'hodographs'
    #v = 'windrose_vel'
    #v = 'dz_dt'
    v = 'subpl_wl_dz_vel'
    v = 'vel-profiles'
    #v = 'wl_fft_all'
    #v = 'vel_fft_pairs'
    #v = 'temp_fft'
    #v = 'wl_fft_pairs'
    #v = 'plot_fft_v_T_wl'
    #v = 'calc_vel_flush'
    #v = 'calc_vel_flush_dh'
    #v = 'temp_fft_all'
    #v = 'conv_wl_delft3d_min'
    #v = 'ctd'
    v = "subplot_lake_harbour"
    #v = 'subplot_Dz_ADCP_T_harbour'
    v = 'subplot_Temp_OH_2015'
    #v = "wct_v_t"
    #v = 'avg-vel-profiles'
    #v = "dT_meas_calc"
    #v = 'calc_delft3d_vel'
    # v = "stddev_T"
    # v = "convert_data_for_lakeStability"

    # map the inputs to the function blocks
    for case in switch(v):
        if case('wl_fft_all'):
            WL_FFT_analysis_all()
            break
        if case('temp_fft_all'):
            Temp_FFT_analysis_all()
            break
        if case('wl_fft_pairs'):
            dates = ['13/08/03 00:00:00', '13/08/04 00:00:00']
            dates = ['13/06/06 00:00:00', '13/10/07 00:00:00']
            WL_FFT_pairs(dateinterval=dates)
            break
        if case('vel_fft_pairs'):
            Vel_FFT_pairs(date, skipRDI=True)
            break
        if case('temp_fft'):
            Temp_FFT(date, skipRDI=True)
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
            factor = 1.  # 1 hour in dayslayer subplot with:
            factor = 6.  # 10 minutes  in days
            dt = 1. / 24. / factor
            #modd = None
            modd = 24 # every nth value of the dataset gets an arrow
            Vel_hodographs(date, dt, modd)
            break
        if  case ('avg-vel-profiles'):
            #location can be 'Cell3' 'EmbC' 'OH'
            location = 'Cell3' 
            location = 'EmbC'
            location = 'OH'
            
            firstbins = {"Cell3":0.3, "EmbC":1.0, "OH":1.0}
            intervals= {"Cell3":1.0, "EmbC":1.0, "OH":1.0}
            Avg_Vel_profiles(date, adptype = location, firstbin = firstbins[location], \
                         interval = intervals[location], eddyvisc = True, ADCPprof= False, grad_valong = True)
            break
        
        if case ('windrose_vel'):
            skipRDI = False
            Vel_windrose(date, skipRDI)
            break
        
        if  case ('vel-profiles'):
            #location can be 'Cell3' 'EmbC' 'OH'
            #location = 'EmbC'
            location = 'OH'
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
            date = ['13/07/24 00:00:00', '13/07/24 23:00:00']    # one day shows hourly alternation in DZ
            #date = ['13/07/19 23:00:00', '13/07/21 13:00:00']    # stormy period
            resampled = True
            hourgrid = True
            DOY = False
            img=True
            location = "EmbC"
            location = "Cell3"
            wl = True
            plot_int = 1./6. # in hours -> 1/6 =10 minutes
            subpl_wl_dz_vel(date, location, resampled, hourgrid, doy = DOY, img=img, cbrange=[-0.4,0.4],
                            wl=wl, plot_interval=plot_int)
            break
        if case ('plot_fft_v_T_wl'):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            location = "OH"
            location ="Cell3"
            #location="EmbC"
            drawslope = False
            plot_FFT_V_T_WL(location, date, scale = 'log', drawslope = drawslope)
            break
        if case ('calc_delft3d_vel'):
            date = ['13/08/01 00:00:00', '13/08/31 00:00:00']    # TOO LONG
            interv = 3 #min
            loc="OH"
            loc="Cell3"
            #loc="EmbC"
            get_delft3d_velocities(date, adptype = loc, interv=interv)
            break
        
        if case ('calc_vel_flush'):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            date = ['13/07/10 00:00:00', '13/07/28 23:00:00']    # Upwelling
            date = ['13/07/29 00:00:00', '13/08/12 23:00:00']    # NO Upwelling
            #date = ['13/07/29 00:00:00', '13/08/11 00:00:00']    # quiet warm
            #date = ['13/06/19 00:00:00', '13/06/30 00:00:00']    # quiet cold
            
            interv = 1 #hours
            avg = True # False calculat the vel per veritcal slice; True = calucalte one avg vel per whole section
            loc="OH"
            loc="Cell3"
            #loc="EmbC"
            level = 75
            #calc_velocities(loc,date, avg=avg, interv=interv, resampled = False)
            calc_flush_time(loc, date, level, avg=avg, interv = interv, resampled =False)
            break
        if case ('calc_vel_flush_dh'):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            date = ['13/07/10 00:00:00', '13/07/28 23:00:00']    # Upwelling
            date = ['13/07/29 00:00:00', '13/08/12 23:00:00']    # NO Upwelling
            #date = ['13/07/29 00:00:00', '13/08/11 00:00:00']    # quiet warm
            #date = ['13/06/19 00:00:00', '13/06/30 00:00:00']    # quiet cold
            
            interv = 1 #hours
            avg = True # False calculat the vel per veritcal slice; True = calucalte one avg vel per whole section
            loc="OH"
            loc="Cell3"
            #loc="EmbC"
            #loc="Cell2"
            #loc="Cell1"
            #loc="EmbB"
            #loc="EmbA"
            level = 74
            ADCP=False
            if ADCP:
                calc_velocities(loc, date, avg=avg, interv=interv, resampled = False)
            trebitz=True
            calc_flush_time_by_dh(loc, date, level, avg=avg, interv = interv, resampled =False, trebitz=trebitz)
            break
        if case ("conv_wl_delft3d_min"):
            date = ['13/06/25 00:00:00', '13/08/16 00:00:00']    # TOO LONG
            start_WL = {'13/06/25 00:00:00':75.16}
            pathwl=paths[1]
            fn="10279444_corr.csv"
            step_min = 60
            ta = calculate_min_since_20010101000000(pathwl, fn, step_min, date)
            print(ta)
            convert_wl_to_delft3d_tim(pathwl, fn, step_min, date, start_WL[date[0]], ta)
            break
        if case ("ctd"):
            plot_ctd()
            break
        if case ("subplot_lake_harbour"):
            #for Liset 2012
            date = ['12/06/15 00:00:00', '12/09/30 00:00:00']
            date = ['15/06/01 12:30:00', '15/11/04 12:30:00']
            dt = datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
            start_num = dates.date2num(dt)
            dt = datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
            end_num = dates.date2num(dt)
            adptype = "EmbC"
            subplot_lake_harbour(date, [start_num, end_num], adptype, just_lake=True)
            break
        if case ("subplot_Dz_ADCP_T_harbour"):
            dt = datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
            start_num = dates.date2num(dt)
            dt = datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
            end_num = dates.date2num(dt)
            adptype = "Cell3"
            #adptype = "OH"
            one_V = True
            subplot_Dz_ADCP_T_harbour(date, [start_num, end_num], adptype, one_V)
            break
        if case ("subplot_Temp_OH_2015"):
            date = ['15/06/11 12:30:00', '15/11/4 12:30:00']
            dt = datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
            start_num = dates.date2num(dt)
            dt = datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
            end_num = dates.date2num(dt)
            Temperature.Temperature.subplot_Temp_OH_2015(date, [start_num, end_num])
            break
        
        if case ("wct_v_t"):
            bin = 0
            adptype = 'OH'
            wct_V_T(date, bin, adptype)
            break

        if case ("convert_data_for_lakeStability"):
            writepath = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/LakeStability'
            harbour_path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/AllHarbour/csv_processed/EGap-JarvisDock'
            date = ['13/04/29 00:00:00', '13/11/10 00:00:00']
            dt = datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
            start_num = dates.date2num(dt)
            dt = datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
            end_num = dates.date2num(dt)
            timeinterv =[start_num, end_num]
            convert_tempchain_data_for_lakeStability(harbour_path=harbour_path,
                                                     writepath=writepath,
                                                     resmpl=1./24,               # 1h in days
                                                     timeinterv=timeinterv)
      
        if case("dT_meas_calc"):
            ''' 
            Station: ea,  mean: 12.364103  max: 16.654000 min: 5.616000
            Station: ohtc3,  mean: 10.264501  max: 17.153000 min: 4.844000
            Station: ec,  mean: 12.833353  max: 17.177000 min: 7.695000
            Station: oh21,  mean: 13.641496  max: 20.627000 min: 4.895000
            Station: eb,  mean: 18.505002  max: 25.574000 min: 8.145000
            Station: c2,  mean: 20.343706  max: 26.940000 min: 14.888000
            Station: lo12,  mean: 9.855198  max: 16.773000 min: 4.350000
            Station: c3,  mean: 19.853106  max: 25.598000 min: 14.098000
            Station: oh02,  mean: 10.955305  max: 19.222000 min: 5.050000
            Station: c1,  mean: 23.986042  max: 30.268000 min: 19.936000
            Chain 0):
               bay: 1,  dt(oh-lake): 1.765236
               bay: 2,  dt(ec-lake): 2.978155
               bay: 3,  dt(c3-lake): 9.997908
               bay: 4,  dt(c2-lake): 10.488507
               bay: 5,  dt(c1-lake): 14.130844
            Chain 1):
               bay: 1,  dt(oh-lake): 1.765236
               bay: 2,  dt(eb-lake): 8.649804
            Chain 2):
               bay: 1,  dt(oh-lake): 1.765236
               bay: 2,  dt(ea-lake): 2.508904
            ----------------------------------------------------------------------------------
            For H=74
            ----------------------------------------------------------------------------------
            Station    Meas. ΔT [ºC]    Res time [days]    Mean  depth [m]    Mean temp [ºC]
            ----------------------------------------------------------------------------------
            OH              3               1.5                 9.5             12.8
            Emb. C         6.4              2.8                 3.5             16.4
            Cell 3         7.7             6/8.6              7.6/6.2           18.1
            Cell 2         10.4            1.18/5.18           1.6/1.2          20.4
            Cell 1 Shidan  13.2             9.7                 0.6             24.2
            Emb  A                          1.25                2.6
            Emb  Ahf                        1.65                2.6
            Emb  Alf                        1.95                2.6      
            Emb  B                          0.6                 1.19            
            '''
            #stations  = ["OH", "Emb C", "Cell 3", "Cell 2", "Cell 1", "Emb B", "Emb A"]
            stations  = ["Cell 1", "Cell 2", "Cell 3", "Emb A", "Emb B", "Emb C", "OH"]
            #depths    = [8.9 , 3.5, 7.2, 1.2, 1.3, ,]
            
            range_x   = 10
            range_y   = 10
            # The valuse below are credible and not fabricated
            # Calculated Residence time, mean depth, delta T meas, mean temp 
            stddev = {"C1": 2.693083-1, 
                     "C2": 2.514186-1,
                     "C3": 2.1890-1,
                     "EA": 3.467984-1,    
                     "Ah": 3.467984-1,
                     "Al": 3.467984-1,
                     "EB": 3.383678-1,
                     "EC": 2.084654-1,
                     "OH": 3.328926-1,
                     "LO": 3.034539-1}

            points = {
                      "C1": [[5.7, 1.3, 14.1, 23.98], stddev["C1"]],          # with the current area/vol - res time =2.17
                      "C2": [[2.1, 0.8, 10.4, 20.34], stddev["C2"]],          # res 1.1 ->2.1 
                      "C3": [[6.9,  6.1, 5.9, 15.3], stddev["C3"]],           # res 6 -> 6.9  dept 6.8 -> 6.1 mean deltaT 9->4.5 temp 19->14.5         
                      "EA": [[1.3, 2.6, 2.5, 12.36], stddev["EA"]],
                      "Ah": [[2.6, 2.6, 2.5, 12.36], stddev["Ah"]],
                      "Al": [[3.0, 2.6, 2.5, 12.36], stddev["Al"]],
                      "EB": [[0.6, 1.19, 8.6, 18.5], stddev["EB"]],           #
                      "EC": [[2.8, 3.3, 3.0, 12.8], stddev["EC"]],
                      "OH": [[3.2, 8.9, 1.7, 11.4], stddev["OH"]],             # mean temp an average of 3  
                      "LO": [[0, 80, 0, 9.85], stddev["LO"]],
                      }
            

            #name and mean temp
            chain1={"LO":points["LO"][0][3], "OH":points["OH"][0][3], "EC":points["EC"][0][3], "C3":points["C3"][0][3], "C2":points["C2"][0][3], "C1":points["C1"][0][3]}
            chain2={"LO":points["LO"][0][3], "OH":points["OH"][0][3], "EB":points["EB"][0][3]}
            chain3={"LO":points["LO"][0][3], "OH":points["OH"][0][3], "EA":points["EB"][0][3]}
            chains=[chain1, chain2, chain3]
            
            #Temperature.Temperature.Dt_Tres_H_depth([1, 10], 100, [1,10], 100, [10,350], 3,\
            Temperature.Temperature.DTemp_Model_vs_Meas([0, 10], 100, [0,10], 100, [150, 500], 1,
                                                     xlabel="Residence time [days]",
                                                     ylabel="Depth [m]",
                                                     cblabel="$\Delta$T [$^\circ$C]",
                                                     points=points,
                                                     lake = {"LO": [0, 80, 0, 9.85]},
                                                     range_x =range_x, range_y=range_y,
                                                     fontsize = 20,  chains = chains, exclusions =["LO", "EB"]  )
            break
        
        if ("stddev_T"):
            filenames=["Cell_1A.csv","Cell2.csv","Cell3.csv","Stn14.csv","Stn13.csv", "EmbC.csv", "Stn20.csv", "LO_10_2393006.csv"]
            
            tinterv = ["13/06/01 00:00:00","13/09/01 00:00:00"]
            temp = Temperature.Temperature(paths[7], filenames, 3, tinterv = tinterv)
            stddev = temp.stddev()
            i=0
            for f in filenames:
                print ("Fname:%s, minSD:%f, maxSD:%f" % (f, stddev[i][0],stddev[i][1]))
                i+=1 
            break
        if case():  # default, could also just omit condition or 'if True'
            
            print ("something else!")
            # No need to break here, it'll stop anyway

    print("Done!")
