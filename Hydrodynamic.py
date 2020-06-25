'''
Created on Nov 20, 2014

@author: bogdan
'''
from rdradcp.readRawBinADCP import readRawBinADCP
import rdradcp.plot_ADCP_vel_FFT as plot_ADCP_vel_FFT
import rdradcp.plot_ADCP_velocity as plot_ADCP_velocity
import rdradcp.writebins as writebins
import rdradcp.Timer
import csv
import pcadcp.rdTxtPcAdp as rdTxtPcAdp

import ufft.spectral_analysis
import ufft.FFTGraphs as fftGraphs
import utools.windows

import utools.display_data as display_data
import ufft.smooth as smooth
import utools.readTempHoboFiles as readTempHoboFiles

import datetime, math
import matplotlib.dates as dates
import numpy, scipy
import Water_Level
import utools.interpolate_data
import matplotlib

class Hydrodynamic(object):
    '''
    classdocs
    '''

    def __init__(self, filename, sitename, path, num_segments, tinterv, ctlname = None, filenames = None):
        '''
        Constructor
        '''
        self.path = path
        self.filename = filename
        self.num_segments = num_segments
        self.adcp = None
        self.ens = None
        self.cfg = None
        self.hdr = None
        self.results_u = None
        self.results_v = None
        self.results_z = None
        self.results_temp = None
        self.time = None
        self.locname = sitename
        self.goodbins = 0
        if tinterv != None:
            self.date = [tinterv[0], tinterv[1]]

        self.ctlname = ctlname
        self.filenames = filenames

        if ctlname != None:
            self.pcadp = rdTxtPcAdp.PcADP(filename, path, ctlname, filenames, num_segments, tinterv)
        else:
            self.pcadp = None

    def getData(self):
        return [self.locname, self.time, self.results_u, self.results_v]

    def getDateIndex(self, date, timearr):
        # select the data by dates
        dt = datetime.datetime.strptime(date, "%y/%m/%d %H:%M:%S")
        datenum = dates.date2num(dt)
        
        sz = len(timearr)
        idx = 0
        found = False
        for j in range(0, sz):
            d = timearr[j]
            #print "Datetime: %s" % d
            if d >= datenum:
                idx =j
                found = True
                break
        if found:    
            return idx
        else:
            raise Exception("Datetime index not found!" )

    @staticmethod
    def get_date_num(date, FORMAT):
        '''
        Convert a string date to a UNIX  date number
        '''
        #dt = datetime.datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
        dt = datetime.datetime.strptime(date, FORMAT)
        start_num = dates.date2num(dt)
        return start_num

    @staticmethod
    def plotPVD(data, dt, bins, tunit, vunit, lunit, fontsize = 20, title = False, modd = 9):
        [name, time, results_u, results_v] = data

        # resample
        ur = []
        vr = []
        for b in bins:
            ur.append(Hydrodynamic.resample(time, results_u, dt, b)[1])  # [0] is time returned by resample
            vr.append(Hydrodynamic.resample(time, results_v, dt, b)[1])

        # display
        display_data.plot_hodograph(name, ur, vr, bins, dt, tunit = tunit, vunit = vunit, \
                                    disp_lunit = lunit, fontsize = fontsize, title = title, modd = modd)


    #-------------------------------------------------------------------------------------------------------------------
    def get_Velocities(self, adptype, date, num_segments, delft3d=False, d3d_int=3, adcp=None):

        if adcp is not None:
            self.adcp = adcp
            self.cfg = adcp.config
            self.ens = None
            self.hdr = None
        elif adptype == 'EmbC' or adptype == 'OH' or adptype == 'WGap' or adptype == 'EGap' or adptype == 'TI' :
            self.readRawBinADCP()
        else:
            self.readPcAdpVel()

        self.select_data_dates(delft3d=delft3d, d3d_int=d3d_int)
            # these ones have dates already selected from reading. no select_data_dates is necessary
        # endif
        return self


    #-------------------------------------------------------------------------------------------------------------------
    def readPcAdpVel(self, velprofiles = False):
        velocities, timevec, goodbins = self.pcadp.readVel()
        self.results_u = velocities['ve']
        self.results_v = velocities['vn']
        self.results_z = velocities['vu']
        self.time = timevec
        self.goodbins = goodbins
        if velprofiles:
            #display_data.display_temperature([self.time[0]], [self.results_z[0]], [self.results_z[0]] , [[1]], fnames = ['Vert. Vel.'])
        #=======================================================================
        #     angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
        #                                                   #+ direction is TO outer harbour and to Lake Ontario    
        #     num_segments = 10
        # 
        #     Theta = angles_from_N[self.locname] #35.1287  # degrees
        #     tet = 2 * math.pi * Theta / 360
        #     up, vp = plot_ADCP_velocity.rotation_transform(self.results_u, self.results_v, tet, clockwise=False)
        #     dz = 1
        #     vp_T = numpy.transpose(vp)
        #     z_T  = numpy.transpose(self.results_z)
        #     Kv_lst = []
        #     for i in range(0, len(z_T)):
        #         Kv = Hydrodynamic.eddy_viscosity(vp_T[i], z_T[i], dz)
        #         Kv_lst.append(Kv)
        #     Kv_arr = numpy.array(Kv_lst)
        #=======================================================================
            
            data = [self.results_u[0],self.results_v[0],self.results_z[0]] #, Kv_arr]
            dateTimes1 = [self.time[0], self.time[0], self.time[0]] #, self.time[0]]
            varnames = ["East Velocity", "North Velocity", "Vert. Velocity"] #, "Eddy Viscosity"]
            limits = [[-0.6,0.6],[-0.6,0.6],[-0.06,0.06]] #, None]
            ylabels1 = [r'E Vel. [$\mathsf{m\cdot s^{-1}}$]',\
                        r'N Vel. [$\mathsf{m\cdot s^{-1}}$]', \
                        r'U Vel. [$\mathsf{m\cdot s^{-1}}$]']#, \
                        #r'$K_\mathsf{\nu}$ [$\mathsf{ m^2\cdot s^{-1}}$]']
            display_data.display_mixed_subplot(dateTimes1 = dateTimes1, data = data, varnames = varnames, ylabels1 = ylabels1, \
                                       yday = True, limits = limits, sharex = True, fontsize = 18)
        # for vel in velocities

    def readPcAdpTemp(self):
        [time, temp]  = self.pcadp.readTemp()
        self.results_temp = temp
        self.time =time
        return time, temp  

    def readRawBinADCP(self):
        try:
            with rdradcp.Timer.Timer() as t:
                [self.adcp, self.cfg, self.ens, self.hdr] = \
                    readRawBinADCP(self.path + '/' + self.filename, 1, [700, 239800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no') \

                #===============================================================
                # [self.adcp, self.cfg, self.ens, self.hdr] = \
                #     readRawBinADCP(self.path + '/' + self.filename, 1, [10000, 19800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no')
                #===============================================================
            self.adcp.goodbins = numpy.int(self.adcp.depth[0][500] - self.adcp.config.bin1_dist)+2
            self.goodbins=self.adcp.goodbins
           
        finally:
            print(('Read took %.03f sec.' % t.interval))
    
    @staticmethod
    def eddy_viscosity(V, W, dz, deltav = False):
        '''
        THIS IS WRONG IT NEEDS COVARIANCE AT THE spot of teh timeseries of a pulse
        NOW IT IS covariace on the veritcal space instead a time series.
        -------------------------------------
        Depth averagededdy viscosity 
        Using the method from: 
                    W. J. Shaw and J. H. Trowbridge 
                    The Direct Estimation of Near-Bottom Turbulent Fluxes in the Presence of Energetic Wave Motions,
                    J. Atmos. Oceanic Technol., 2001
                    
                    numpy.cov(a, b) retuns
                        |cov(a,a)  cov(a,b)|
                        |cov(a,b)  cov(b,b)|
        '''
        Kv = []
        V_T  = numpy.transpose(V)# The transposed
        grad = numpy.absolute(numpy.gradient(numpy.mean(V_T, axis = 0),dz))
                
        for k in range(0, len(V)):
            if deltav:
                "THIS NEEDS TO BE VERIFIED"
                dV = []
                dW = []
                for i in range(0, len(V[k])-1):
                     dV.append((V[k][i]-V[k][i+1]))
                dV.append((V[i]-V[i+1]))
                 
                for i in range(0, len(W[k])-1):
                     dW.append((W[k][i]-W[k][i+1]))
                dW.append((W[k][i]-W[k][i+1]))
                 
                covariance = numpy.absolute(numpy.cov(numpy.array(dV),numpy.array(dW))[0][1])
                Kv.append(covariance/grad[k])
            else:
                covariance = numpy.absolute(numpy.cov(V[k],W[k])[0][1])
                Kv.append(covariance/grad[k])
        return Kv

    def select_data_dates(self, smooth = False, dateint = None, velprofile = False, savetofile = False, \
                          delft3d=False, d3d_int=3):

        if dateint == None:
            date = self.date 
        else:
            date = dateint
            
        if date == None:
            print("No selection is possible 'date' is none")
            self.results_u = self.adcp.east_vel
            self.results_v = self.adcp.north_vel
            self.results_z = self.adcp.vert_vel
            
            self.time = self.adcp.mtime[0, :]
            return

        # select the data by dates
        dt = datetime.datetime.strptime(date[0], "%y/%m/%d %H:%M:%S")
        start_num = dates.date2num(dt)
        dt = datetime.datetime.strptime(date[1], "%y/%m/%d %H:%M:%S")
        end_num = dates.date2num(dt)

        if self.adcp is not None:
            time1, evel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.east_vel)
            time2, nvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.north_vel)
            time3, zvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.vert_vel)
        else:
            time1, evel = plot_ADCP_velocity.select_dates(start_num, end_num, self.time[0, :], self.results_u)
            time2, nvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.time[0, :], self.results_v)
            time3, zvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.time[0, :], self.results_z)

        evel = numpy.array(evel)
        nvel = numpy.array(nvel)
        zvel = numpy.array(zvel)
        time1 = numpy.array(time1)
        time2 = numpy.array(time2)
        time3 = numpy.array(time3)
        
        evel[numpy.isnan(evel)] = 0  # set to zero NAN values
        evel[numpy.isinf(evel)] = 0  # set to 0 infinite values
        nvel[numpy.isnan(nvel)] = 0  # set to zero NAN values
        nvel[numpy.isinf(nvel)] = 0  # set to 0 infinite values
        zvel[numpy.isnan(zvel)] = 0  # set to zero NAN values
        zvel[numpy.isinf(zvel)] = 0  # set to 0 infinite values

        if savetofile or delft3d:
            if delft3d:
                #for tm in time1:
                #dtime = datetime.datetime.strptime(tm, "%y/%m/%d %H:%M:%S")
                #st_strg = dtime.strftime("%y%m%d")
                dt0 = (time1[0][2] - time1[0][1])             #days
                d3d_int_days = d3d_int/24./60.                #convert into days
                print("dt0=%f" % dt0)
                if abs(d3d_int_days - dt0) > 1e-5: 
                    ratio = (time1[0][2]-time1[0][1]) / d3d_int_days
                    newlen = round(ratio* len(time1[0])) 
                    rtime1 = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    rtime2 = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    rtime3 = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    r_evel = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    r_nvel = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    r_uvel = numpy.zeros((len(time1),newlen), dtype = numpy.float)
                    for bin in range(0, len(time1)):
                        [r_evel[bin],rtime1[bin]] = utools.interpolate_data.interpolateData(d3d_int_days, evel[bin], time1[bin])
                        [r_nvel[bin],rtime2[bin]] = utools.interpolate_data.interpolateData(d3d_int_days, nvel[bin], time2[bin])
                        [r_uvel[bin],rtime3[bin]] = utools.interpolate_data.interpolateData(d3d_int_days, zvel[bin], time3[bin])
                else:
                    r_evel = evel
                    r_nvel = nvel
                    r_uvel = zvel
                    rtime1 = time1
                    rtime2 = time2
                    rtime3 = time3
                        
                writebins.writeBins(rtime1, r_evel, self.path, "Cell3-Aug2013-eastVel.tim", delft3d)
                writebins.writeBins(rtime2, r_nvel, self.path, "Cell3-Aug2013-northVel.tim", delft3d)
                writebins.writeBins(rtime3, r_uvel, self.path, "Cell3-Aug2013-vertVel.tim", delft3d)    
            #endif
            
            writebins.writeBins(time1, evel, self.path, "eastVel.csv")
            writebins.writeBins(time2, nvel, self.path, "northVel.csv")
            writebins.writeBins(time3, zvel, self.path, "vertVel.csv")


        # smoothfit requires list so we need to convert velocities back to list
        if smooth:
            span_window = utools.windows.window_hour
            smooth_window = utools.windows.windows[1]
            results_u = []
            results_v = []
            results_z = []
            i = 0
            for ev in evel:
                results_u.append(smooth.smoothfit(time1[i], ev.tolist(), span_window, smooth_window)['smoothed'])
                i += 1
            # end for
            i = 0
            for nv in nvel:
                results_v.append(smooth.smoothfit(time2[i], nv.tolist(), span_window, smooth_window)['smoothed'])
                i += 1
            for zv in zvel:
                results_z.append(smooth.smoothfit(time3[i], zv.tolist(), span_window, smooth_window)['smoothed'])
                i += 1
            # end for
            if dateint == None:
                self.results_u = numpy.array(results_u)
                self.results_v = numpy.array(results_v)
                self.results_z = numpy.array(results_z)
            else:
                tresults_u = numpy.array(results_u)
                tresults_v = numpy.array(results_v)
                tresults_z = numpy.array(results_z)
        else:
            if dateint == None:
                self.results_u = evel
                self.results_v = nvel
                self.results_z = zvel
            else:
                tresults_u = evel
                tresults_v = nvel
                tresults_z = zvel
        # end if
        if dateint == None:
            self.time = time1
        if self.adcp is not None:
            self.results_temp =  self.adcp.temperature
   
        if dateint != None:
            return [time1,  tresults_u, tresults_v, tresults_z] 
        
        if velprofile:
            
        #=======================================================================
        #     angles_from_N = {'OH':143, 'EmbC':53,'Cell3':43}      # for counter clockwise  
        #                                                   #+ direction is TO outer harbour and to Lake Ontario    
        #     num_segments = 10
        # 
        #     Theta = angles_from_N[self.locname] #35.1287  # degrees
        #     tet = 2 * math.pi * Theta / 360
        #     up, vp = plot_ADCP_velocity.rotation_transform(self.results_u, self.results_v, tet, clockwise=False)
        #     dz = 1
        #     vp_T = numpy.transpose(vp)
        #     z_T  = numpy.transpose(self.results_z)
        #     Kv_lst = []
        #     for i in range(0, len(z_t)):
        #         Kv = Hydrodynamic.eddy_viscosity(vp_T[i], z_T[i], dz)
        #         Kv_lst.append(Kv)
        #     Kv_arr = numpy.array(Kv_lst)
        #=======================================================================
            
            data = [self.results_u[5],self.results_v[5],self.results_z[3]]#, Kv_arr]
            dateTimes1 = [self.time[0], self.time[0], self.time[0]]#, self.time[0]]
            varnames = ["East Velocity", "North Velocity", "Vert. Velocity"]#, "Eddy Viscosity"]
            limits = [[-0.8,0.8],[-0.8,0.8],[-0.06,0.06]] #, [0, 1e-3]]
            
            ylabels1 = [r'E Vel. [$\mathsf{m\cdot s^{-1}}$]',\
                        r'N Vel. [$\mathsf{m\cdot s^{-1}}$]', \
                        r'U Vel. [$\mathsf{m\cdot s^{-1}}$]'] #, \
                        #r'$K_\mathsf{\nu}$ [$\mathsf{ m^2\cdot s^{-1}}$]']
                                    
            display_data.display_mixed_subplot(dateTimes1 = dateTimes1, data = data, varnames = varnames, ylabels1 = ylabels1, \
                                       yday = True, limits = limits, sharex = True, fontsize = 18)


    def plot_uvt(self, smooth = False, sel_dates = True):
        self.select_data_dates(smooth)
        plot_ADCP_velocity.plot_temp_u_v(adcp, time1, results_u, time2, results_v, TH_dateTimeArr, TH_resultsArr, interp = interp)


    def plot_wavelet(self, smooth = False, sel_dates = True):
        self.select_data_dates(smooth)
        scaleunit = 'day'

        # must down sample here
        if resample:
            if len(TH_dateTimeArr) != len(time1):
                res_time1 = scipy.ndimage.interpolation.zoom(time1[bin], float(len(TH_dateTimeArr[tlogno])) / float(len(time1[bin])))
                res_u = scipy.ndimage.interpolation.zoom(evel[bin], float(len(TH_dateTimeArr[tlogno])) / float(len(time1[bin])))
                res_v = scipy.ndimage.interpolation.zoom(nvel[bin], float(len(TH_dateTimeArr[tlogno])) / float(len(time1[bin])))

                print("plot vel wavelet")
                plot_ADCP_vel_FFT.plot_velocity_wavelet_spectrum(res_time1, res_u, scaleunit = scaleunit)

                print("plot temp wavelet")
                plot_ADCP_vel_FFT.plot_velocity_wavelet_spectrum(TH_dateTimeArr[tlogno][1:], TH_resultsArr[tlogno][1:], scaleunit = scaleunit)

                print("plot cross wavelet")
        else:
            print("plot cross wavelet")
            plot_ADCP_vel_FFT.plot_cross_spectogram_w_T(time1[bin], results_u[bin], results_v[bin], TH_dateTimeArr[tlogno], TH_resultsArr[tlogno], scaleunit = scaleunit)
        # endif

    def plot_FFT(self, locname, bins, smooth = False, tunits = "day", funits = "Hz", \
                 log = False, grid = False, type = 'power', withci = True, sel_dates = True):
        window = "hanning"

        if sel_dates:
            self.select_data_dates(smooth)

        showLevels = False
        detrend = False
        show = False
        draw = False  # do not draw the series

        fftsa = numpy.zeros(len(bins), dtype = numpy.object)
        data = []
        ci05 = []
        ci95 = []
        freq = []
        names = []

        for i in range(0, len(bins)):
            fftsa[i] = fftGraphs.FFTGraphs(path = None, file1 = None, file2 = None, show = None, \
                                           data = [self.time[bins[i]], self.results_u[bins[i]] + 1j * self.results_v[bins[i]]], data1 = None)

            date1st = True
            names.append(locname + ": bin " + str(bins[i]))
            fftsa[i].doSpectralAnalysis(showLevels, draw, tunits, window, self.num_segments, filter = None, log = log, date1st = date1st)
            if type == 'power':
                data.append(fftsa[i].power)
            else:
                data.append(fftsa[i].mx)
            ci05.append(fftsa[i].x05)
            ci95.append(fftsa[i].x95)
            freq.append(fftsa[i].f)
        # end for
        lake_name = ""
        if type == 'power':
            ylabel = "Spectral Power [$m^2 s^{-4} Hz^{-1}$]"
        else:
            ylabel = "Velocity [$m s^{-1}$]"
        if withci:
            ci = [ci05, ci95]
        else:
            ci = None
        fftGraphs.plotSingleSideSpectrumFreqMultiple(lake_name, names, data, freq, ci, type, \
                                                             self.num_segments, funits, y_label = ylabel, title = None, \
                                                             log = log, fontsize = 20, tunits = tunits)


    #---------------------------------------------------------------------------------------------------------------
    def Vel_profiles_2016_ADCP(self, date, adptype='OH', firstbin=1, interval=1,
                               ADCPprof=False, grad_valong=False, clockwise=False, adcp=None):
        '''
         Based on the rotation transformation with counter clockwise
          (+) is going to the bay
          (-) in coming out of the bay into Outer harbour

        '''
        print("Processing graphs for ADCP: %s" % adptype)

        matplotlib.mathtext.SHRINK_FACTOR = 0.9


        # + direction is TO outer harbour and to Lake Ontario
        num_segments = self.num_segments

        factor = 1.  # 1 hour in days
        factor = 6.  # 10 minutes  in days

        dt = 1. / 24. / factor
        labels = [' velocity [m/s]']

        self.get_Velocities(adptype, date, num_segments, adcp=adcp)

        if adptype == 'OH':
            Yrange = [0.0, 7.5]
            Xrange = [-0.70, 0.70]
        elif adptype == "EmbC":
            Yrange = [0.0, 5.5]
            Xrange = [-0.70, 0.70]
        elif adptype == "TI":
            Xrange = [-0.7, 0.9]
            Yrange = [0.0, 5.5]
        elif adptype == "EGap":
            Xrange = [-0.7, 0.9]
            Yrange = [0.0, 5.5]
        elif adptype == "WGap":
            Xrange = [-0.7, 0.9]
            Yrange = [0.0, 5.5]

        if clockwise:
            angles_from_N = {'OH':37, 'EmbC':127,'TI':165, 'EGap': 129, 'WGap': 56}    # for clockwise
            # + sign is into Cell3   and to Cherry beach and into the TI channels and out from the EGap
        else:
            # for counter clockwise  and vp the along dir
            angles_from_N = {'OH': 143, 'EmbC': 53, 'TI': 15, 'EGap': 51, 'WGap': 124}

        # reproject on the direction of thde flow
        Theta = angles_from_N[adptype]  # 35.1287  # degrees
        tet = 2 * math.pi * Theta / 360
        up, vp = self.rotation_transform(tet, clockwise=clockwise)

        # if counterclockwise vp is the along the channel veloccity and + (PLUS) is toward the lake
        # if clockwise up is the along the channel veloccity and + (PLUS) is toward the embayment
        print(
            "htime-0=%f, htime0-1=%f, size=%d" % (self.time[0][0], self.time[0][len(self.time[0]) - 1], len(self.time[0])))

        upw_str = '16/07/20 12:00:00'
        dt = datetime.datetime.strptime(upw_str, "%y/%m/%d %H:%M:%S")
        upwell = dates.date2num(dt)

        dw_str = '16/07/26 12:00:00'
        dt = datetime.datetime.strptime(dw_str, "%y/%m/%d %H:%M:%S")
        downwell = dates.date2num(dt)

        stab_str = '16/08/05 12:00:00'
        dt = datetime.datetime.strptime(stab_str, "%y/%m/%d %H:%M:%S")
        stable = dates.date2num(dt)
        profiledates = [stable, upwell, downwell]

        # profiledates = None
        revert = True  # True 0 bin bottom of graph because it is reversing limits in the routine below

        diff = 0 #self.goodbins - len(vp)
        if diff == 0:
            timearr = self.time[:]
            velarr = vp[:]
        else:
            timearr = self.time[:diff]
            velarr = vp[:diff]
        utools.display_data.display_avg_vertical_temperature_profiles_err_bar_range([timearr], [velarr],
                                                                                    startdeptharr=[0],
                                                                                    profiledates=profiledates, doy=True,
                                                                                    revert=revert,
                                                                                    legendloc=4,
                                                                                    xlabel=' Velocity [$m s^{-1}$]',
                                                                                    errbar=True, rangebar=True,
                                                                                    debug=False)

        if grad_valong or ADCPprof:
            # 3) Mixed water, air ,img data
            # ToDO: Add short and long radiation
            print("Grad Valong Start display mixed subplots  ")
            grad = []
            V_T = numpy.transpose(vp)
            U_T = numpy.transpose(up)

            # TODO get the config value
            dz = adcp.config.cell_size
            for i in range(0, len(self.time[0])):
                dvdz = []
                for j in range(0, len(self.time) - 1):
                    W1 = math.sqrt(V_T[i][j + 1] ** 2 + U_T[i][j + 1] ** 2)
                    W2 = math.sqrt(V_T[i][j] ** 2 + U_T[i][j] ** 2)
                    dvdz.append((W1 - W2) / dz / 2)

                dvdz.insert(0, (0 + dvdz[0]) / 2.)  # insert fist item to match the velocity array size
                grad.append(dvdz)

            gradarr = numpy.array(grad)
            grad_T = numpy.transpose(gradarr)

            if adptype == "WGap":
                dd = 1  # for Cell 3
            elif adptype == "EmbC":
                dd = 2
            elif adptype == "OH":
                dd = 3
            else:
                dd = 1

            dateTimes3 = [self.time[:-dd], self.time[:-dd], self.time[:-dd]]
            ylabels3 = ["Depth [m]", "Depth [m]", "Depth [m]"]

            imgs = [vp[:-dd], up[:-dd], grad_T[:-dd]]
            # imgs = [self.results_u[:-dd], self.results_v[:-dd]] #Harbour first

            if adptype == 'OH' or adptype == 'EGap':
                t11 = ['0', '3', '6', '10']
                t12 = [10, 6, 3, 0]

                t21 = ['0', '3', '6', '10']
                t22 = [10, 6, 3, 0]

                t31 = ['0', '3', '6', '10']
                t32 = [10, 6, 3, 0]

                # maxdepth = [9, 27] # Harbour first
                maxdepth = [10, 10, 10]
                # firstlogdepth = [0, 3] Harbour first
                firstlogdepth = [0, 0, 0]
                mindepths = [0, 0, 0]
                maxtemp = [0.2, 0.2, 0.05]
                mintemps = [-0.2, -0.2, -0.05]
            elif adptype == 'WGap':
                t11 = ['0', '3', '9', '9']
                t12 = [9, 6, 3, 0]

                t21 = ['0', '3', '6', '9']
                t22 = [9, 6, 3, 0]

                t31 = ['0', '3', '6', '9']
                t32 = [9, 6, 3, 0]

                # maxdepth = [9, 27] # Harbour first
                maxdepth = [9, 9, 9]
                # firstlogdepth = [0, 3] Harbour first
                firstlogdepth = [0, 0, 0]
                mindepths = [0, 0, 0]
                maxtemp = [0.2, 0.2, 0.05]
                mintemps = [-0.2, -0.2, -0.05]
            elif adptype == "EmbC":
                t11 = ['0', '1.5', '3.0', '4.5']
                t12 = [4.5, 3, 1.5, 0]

                t21 = ['0', '1.5', '3.0', '4.5']
                t22 = [4.5, 3, 1.5, 0]

                t31 = ['0', '1.5', '3.0', '4.5']
                t32 = [4.5, 3, 1.5, 0]

                # maxdepth = [9, 27] # Harbour first
                maxdepth = [4.5, 4.5, 4.5]
                # firstlogdepth = [0, 3] Harbour first
                firstlogdepth = [0, 0, 0]
                mindepths = [0, 0, 0]
                maxtemp = [0.2, 0.2, 0.05]
                mintemps = [-0.2, -0.2, -0.05]
            elif adptype == "TI":
                t11 = ['0', '2', '4', '6.5']
                t12 = [6.5, 4, 2, 0]

                t21 = ['0', '2', '4', '6.5']
                t22 = [6.5, 4, 2, 0]

                t31 = ['0', '2', '4', '6.5']
                t32 = [6.5, 4, 2, 0]

                # maxdepth = [9, 27] # Harbour first
                maxdepth = [6.5, 6.5, 6.5]
                # firstlogdepth = [0, 3] Harbour first
                firstlogdepth = [0, 0, 0]
                mindepths = [0, 0, 0]
                maxtemp = [0.2, 0.2, 0.05]
                mintemps = [-0.2, -0.2, -0.05]

            tick = [[t11, t12], [t21, t22], [t31, t32]]

            # limits = [None, None, [-1, 10], None, None ] <= this datesscrews up the tickers
            limits = None

            if grad_valong:
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
                imgs = [self.results_u[:-dd], self.results_v[:-dd], self.results_z[:-dd]]
                utools.display_data.display_mixed_subplot(dateTimes1=[], data=[], varnames=[], ylabels1=[],
                                                          dateTimes2=[], groups=[], groupnames=[], ylabels2=[],
                                                          dateTimes3=dateTimes3, imgs=imgs, ylabels3=ylabels3, ticks=tick,
                                                          maxdepths=maxdepth,
                                                          mindepths=mindepths, mintemps=mintemps, firstlogs=firstlogdepth,
                                                          maxtemps=maxtemp,
                                                          fnames=None, revert=True, custom=None, maxdepth=None, tick=None,
                                                          firstlog=None, yday=True,
                                                          title=False, grid=False, limits=limits, sharex=True, fontsize=18,
                                                          group_first=False, interp=2,
                                                          cblabel=["$V_E$ [$m s^{-1}$]", "$V_N$ [$m s^{-1}$]",
                                                                   "$V_U$ [$m s^{-1}$]"])



    #-------------------------------------------------------------------------------------------------------------------
    def plot_Temp_FFT(self, locname, smooth = False, tunits = "day", funits = "Hz",
                 log = False, grid = False, type = 'power', withci = True, sel_dates = True):
        window = "hanning"

        if sel_dates:
            self.select_data_dates(smooth)

        showLevels = False
        detrend = False
        show = False
        draw = False  # do not draw the series

        
        
        fftsa = fftGraphs.FFTGraphs(path = None, file1 = None, file2 = None, show = None,
                                       data = [self.time, self.results_temp], data1 = None)

        date1st = True
        names=locname
        fftsa.doSpectralAnalysis(showLevels, draw, tunits, window, self.num_segments, filter = None, log = log, date1st = date1st)
        if type == 'power':
            data=fftsa.power
        else:
            data=fftsa.mx
        ci05=[fftsa.x05]
        ci95=[fftsa.x95]
        freq=[fftsa.f]
        lake_name = ""
        if type == 'power':
            ylabel = "Spectral Power [$m^2 s^{-4} Hz^{-1}$]"
        else:
            ylabel = 'Temperature [$^\circ$C]'
            
        
        if withci:
            ci = [ci05, ci95]
        else:
            ci = None
            
        fftGraphs.plotSingleSideSpectrumFreqMultiple(lake_name, [names], [data], freq, ci, type, \
                                                             self.num_segments, funits, y_label = ylabel, title = None, \
                                                             log = log, fontsize = 20, tunits = tunits)
        
    def rotation_transform(self, tet, clockwise=False, u=None, v=None):
        if self.results_u is not None and u is None:
            return plot_ADCP_velocity.rotation_transform(self.results_u, self.results_v, tet, clockwise=clockwise)
        elif u is not None:
            return plot_ADCP_velocity.rotation_transform(u, v, tet, clockwise=clockwise)
        print('Velocity values not read!')
        raise Exception("Velocity values not read!")

    @staticmethod
    def resample(time, up, dt, bin):
        dt0 = (time[0][2] - time[0][1])
        print("dt0=%f" % dt0)
        winlen = int(dt / dt0)
        #y, dy = smooth.smooth(up[bin], winlen)
        y = up[bin]
        rup = scipy.ndimage.interpolation.zoom(numpy.real(y), float(dt0 / dt))
        rtime = scipy.ndimage.interpolation.zoom(time[0], float(dt0 / dt))
        return rtime, rup
    
    @staticmethod
    def mov_average(time, up, dt, bin):
        dt0 = (time[0][2] - time[0][1])
        print("dt0=%f" % dt0)
        window = int(dt / dt0)
        weigths = numpy.repeat(1.0, window)/window
        #including valid will REQUIRE there to be enough datapoints.
        #for example, if you take out valid, it will start @ point one,
        #not having any prior points, so itll be 1+0+0 = 1 /3 = .3333
        smas = numpy.convolve(up[bin], weigths, 'valid')
        rtime = scipy.ndimage.interpolation.zoom(time[0], float(dt0 / dt))
        return rtime, smas # as a numpy array

    @staticmethod
    def calculateVector():
        curdir = None
        curspd = None
        return curdir, curspd

    @staticmethod
    def draw_windrose(wd, ws, type, loc = 'best', fontsize = 10, unit = "[m/s]", angle = 67.5):
        return plot_ADCP_velocity.draw_windrose(wd, ws, type, loc , fontsize, unit, angle = angle)
    
    def printVelProfiles(self, name, profiles, datetimes, firstbin, interval, Xrange= None, Yrange=None, save=False, dzdt = None):
        for i in range(0, len(profiles)):
            plot_ADCP_velocity.print_profile(name, profiles[i], datetimes[i], firstbin, \
                                             Xrange= Xrange, Yrange=Yrange, interval = interval, spline=True, save = save, dzdt = dzdt[i])
            
    def display_subplots(self, date, dateIg, dataarr, dnames = None, yday = None, tick = None, legend = None,
                         hourgrid = False, img=False, cbrange = [-1,1], maxdepth=5.0):
        plot_ADCP_velocity.display_subplots(date, dateIg, dataarr, dnames = dnames, yday = yday,
                                            tick = tick, legend = legend, hourgrid = hourgrid, \
                                            img=img, cbrange = cbrange, maxdepth= maxdepth, minorlabel=False)  
        
    def plot_FFT_V_T_WL(self, time, V, Ttime, T, WL, WTime, scale='log', drawslope=False, plotTemp=False ):
        plot_ADCP_vel_FFT.plot_FFT_Three_V_T_WL(time, V, Ttime, T, WL, WTime, scale = scale,
                                                drawslope=drawslope, plotTemp=plotTemp)
    
    @staticmethod
    def get_data_from_file(fname, timeinterv, rpath = None):
        span_window = utools.windows.window_hour
        smooth_window = utools.windows.windows[1]
        return readTempHoboFiles.get_data_from_file(fname, span_window, smooth_window, timeinterv , rpath = rpath)
    
    @staticmethod
    def estimate_vel(ua, a, y):
        '''
        Engineering, 2013, 5, 933-942
        Published Online December 2013 (http://www.scirp.org/journal/eng)
        http://dx.doi.org/10.4236/eng.2013.512114
        Open Access ENG
        Power Law Exponents for Vertical Velocity Distributions in Natural Rivers
        Hae-Eun Lee, Chanjoo Lee, Youg-Jeon Kim,, Ji-Sung Kim2, Won Kim2
        
        estimate the velofity at a 
        @param y: distance form bottom when
        @param ua : konw velocity ar
        @param a: known distance above the bottom
        
        '''
        u = ua * (y/a)**(1/1.2)
        return u
    
    @staticmethod
    def writeVeltoCSV(fname, data, append =False, delft3d=False):
        if delft3d:
            dtime = datetime.datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
            st_strg = dtime.strftime("%y%m%d")
            [dates, depths] = self.read_press_corr_file(path, fn)
            
            dt0 = (dates[2] - dates[1])           #days
            print("dt0=%f" % dt0)
            dt = step_min/24./60.                 #convert into days
            winlen = int(dt / dt0)
            if not isinstance(depths, np.ndarray):
                depths = np.array(depths)
                dates = np.array(dates)
    
            weigths = np.repeat(1.0, winlen)/winlen
            #including valid will REQUIRE there to be enough datapoints.
            #for example, if you take out valid, it will start @ point one,
            #not having any prior points, so itll be 1+0+0 = 1 /3 = .3333
            rdepths = np.convolve(depths, weigths, 'valid')
            rtime = sp.ndimage.interpolation.zoom(dates, float(dt0 / dt))
        plot_ADCP_velocity.writeVeltoCSV(fname, data, append)
    
    @staticmethod    
    def ReadVolfromCSV(fname, level):
        return plot_ADCP_velocity.ReadVolfromCSV(fname, level)
    
    @staticmethod
    def read_Loboviz2datefile(path, filename):
        ifile = open(path + '/' + filename, 'rt', encoding="ISO-8859-1")
        reader = csv.reader(ifile, delimiter=',', quotechar='"')
        rownum = 0
        temp = []
        east = []
        north = []
        up = []
        dateTime = []
        datenum = []
        printHeaderVal = True

        for row in reader:
            # skip comments with @
            if row[0][:1] == '@':
                continue

            # Save header row.
            strg = row[0]
            if strg[0] == '#':
                header = row
                if printHeaderVal == True:
                    colnum = 0
                    for col in row:
                        print('%-8s: %s' % (header[colnum], col))
                        colnum += 1
            else:
                newrow = True
                for ir in range(0, len(row)):
                    if ir == 0:
                        dateTime.append(row[ir])
                    elif ir == 1:
                        datenum.append(row[ir])
                    elif ir == 2:
                        temp.append(float(row[ir]))
                    elif header[ir][:4] == 'east':
                        if newrow:
                            east.append([])
                        east[rownum - 1].append(float(row[ir]))
                        newrow = False
                    elif header[ir][:4] == 'nort':
                        if newrow or len(north) <= rownum - 1:
                            north.append([])
                        north[rownum - 1].append(float(row[ir]))
                        newrow = False
                    elif header[ir][:2] == 'up':
                        if newrow or len(up) <= rownum - 1:
                            up.append([])
                        up[rownum - 1].append(float(row[ir]))
                        newrow = False
            rownum += 1
        ifile.close()
        return [temp, dateTime, datenum, east, north, up]
  
        