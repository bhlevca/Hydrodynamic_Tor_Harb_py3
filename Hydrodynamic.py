'''
Created on Nov 20, 2014

@author: bogdan
'''
from rdradcp.readRawBinADCP import readRawBinADCP
import rdradcp.plot_ADCP_vel_FFT as plot_ADCP_vel_FFT
import rdradcp.plot_ADCP_velocity as plot_ADCP_velocity
import rdradcp.writebins as writebins
import rdradcp.Timer

import pcadcp.rdTxtPcAdp as rdTxtPcAdp

import ufft.spectral_analysis
import ufft.FFTGraphs as fftGraphs
import utools.windows

import utools.display_data as display_data
import ufft.smooth as smooth
import utools.readTempHoboFiles as readTempHoboFiles

import datetime,math
import matplotlib.dates as dates
import numpy, scipy

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
                    readRawBinADCP(self.path + '/' + self.filename, 1, [700, 239800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no')

                #===============================================================
                # [self.adcp, self.cfg, self.ens, self.hdr] = \
                #     readRawBinADCP(self.path + '/' + self.filename, 1, [10000, 19800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no')
                #===============================================================
            self.adcp.goodbins = numpy.int(self.adcp.depth[0][500] - self.adcp.config.bin1_dist)+2
            self.goodbins=self.adcp.goodbins
           
        finally:
            print('Read took %.03f sec.' % t.interval)
    
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

    def select_data_dates(self, smooth = False, dateint = None, velprofile = False, savetofile = False):

        if dateint == None:
            date = self.date 
        else:
            date = dateint
            
        if date == None:
            print "No selection is possible 'date' is none"
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

        time1, evel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.east_vel)
        time2, nvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.north_vel)
        time3, zvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.vert_vel)

        evel = numpy.array(evel)
        nvel = numpy.array(nvel)
        zvel = numpy.array(zvel)
        evel[numpy.isnan(evel)] = 0  # set to zero NAN values
        evel[numpy.isinf(evel)] = 0  # set to 0 infinite values
        nvel[numpy.isnan(nvel)] = 0  # set to zero NAN values
        nvel[numpy.isinf(nvel)] = 0  # set to 0 infinite values
        zvel[numpy.isnan(zvel)] = 0  # set to zero NAN values
        zvel[numpy.isinf(zvel)] = 0  # set to 0 infinite values

        if savetofile:
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

                print "plot vel wavelet"
                plot_ADCP_vel_FFT.plot_velocity_wavelet_spectrum(res_time1, res_u, scaleunit = scaleunit)

                print "plot temp wavelet"
                plot_ADCP_vel_FFT.plot_velocity_wavelet_spectrum(TH_dateTimeArr[tlogno][1:], TH_resultsArr[tlogno][1:], scaleunit = scaleunit)

                print "plot cross wavelet"
        else:
            print "plot cross wavelet"
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
                                                             log = log, fontsize = 24, tunits = tunits)


    def plot_Temp_FFT(self, locname, smooth = False, tunits = "day", funits = "Hz", \
                 log = False, grid = False, type = 'power', withci = True, sel_dates = True):
        window = "hanning"

        if sel_dates:
            self.select_data_dates(smooth)

        showLevels = False
        detrend = False
        show = False
        draw = False  # do not draw the series

        
        
        fftsa = fftGraphs.FFTGraphs(path = None, file1 = None, file2 = None, show = None, \
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
                                                             log = log, fontsize = 24, tunits = tunits)
        
    def rotation_transform(self, tet, clockwise=False, u=None, v=None):
        if self.results_u is not None and u is None:
            return plot_ADCP_velocity.rotation_transform(self.results_u, self.results_v, tet, clockwise=clockwise)
        elif u is not None:
            return plot_ADCP_velocity.rotation_transform(u, v, tet, clockwise=clockwise)
        print 'Velocity values not read!'
        raise Exception("Velocity values not read!")

    @staticmethod
    def resample(time, up, dt, bin):
        dt0 = (time[0][2] - time[0][1])
        print "dt0=%f" % dt0
        winlen = int(dt / dt0)
        #y, dy = smooth.smooth(up[bin], winlen)
        y = up[bin]
        rup = scipy.ndimage.interpolation.zoom(numpy.real(y), float(dt0 / dt))
        rtime = scipy.ndimage.interpolation.zoom(time[0], float(dt0 / dt))
        return rtime, rup
    
    @staticmethod
    def mov_average(time, up, dt, bin):
        dt0 = (time[0][2] - time[0][1])
        print "dt0=%f" % dt0
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
            
    def display_subplots(self, date, dateIg, dataarr, dnames = None, yday = None, tick = None, legend = None, hourgrid = False, img=False, cbrange = [-1,1], maxdepth=5.0):    
        plot_ADCP_velocity.display_subplots(date, dateIg, dataarr, dnames = dnames, yday = yday, tick = tick, legend = legend, hourgrid = hourgrid, \
                                            img=img, cbrange = cbrange, maxdepth= maxdepth, minorlabel=False)  
        
    def plot_FFT_V_T_WL(self, time, V, Ttime, T, WL, WTime, scale = 'log', drawslope = False ):
        plot_ADCP_vel_FFT.plot_FFT_Three_V_T_WL(time, V, Ttime, T, WL, WTime, scale = scale, drawslope = drawslope)
    
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
        @param a: known distance abouve the bottom
        
        '''
        u = ua * (y/a)**(1/1.2)
        return u
    
    @staticmethod
    def writeVeltoCSV(fname, data, append =False):
        plot_ADCP_velocity.writeVeltoCSV(fname, data, append)
    
    @staticmethod    
    def ReadVolfromCSV(fname, level):
        return plot_ADCP_velocity.ReadVolfromCSV(fname, level)
    
    
  
        