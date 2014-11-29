'''
Created on Nov 20, 2014

@author: bogdan
'''
from rdradcp.readRawBinADCP import readRawBinADCP
import rdradcp.plot_ADCP_vel_FFT as plot_depth_averaged_analysis
import rdradcp.plot_ADCP_vel_FFT as plot_velocity_wavelet_spectrum
import rdradcp.plot_ADCP_vel_FFT as plot_cross_spectogram_w_T
import rdradcp.plot_ADCP_vel_FFT as plot_FFT_twinx_W_T
import rdradcp.plot_ADCP_vel_FFT as plot_fft_analysis
import rdradcp.plot_ADCP_velocity as plot_ADCP_velocity
import rdradcp.writebins as writebins
import rdradcp.Timer

import pcadcp.rdTxtPcAdp as rdTxtPcAdp

import ufft.spectral_analysis
import ufft.FFTGraphs as fftGraphs
import utools.windows

import datetime
import matplotlib.dates as dates
import numpy

class Hydrodynamic(object):
    '''
    classdocs
    '''

    def __init__(self, name, path, num_segments, tinterv, ctlname = None, filenames = None):
        '''
        Constructor
        '''
        self.path = path
        self.filename = name
        self.num_segments = num_segments
        self.adcp = None
        self.ens = None
        self.cfg = None
        self.hdr = None
        self.results_u = None
        self.results_v = None
        self.results_z = None
        self.time = None

        if tinterv != None:
            self.date = [tinterv[0], tinterv[1]]

        self.ctlname = ctlname
        self.filenames = filenames

        if ctlname != None:
            self.pcadp = rdTxtPcAdp.PcADP(name, path, ctlname, filenames, num_segments, tinterv)
        else:
            self.pcadp = None

    def readPcAdpVel(self):
        velocities, timevec = self.pcadp.readVel()
        self.results_u = velocities['ve']
        self.results_v = velocities['vn']
        self.results_z = velocities['vu']
        self.time = timevec
        # for vel in velocities


    def readRawBinADCP(self):
        try:
            with rdradcp.Timer.Timer() as t:
                [self.adcp, self.cfg, self.ens, self.hdr] = \
                    readRawBinADCP(self.path + '/' + self.filename, 1, [700, 239800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no')
                #===============================================================
                # [self.adcp, self.cfg, self.ens, self.hdr] = \
                #     readRawBinADCP(self.path + '/' + self.filename, 1, [10000, 19800], 'info', 'yes', 'baseyear', 2000, 'despike', 'yes', 'debug', 'no')
                #===============================================================
        finally:
            print('Read took %.03f sec.' % t.interval)

    def select_data_dates(self, smooth = False):

        if self.date == None:
            print "No selection is possible 'date' is none"
            self.results_u = self.adcp.east_vel
            self.results_v = self.adcp.north_vel
            self.time = self.adcp.mtime[0, :]
            return

        # select the data by dates
        dt = datetime.datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
        start_num = dates.date2num(dt)
        dt = datetime.datetime.strptime(self.date[1], "%y/%m/%d %H:%M:%S")
        end_num = dates.date2num(dt)

        time1, evel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.east_vel)
        time2, nvel = plot_ADCP_velocity.select_dates(start_num, end_num, self.adcp.mtime[0, :], self.adcp.north_vel)

        evel = numpy.array(evel)
        nvel = numpy.array(nvel)
        evel[numpy.isnan(evel)] = 0  # set to zero NAN values
        evel[numpy.isinf(evel)] = 0  # set to 0 infinite values
        nvel[numpy.isnan(nvel)] = 0  # set to zero NAN values
        nvel[numpy.isinf(nvel)] = 0  # set to 0 infinite values

        writebins.writeBins(time1, evel, self.paths['path_out'], "eastVel.csv")
        writebins.writeBins(time2, nvel, self.paths['path_out'], "northVel.csv")


        # smoothfit requires list so we need to convert velocities back to list
        if smooth:
            span_window = utools.interp_windows.window_hour
            smooth_window = utools.interp_windows.windows[1]
            results_u = []
            results_v = []
            i = 0
            for ev in evel:
                results_u.append(smooth.smoothfit(time1[i], ev.tolist(), span_window, smooth_window)['smoothed'])
                i += 1
            # end for
            i = 0
            for nv in nvel:
                results_v.append(smooth.smoothfit(time2[i], nv.tolist(), span_window, smooth_window)['smoothed'])
                i += 1
            # end for
            self.results_u = numpy.array(results_u)
            self.results_v = numpy.array(results_v)
        else:
            self.results_u = evel
            self.results_v = nvel
        # end if
        self.time = time1

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
                plot_ADCP_vel_FFT.plot_cross_spectogram_w_T(res_time1, res_u, res_v, res_time1, TH_resultsArr[tlogno], scaleunit = scaleunit)
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


