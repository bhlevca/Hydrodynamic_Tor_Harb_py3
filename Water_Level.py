'''
Class handling the water level processing from HOPB pressure loggers.
Data is supposed to be processed with the right interval and in m from kPa
'''

import scipy as sp
import numpy as np
import math
import matplotlib.mlab as mlab
import csv
from scipy import fftpack
from datetime import datetime
import matplotlib.dates as dates
# local imports
import ufft.fft_utils as fft_utils
import wavelets.kCwt
import ufft.filters
import ufft.FFTGraphs as FFTGraphs
import utools.display_data as display_data
import ufft.smooth as smooth

class WaterLevelAnalysis(object):
    '''
    classdoc
    '''

    def __init__(self, path, filenames, numsegments, tinterv = None):
        self.path = path
        self.dict = filenames
        self.num_segments = numsegments
        if tinterv != None:
            self.date = [tinterv[0], tinterv[1]]
        else:
            self.date = None
    def getDict(self):
        return self.dict


    def read_press_corr_file(self, path, fname):
        '''
        :param path: -- path to the file location
        :param fname: -- file name
        '''

        if self.date != None:
            dtime = datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
            startt = dates.date2num(dtime)
            dtime = datetime.strptime(self.date[1], "%y/%m/%d %H:%M:%S")
            endt = dates.date2num(dtime)

        ifile = open(path + '/' + fname, 'rt')
        reader = csv.reader(ifile, delimiter = ',', quotechar = '"')
        dateTime = []
        depths = []

        for row in reader:
            try:
                time = float(row[0])

                if self.date != None:
                    if time < startt or time > endt:
                        continue
                dateTime.append(time)
                depths.append(float(row[1]))
            except:
                print("Error:read_press_corr_file")
            # end try

        ifile.close()

        return [dateTime, depths]

    def calculateFFT(self, Time, SensorDepth, tunits = 'sec', log = False):
        '''
        Clearly, it is difficult to identify the frequency components from looking at this signal;
        that's why spectral analysis is so popular.

        Finding the discrete Fourier transform of the noisy signal y is easy; just take the
        fast-Fourier transform (FFT).

        Compute the amplitude spectral density, a measurement of the amplitude at various frequencies,
        using module (abs)

        Compute the power spectral density, a measurement of the energy at various frequencies,
        using the complex conjugate (CONJ).

        nextpow2 finds the exponent of the next power of two greater than or equal to the window length (ceil(log2(m))), and pow2 computes the power. Using a power of two for the transform length optimizes the FFT algorithm, though in practice there is usually little difference in execution time from using n = m.
        To visualize the DFT, plots of abs(y), abs(y).^2, and log(abs(y)) are all common. A plot of power versus frequency is called a periodogram:
        @param Time : the time points
        @param SensorDepth: the depth data timeseries
        @param tunits:
        @param window: = 'blackman' #NO; 'bartlett' #OK; 'hamming' #OK; 'hann' #BEST default; flattop; gaussian; blackmanharris; barthann; bartlett;
        @param num_segments = 1 default represents tne number of non overlapping segments used for the overlapping Welch method.
                              The number for the total number of ssegments is M = 2* num_segments-1: For example if num_segments=2 => M=3
        @param filter: defaul = None to avoid filtering twice in a recursive method. filtes is of type fft.Filter

        @return: y             - detrended water levels
                 Time          - Time data points
                 fftx          - unique FFT values for the series
                 NumUniquePts  - size of the unique FFT values
                 mx            - the value of the single-sided FFT amplitude
                 f             - linear frequency vector for the mx points
        '''
        eps = 1e-3

        L = len(Time)


        # prepare for the amplitude spectrum analysis
        if tunits == 'day':
            factor = 86400
        elif tunits == 'hour':
            factor = 3600
        else:
            factor = 1

        dt_s = (Time[2] - Time[1]) * factor  # Sampling period [s]
        Fs = 1 / dt_s  # Samplig freq    [Hz]

        # nextpow2 = This function is useful for optimizing FFT operations, which are most efficient when sequence length is an exact power of two.
        #  does seem to affect the value of the amplitude
        # # NFFT = fft_utils.nextpow2(L)   # Next power of 2 from length of the original vector, transform length
        #
        NFFT = L
        # DETREND THE SIGNAL is necessary to put all signals oscialltions around 0 ant the real level in spectral analysis
        yd = sp.signal.detrend(SensorDepth)

        # Take fft, padding with zeros so that length(fftx) is equal to nfft
        fftx = fftpack.fft(yd, NFFT)  # DFT
        sFreq = np.sum(abs(fftx) ** 2) / NFFT
        sTime = np.sum(yd ** 2)

        # This is a sanity check
        assert abs(sFreq - sTime) < eps

        # What is power of the DFT and why does not show anything for us?
        # The FFT of the depth non CONJ shows good data except the beginning due to the fact
        # that the time series are finite.
        power = (np.abs(fftx)) ** 2  # Power of the DFT

        # TEST the Flat Top
        # amp = np.sqrt(2 * power / NFFT / Fw)
        # amp = 2 * np.abs(self.Wind_Flattop(fftx)) / NFFT

        # Calculate the number of unique points
        NumUniquePts = int(math.ceil((NFFT / 2) + 1))

        # FFT is symmetric, throw away second half
        fft2 = fftx[0:NumUniquePts]
        power = power[0:NumUniquePts]
        # amp = amp[0:NumUniquePts]


        # Take the magnitude of fft of x and scale the fft so that it is not a function of % the length of x
        # NOTE: If the frequency of interest is not represented exactly at one of the discrete points
        #       where the FFT is calculated, the FFT magnitude is lower.
        #
        # mx = np.abs(fftx.real) # was NumUniquePts or L but Numpy does normalization #
        # Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
        # The DC component and Nyquist component, if it exists, are unique and should not
        # be multiplied by 2.
        mx = 2 * np.abs(fft2) / NumUniquePts

        # This is an evenly spaced frequency vector with NumUniquePts points.
        # generate a freq spectrum from 0 to Fs / 2 (Nyquist freq) , NFFT / 2 + 1 points
        # The FFT is calculated for every discrete point of the frequency vector described by
        freq = np.array(list(range(0, NumUniquePts)))
        freq = freq * Fs / NFFT  # 2
        # same as
        # freq = np.fft.fftfreq(NFFT, d = dt_s)[:NumUniquePts]

        return [yd, Time, fft2, NumUniquePts, mx, freq, power]
    # end


    def WelchFourierAnalysis_overlap50pct(self, Time, SensorDepth, tunits = "sec", \
                                           window = 'hanning', \
                                           nseg = 1, log = False):
        '''
        :param nseg: number of non overlapping segments. In calculation the total number of 50% overlapping segments
                      K = N/M = 2*nseg-1
        :param Time : the time points
        :param SensorDepth: the depth data timeseries
        :param tunits: time units
        :param window: = 'blackman' #NO; 'bartlett' #OK; 'hamming' #OK; 'hann' #BEST default; flattop; gaussian;
                         blackmanharris; barthann; bartlett;

        :return: y             - detrended water levels
                 Time          - Time data points
                 fftx          - unique FFT values for the series
                 NumUniquePts  - size of the unique FFT values
                 mx            - the value of the single-sided FFT amplitude
                 f             - linear frequency vector for the mx points
                 ph            - phase vector for the mx points
        '''
        den = 2 * nseg
        N = len(Time)
        M = int(N / nseg)
        t = np.zeros((den - 1, M))  # time segments
        x = np.zeros((den - 1, M))  # data segments

        for i in range(0, den - 1):
            fct = int(N / den)
            LInt = i * fct
            RInt = LInt + M
            tt = Time[LInt:RInt]
            xx = SensorDepth[LInt:RInt]
            t[i] = tt
            x[i] = xx
        # end for

        # perform FFT
        y = np.zeros((den - 1, M), dtype = np.float)  # data segments
        Tm = np.zeros((den - 1, M), dtype = np.float)  # time segments
        fftx = np.zeros((den - 1, M / 2 + 1), dtype = np.complex)  # transform segments
        NumUniquePts = np.zeros(den - 1)  # point segments
        amplit = np.zeros((den - 1, M / 2 + 1), dtype = np.float)  # amplit segments
        f = np.zeros((den - 1, M / 2 + 1), dtype = np.float)  # freq segments
        power = np.zeros((den - 1, M / 2 + 1), dtype = np.complex)  # power segments

        for i in range(0, den - 1):
            a = self.calculateFFT(t[i], x[i], tunits, window)
            [y[i], Tm[i], fftx[i], NumUniquePts[i], amplit[i], f[i], power[i]] = a

        avg_amplit = 0
        avg_power = 0
        avg_y = 0
        avg_fftx = 0
        # calculate the average values
        for i in range(0, den - 1):
            avg_amplit += amplit[i]
            avg_power += power[i]
            # avg_y += y[i]
            avg_fftx += fftx[i]
        avg_amplit /= den - 1
        avg_power /= den - 1
        # avg_y /= den - 1 <= pointless
        avg_fftx /= den - 1

        interval_len = len(Time) / nseg
        data_len = len(Time)
        edof = fft_utils.edof(avg_amplit, data_len, interval_len, window)  # one dt chunk see Hartman Notes ATM 552 page 159 example

        (x05, x95) = fft_utils.confidence_interval(avg_amplit, edof, 0.95, log)

        return [f[0], avg_fftx, avg_amplit, avg_power, x05, x95]
    # end


    def filter(self, Time, SensorDepth, tunits):
        '''
        '''
        # prepare for the amplitude spectrum analysis
        if tunits == 'day':
            factor = 86400
        elif tunits == 'hour':
            factor = 3600
        else:
            factor = 1
        dt_s = (dates[2] - dates[1]) * factor  # Sampling period [s]
        Fs = 1 / dt_s  # Samplig freq    [Hz]

        # Filter data here on the whole lenght
        lowcut = filter[0]
        highcut = filter[1]
        btype = 'band'
        y, w, h, b, a = fft.filters.butterworth(SensorDepth, btype, lowcut, highcut, Fs, recurse = True)
        # filter.butterworth(SensorDepth, Fs)

        # ToDo:
        # here we can enable some filter display and test
        return y


    def doFFTSpectralAnalysis(self, dates, depths, tunits = 'sec', window = 'hanning', filter = None, log = False):
        '''
        FFT for Spectral Analysis
        =========================

        This example shows the use of the FFT function for spectral analysis.
        A common use of FFT's is to find the frequency components of a signal buried in a noisy
        time domain signal.

        First create some data. Consider data sampled at every 5 min for Lake Ontario and every  3 min in Frechman's bay.
        '''

        SensorDepth = np.array(depths)
        Time = np.array(dates)

        if filter:
            SensorDepth = self.filter(Time, SensorDepth, tunits)

        x05 = None
        x95 = None

        if Time[0] < 695056:
            Time += 695056
        if self.num_segments == 1:
            [y, Time, fftx, NumUniquePts, mx, f, power] = self.calculateFFT(Time, SensorDepth, tunits)  # , window, filter)
        else:
            [f, avg_fftx, avg_amplit, avg_power, x05, x95] = self.WelchFourierAnalysis_overlap50pct(Time, SensorDepth, tunits, \
                                                                                                     window, self.num_segments, log)
            fftx = avg_fftx
            mx = avg_amplit
            power = avg_power
            y = sp.signal.detrend(SensorDepth)
            NFFT = len(Time)
            NumUniquePts = int(math.ceil((NFFT / 2) + 1))

        return [y, Time, fftx, NumUniquePts, mx, f, power, x05, x95]

    # end

    def plotSingleSideAplitudeSpectrumFreq(self, data, funits = "Hz", y_label = None, title = None, log = False, \
                                           fontsize = 20, tunits = None, plottitle = False, grid = False, \
                                           ymax = None):
        '''
        :param data:
                   mx_arr = data[0]
                   name_arr = data[1]
        '''
        mx_arr = data[0]  # this is an aarray
        name_arr = data[1]
        ci05_arr = data[2]
        ci95_arr = data[3]
        f_arr = data[4]
        sSeries = []
        ci05 = []
        ci95 = []
        xa = []

        # smooth only if not segmented
        for mx in mx_arr:
            if self.num_segments == 1:
                sSeries.append(fft_utils.smoothSeries(mx, 5))
            else:
                sSeries.append(mx)

        # Plot single - sided amplitude spectrum.
        if title == None:
            title = 'Single-Sided Amplitude spectrum'

        if funits == 'Hz':
            xlabel = 'Frequency [Hz]'
            f = f_arr[0]
        elif funits == 'cph':
            xlabel = 'Frequency [cph]'
            f = f_arr[0] * 3600
        # end if

        if y_label == None:
            ylabel = 'Z(t) [m]'
        else :
            ylabel = y_label

        ya = sSeries
        if self.num_segments != 1:
            for ci in ci05_arr:
                ci05.append(ci)
            for ci in ci95_arr:
                ci95.append(ci)

        legend = []
        for i in range(0, len(data[0])):
            xa.append(f)
            legend.append(name_arr[i])

        if self.num_segments == 1:
            fft_utils.plot_n_Array(title, xlabel, ylabel, xa, ya, legend = legend, log = log, plottitle = plottitle, ymax_lim = ymax)
        else:
            fft_utils.plot_n_Array_with_CI(title, xlabel, ylabel, xa, ya, ci05, ci95, legend = legend, \
                                            log = log, fontsize = fontsize, plottitle = plottitle, grid = grid, ymax_lim = ymax)
    # end plotSingleSideAplitudeSpectrumFreq

    def plotWaveletScalogram(self, dates, data, tunits, title = None, debug = False):
        '''
        '''

        tunits = 'day'
        slevel = 0.95
        # range 0-65000 is good to catch the high frequencies
        #       0-600000 if you need to catch internal waves with much lower frequencies and large periods
        val1, val2 = (0, 130000)  # Range of sc_type (ex periods) to plot in spectogram
        avg1, avg2 = (0, 130000)  # Range of sc_type (ex periods) to plot in spectogram
        kwavelet = wavelets.kCwt.kCwt(dates, data, tunits, False)
        dj = 0.05  # Four sub-octaves per octaves
        s0 = -1  # 2 * dt                      # Starting scale, here 6 months
        J = -1  # 7 / dj                      # Seven powers of two with dj sub-octaves
        alpha = 0.5  # Lag-1 autocorrelation for white noise

        kwavelet.doSpectralAnalysis(title, "morlet", slevel, avg1, avg2, dj, s0, J, alpha)
        if debug:
            print("fftfreq=", kwavelet.wpar1.fftfreqs)
            print("amplit=", kwavelet.wpar1.amplitude)
            print("phase=", kwavelet.wpar1.phase)

        ylabel_ts = "amplitude"
        yunits_ts = 'm'
        xlabel_sc = ""
        ylabel_sc = 'Period (%s)' % kwavelet.wpar1.tunits
        # ylabel_sc = 'Freq (Hz)'
        sc_type = "period"
        # sc_type = "freq"
        # x_type = 'date'
        x_type = 'dayofyear'
        kwavelet.plotSpectrogram(kwavelet.wpar1, ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2)
    # end plotWaveletScalogram

    def plotTimeSeries(self, title, xlabel, ylabel, dates, depths, ts_legend):
        '''
        Wrapper function around fft_utils function
        '''
        fft_utils.plot_n_TimeSeries(title, xlabel, ylabel, [dates], [depths], legend = ts_legend, doy = True)
    # end

    def doDualSpectralAnalysis(self, path, filenames, names, b_wavelets = False, window = "hanning", num_segments = 1, tunits = 'day', \
                         funits = "Hz", filter = None, log = False, doy = False, grid = False, dateinterval=None, d3d=False):

        # show extended calculation of spectrum analysis
        show = True

        fftsa = FFTGraphs.FFTGraphs(path, filenames[0], filenames[1], show, tunits, )
        lake_name = names[0]
        bay_name = names[1]

        showLevels = False
        detrend = False
        draw = False
        [Time, y, x05, x95, fftx, freq, mx] = fftsa.doSpectralAnalysis(showLevels, draw, tunits, window, num_segments, filter, log, dateinterval=dateinterval, d3d=d3d)
        phase = np.zeros(len(fftx), dtype = np.float)
        deg = True
        phase = np.angle(fftx, deg)
        print("*****************************")
        print(" PHASEs")
        for i in range(0, len(fftx)):
            print("Period %f  phase:%f  amplit:%f" % (1. / freq[i] / 3600, phase[i], mx[i]))
        print("*****************************")

        plottitle = False
        ylabel = "Water level [m]"
        ylabel = "Discharge [$m^3 s^{-1}$]"
        xlabel = "Time [days]"
        
        fftsa.plotLakeLevels(lake_name, bay_name, detrend = detrend, x_label = xlabel, y_label = ylabel, plottitle = plottitle, doy = doy, grid = grid, dateinterval=dateinterval)

        fftsa.plotSingleSideAplitudeSpectrumFreq(lake_name, bay_name, funits, y_label = None, title = None, log = log, \
                                                         fontsize = 20, tunits = tunits, plottitle = plottitle, grid = grid, ymax = None)
        grid = False

        fftsa.plotPowerDensitySpectrumFreq(lake_name, bay_name, funits, plottitle = plottitle, grid = grid)
        fftsa.plotSingleSideAplitudeSpectrumTime(lake_name, bay_name, plottitle = plottitle, grid = grid)

        # fftsa.plotCospectralDensity(log = log)
        # fftsa.plotPhase()

        
    # end SpectralAnalysis

    def delta_z_resample(self, dates, depths, dt):
        dt0 = (dates[2] - dates[1])
        print("dt0=%f" % dt0)

        winlen = int(dt / dt0)
        if not isinstance(depths, np.ndarray):
            depths = np.array(depths)
            dates = np.array(dates)

        if winlen != 0 :
            y, yd = smooth.smooth(depths, winlen)
        else:
            y = depths
        y = depths
        rdepths = sp.ndimage.interpolation.zoom(y, float(dt0 / dt))
        rtime = sp.ndimage.interpolation.zoom(dates, float(dt0 / dt))
        rdzdt = [(rdepths[j + 1] - rdepths[j]) for j in range(0, len(rdepths) - 1) ]
        rdzdt.append((rdzdt[-1]+rdzdt[-2])/2)
        dzdt = [(depths[j + 1] - depths[j]) for j in range(0, len(depths) - 1) ]
        dzdt.append((dzdt[-1]+dzdt[-2])/2)
        
        return rtime, rdepths, rdzdt, dzdt

    def delta_z_mov_average(self, dates, depths, dt):
        dt0 = (dates[2] - dates[1])
        print("dt0=%f" % dt0)

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
        dzdt = [(rdepths[j + 1] - rdepths[j]) for j in range(0, len(rdepths) - 1) ]
        dzdt.append(dzdt[-1])
        dzdta = np.array(dzdt)
        return rtime, dzdta


    def plot_dzdt_up_line(self, up, dzdt, labels):
        display_data.display_scatter(up, dzdt, slope = True, type = 'stats' , labels = labels)


    def plot_cross_spectogram_u_dz(self, t, w, tt, dz, scaleunit = 'day', da = [6, 300]):
        '''
        Plot the  cross wavelet spectrogram between the ADCP velocity and Temperature
        :param t : velocity time array
        :param w : longitudinal velocity
        :param tt : wl time array
        :param dz : dz/dt - water level variation
        '''
        print("Start cross spectrum analysis")


        # w = u + 1j * v
        if len(w) != len(dz):
            ratio = float(len(w)) / len(tt)
            print("len w: %d | len t: %d | len tt: %d | len T: %d" % (len(w), len (t), len(tt), len(dz)))
            print("ratio: %f" % ratio)
            #oldrange = range(0, ratio * len(t), ratio)
            #print "len(oldrange) :%d" % len(oldrange)
            if ratio < 1 :
                idz = np.interp(t, tt, dz)
                # iT2 = numpy.interp(range(0, ratio * len(t)), oldrange), T)
                print("len w: %d | len iT: %d |  len T: %d" % (len(w), len(idz), len(dz)))
                tt = t
                iw = w
            else:
                iw= np.interp(tt, t, w)
                t = tt
                idz = dz
        else:
            idz = dz
            iw  = w 

        kwavelet = wavelets.kCwt.kCwt(t, iw, tt, idz, scaleunit, False, True)
        dj = 0.05  # Four sub-octaves per octaves
        s0 = -1  # 2 * dt                              # Starting scale, here 6 months
        J = -1  # 7 / dj                               # Seven powers of two with dj sub-octaves
        alpha = 0.5  # Lag-1 autocorrelation for white noise
        slevel = 0.95
        val1, val2 = (0, 520000)  # Range of sc_type (ex periods) to plot in spectogram
        avg1, avg2 = (0, 520000)  # Range of sc_type (ex periods) to plot in spectogram
        title = "Cross vel-dz/dt spectrum"
        print(title + "start: Monte Carlo")
        kwavelet.doCrossSpectralAnalysis(title, "morlet", slevel, avg1, avg2, dj, s0, J, alpha)
        print(title + "end: Monte Carlo")

        ylabel_ts = "PSD"
        yunits_ts = 'm/s'
        xlabel_sc = ""
        ylabel_sc = 'Period [%s]' % kwavelet.wpar1.tunits
        sc_type = "period"
        # sc_type = "freq"
        # x_type = 'date'
        x_type = 'dayofyear'


        print("#plot cross wavelet spectrogram")
        kwavelet.plotXSpectrogram(kwavelet.get_xwt(), extend = 'both', x_type = x_type, ylabel_sc = ylabel_sc, da = da,
                      tfactor = kwavelet.wpar1.tfactor, scale = 'linear')

        print("plot coherence wavelet spectrogram")
        kwavelet.plotXSpectrogram(kwavelet.get_wct(), extend = 'neither', x_type = x_type, ylabel_sc = ylabel_sc, da = da,
                      tfactor = kwavelet.wpar1.tfactor, crange = np.arange(0, 1.1, 0.1), scale = 'linear', angle = kwavelet.get_wct().angle)



    def interpolate_to_delft3d(self, fn, step_min, start_WL, timesince2001):
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
        return [rdepths, rtime]

    def convert_wl_to_delft3d_tim(self, path, fn, step_min, start_WL, timesince2001):
        dtime = datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
        st_strg = dtime.strftime("%y%m%d")
        [rdepths, rtime] = interpolate_to_delft3d(fn, step_min, start_WL, timesince2001)     
        #write to 
        outfn=path+"/"+ st_strg + fn+".tim"
        ofile = open(outfn, "wt")
        writer = csv.writer(ofile, delimiter = '\t', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
        #convert time to min starting from 0
        start = rtime[0]
        deltawl = start_WL - rdepths[0]
        for i in range(len(rtime)-1):
            rtime[i] = int((rtime[i]-start)*24*60)           #convett to min
            tm = "%d" % rtime[i]
            rdp = "%.3f" % (rdepths[i] + deltawl)
            writer.writerow([tm, timesince2001[i], rdp])
        ofile.close    
       
    def calculate_min_since_20010101000000(self, step_min):
        d1 = datetime.strptime(self.date[0], "%y/%m/%d %H:%M:%S")
        d2001 = datetime(2001,0o1,0o1,00,00,00)
        d2 = datetime.strptime(self.date[1], "%y/%m/%d %H:%M:%S")
        dt= d1-d2001
        totmin1 = int(dt.total_seconds()/60)
        dt= d2-d2001
        totmin2 = int(dt.total_seconds()/60)
        timearr = []
        for i in range(totmin1,totmin2, step_min):
            timearr.append(i)
        return timearr
        
       
if __name__ == '__main__':
    path = '/home/bogdan/Documents/UofT/PhD/Data_Files/2013/Hobo-Apr-Nov-2013/WL/csv_press_corr'
    filenames = {'Emb A':'10279443_corr.csv',
                 'Emb B':'1115681_corr.csv',
                 'Emb C':'10238147_corr.csv',
                 'Cell 1':'10279696_corr.csv',
                 'Cell 2':'10279693_corr.csv',
                 'Cell 3':'10279699_corr.csv',
                 'Out Harb':'10279444_corr.csv'}

    num_segments = 10
    wla = WaterLevelAnalysis(path, filenames, num_segments)

    for key in sorted(wla.getDict().keys()):
        fname = wla.getDict()[key]

        [dates, depths] = wla.read_press_corr_file(path, fname)

        # plot the original Lake oscillation input
        xlabel = 'Time [days]'
        ylabel = 'Z(t) [m]'
        ts_legend = [key + ' - water levels [m]']
        fft_utils.plotTimeSeries("Lake levels", xlabel, ylabel, dates, depths, ts_legend)
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

        wla.plotWaveletScalogram(dates, depths, tunits)

    print("Done !!")
