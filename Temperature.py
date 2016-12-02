'''
Created on Nov 20, 2014

@author: bogdan
'''
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import math
import csv
from scipy import fftpack
from datetime import datetime
import matplotlib.dates as dates
# local imports
import ufft.fft_utils as fft_utils
import wavelets.kCwt
import ufft.FFTGraphs as FFTGraphs
import utools.display_data as display_data
import utools
import utools.stats as stat
import ufft.smooth as smooth
import matplotlib

class Temperature(object):
    '''
    classdocs
    '''


    def __init__(self, path, filenames, numsegments, tinterv = None):
        self.path = path
        self.dict = filenames
        self.num_segments = numsegments
        if tinterv != None:
            self.date = [tinterv[0], tinterv[1]]
        else:
            self.date = None
            
            
       
        
    @staticmethod
    def DTemp_Model_vs_Meas(Tres_lim, Tres_step, D_lim, D_step, H_lim, H_step, xlabel, ylabel, cblabel, \
                            points, lake, range_x=None, range_y=None, \
                            fontsize = 20, chains = None, exclusions = None):
        
        
        Cp = 4185.5   #J/(kg*K) (15 C, 101.325 kPa)
        rho = 999.1026 #kg/m^3
        tres = np.linspace(Tres_lim[0], Tres_lim[1], Tres_step)
        H = np.linspace(H_lim[0], H_lim[1], H_step)
        D = np.linspace(D_lim[0], D_lim[1], D_step)
        
        X, Y = np.meshgrid(tres,D)
        
        secs_in_day=24*3600
        
        fig = plt.figure(facecolor = 'w', edgecolor = 'k')
        format = 100 * H_step + 10
        ax = np.zeros(H_step, dtype = matplotlib.axes.Subplot)
        for h in range(0, H_step):
            dT= np.zeros((D_step, Tres_step))
            if (h == 0) :
                ax[h] = fig.add_subplot(format + h + 1)
            else:
                ax[h] = fig.add_subplot(format + h + 1, sharex = ax[0])
            
        
            for j in range(0, len(D)):
                for i in range(0, len(tres)):
                    dT[j][i] = H[h]/(D[j]*rho*Cp)*(tres[i]*secs_in_day) 
    
            #im = ax[h].pcolormesh(X, Y, dT, shading = 'gouraud', vmin = -30, vmax = 0)
            im = ax[h].pcolormesh(X, Y, dT, shading = 'gouraud', vmin = 0, vmax = 30)
            
            cb = fig.colorbar(im, ax = ax[h], aspect = 18)

            #cb.set_clim(mintemp, maxtemp)
            labels = cb.ax.get_yticklabels()
            for t in labels:
                t.set_fontsize(fontsize - 4)
           
            from matplotlib import ticker
            tick_locator = ticker.MaxNLocator(nbins = 5)
            cb.locator = tick_locator
            cb.update_ticks()
            cb.set_label(cblabel)
            text = cb.ax.yaxis.label
            font = matplotlib.font_manager.FontProperties(size = fontsize-3)
            text.set_font_properties(font)
            
            if h == H_step - 1:
                ax[h].set_xlabel(xlabel).set_fontsize(fontsize)
                plt.setp(ax[h].get_xticklabels(), visible = True)
            else:
                plt.setp(ax[h].get_xticklabels(), visible = False)
            
            ax[h].set_ylabel(ylabel).set_fontsize(fontsize)
            
            for tick in ax[h].xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize - 3)
            for tick in ax[h].yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize - 3)
        
            ax[h].set_yticks(np.arange(min(D), max(D)+1, 2.0))
        
            #set the points and the text
            xmin = 100; ymin = 100
            for k,v in points.items():
                xmin = min(v[0][0], xmin)
                ymin = min(v[0][1], ymin)
             
            ii=0
            for k,v in points.items():
                text = k
                if text == "LO":
                    continue
                location = v[0]
                if text == 'Ah' or text == 'Al':
                    ax[h].plot(location[0], location[1], 'dr')
                else:
                    ax[h].plot(location[0], location[1], 'dw')
                
                ax[h].annotate(text, xy=(location[0], location[1]), xytext = (location[0] + xmin/5., location[1] + ymin/15.),\
                               color='white', fontsize = fontsize-5)
                ii+=1
    
            #testing images
            #===================================================================
            #  
            # D2 = im.get_array().reshape(len(X), len(Y))
            # #fig.close()
            # #assert  np.array_equal(D, D2)
            # # get the values from the
            # calcDTemp = []
            # measTemp = []
            # point_labels=[]
            # print ("------------------------------------")
            # for k, v in points.items():
            #     # WHY IS REVERSED HERE?
            #     name = k 
            #     if name == "LO":
            #         continue
            #     point_labels.append(name)
            #     xp= int((v[1])*len(Y)/range_y)
            #     yp=int((v[0])*len(X)/range_x)
            #     calcDTemp.append(D2[xp,yp])
            #     measTemp.append(v[2])  
            #     print ("bay:%s dtMeas:%f dtCalc:%f " % (name, v[2], D2[xp,yp]))
            # print ("------------------------------------")
            # 
            # #create the regression
            # [r2, slope, intercept, r_value, p_value, std_err] = stat.rsquared(measTemp, calcDTemp)
            # stat.plot_regression(np.array(measTemp), np.array(calcDTemp), slope = slope, intercept = intercept, point_labels = point_labels, \
            #                      x_label = "Measured $\Delta$T [$^\circ$C]", y_label = "Calculated $\Delta$T [$^\circ$C]", \
            #                      r_value = r_value, p_value = p_value, show=True, \
            #                      title = "Measured vs. Calculated $\Delta$T")
            #===================================================================
            
            #use the chains
            
            if chains is not None:
                deltaT = []
                point_labels=[]
                measTemp = []
                stddev = []
                
                loTemp=lake["LO"][3]
                if exclusions != None:
                    deltaT_excl = []
                    measTemp_excl = []
                    stddev = []
                for name, data in points.items():
                    if name == "LO" or name == "EB" or name == "Al" or name == "Ah": # if exclude EB as well
                        continue
                    D = data[0][1]
                    resTime = data[0][0]
                    meanT = data[0][3]
                    Heat=H[h]
                    cT=Heat/(D*rho*Cp)*(resTime*secs_in_day)
                    deltaT.append(cT)
                    measDT=meanT-loTemp
                    measTemp.append(measDT)
                    stddev.append(data[1])
                    if exclusions != None and name not in exclusions:
                        deltaT_excl.append(cT) 
                        measTemp_excl.append(measDT)
                        
                    print("name:%s, measDT:%f  calcDT:%f" %(name, measDT,cT))
                    point_labels.append(name)
                if exclusions != None:
                    [r2, slope, intercept, r_value, p_value, std_err] = stat.rsquared(measTemp_excl, deltaT_excl)
                else:
                    [r2, slope, intercept, r_value, p_value, std_err] = stat.rsquared(measTemp, deltaT)
                stat.plot_regression(np.array(measTemp), np.array(deltaT), stddev, slope = slope, intercept = intercept,
                                     point_labels = point_labels,
                                     x_label = "Measured $\Delta$T [$^\circ$C]",
                                     y_label = "Calculated $\Delta$T [$^\circ$C]",
                                     r_value = r_value, p_value = p_value, show=True, exclusions = exclusions) #,
                                 #title = "Measured vs. $\Delta$T Calculated Directly")
            
    def getDict(self):
        return self.dict

    @staticmethod
    def read_ctd_file(path, fname):
        ifile = open(path + '/' + fname, 'rb')
        reader = csv.reader(ifile, delimiter = ',', quotechar = '"')
        temps = []
        depths = []

        for row in reader:
            try:
                temp = float(row[2])
                temps.append(temp)
                depths.append(float(row[3]))
            except:
                print("Error:read_temp_file")
            # end try

        ifile.close()

        return [temps, depths]

    @staticmethod
    def plot_ctds(path, files):
        flist = ["CTD01.csv","CTD02.csv","CTD03.csv","CTD04.csv","CTD05.csv","CTD06.csv","CTD07.csv","CTD15.csv"]
        dlist = []
        tlist = []
        for f in files:
            [temps, depths] = Temperature.read_ctd_file(path, f)
            
            #converf from decibar totakl pressure in water depth
            depths = [d -10 for d in depths]
            display_data.display_simple_temperature_profiles([depths], [temps], f, revert = False, legendloc = 4)
            if f in flist:
                dlist.append(depths[::-1])
                tlist.append(temps)
        # 3) Mixed water, air ,img data
        custom = np.array(["Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]", "Depth [m]" ])
        # ToDO: Add short and long radiation
        print("ctd Start display mixed subplots  ")
         
        data1 = dlist
        dateTimes1 = tlist
        ylabels = custom
        limits1 = [[0,4],[0,4],[0,4],[0,4],[0,4],[0,4],[0,4],[0,4]]
        utools.display_data.display_mixed_subplot(dateTimes1 = dateTimes1, data = data1, varnames = flist, ylabels1 = ylabels, limits1 = limits1,\
                                           dateTimes2 = [], groups = [], groupnames = [], ylabels2 = [], \
                                           dateTimes3 = [], imgs = [], ylabels3 = [], ticks = [], maxdepths = None, \
                                           mindepths = None, mintemps = None, firstlogs = None, maxtemps = None, \
                              fnames = None, revert = True, custom = None, maxdepth = None, tick = None, firstlog = None, yday = True, \
                              title = False, grid = False, limits = None, sharex = True, fontsize = 20, group_first = False, interp = None)

            
    def read_temp_file(self, path, fname):
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
                time = float(row[1])

                if self.date != None:
                    if time < startt or time > endt:
                        continue
                dateTime.append(time)
                depths.append(float(row[2]))
            except:
                print("Error:read_temp_file")
            # end try

        ifile.close()

        return [dateTime, depths]

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
            ylabel = 'Temperature [$^\circ$C]'
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
            LInt = (i * fct)
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
    
    def stddev(self):
        out = []
        for fn in self.dict:
            [date,temps] = self.read_temp_file(self.path, fn)
            stddev = np.std(temps)
            mean = np.average(temps)
            MinSD = stddev
            MaxSD = stddev
            out.append([MinSD,MaxSD])
        return out
        
