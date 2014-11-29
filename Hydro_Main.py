'''
Created on Nov 20, 2014

@author: bogdan
'''

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
       }

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
        xlabel = 'Time (days)'
        ylabel = 'Z(t) (m)'
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

        y_label = '|Z(t)| (m)'
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
                 'Out Harb':'10279444_corr.csv'}

    for key, value in filenames.iteritems():
        fnames = [filenames['Out Harb'], value]
        names = ['Out Harb', key]
        wla = Water_Level.WaterLevelAnalysis(paths[1], fnames, num_segments)
        wla.doDualSpectralAnalysis(paths[1], fnames, names, b_wavelets = False, window = "hanning", \
                                   num_segments = num_segments, tunits = 'day', \
                                   funits = "cph", filter = None, log = False, doy = True, grid = False)

def Vel_FFT_pairs(date):

    # Process RDI-Teledyne
#===============================================================================
#     num_segments = 10
#     filenames = {'OH':'600mhz-DPL_002.000',
#                  'EmbC':'1200mhz-EMBC_004.000'}
#
#     for key, value in filenames.iteritems():
#         hyd = Hydrodynamic.Hydrodynamic(value, paths[2], num_segments, date)
#         hyd.readRawBinADCP()
#         bins = [0, 4]
#         # type can be only ampl here
#         type = 'ampl'
#         # type = 'power'
#         hyd.plot_FFT(key, bins, tunits = "day", funits = "cph", log = False, grid = False, type = type, withci = True)
#===============================================================================

    # process PC ADP
    ctlfile = 'MODE006.ctl'
    fnames = ['MODE006.ve', 'MODE006.vn', 'MODE006.vu']
    num_segments = 10
    name = 'Cell 3'
    hydpcadp = Hydrodynamic.Hydrodynamic(name, paths[3], num_segments, date, ctlname = ctlfile, filenames = fnames)
    hydpcadp.readPcAdpVel()
    bins = [0, 3]
    hydpcadp.plot_FFT(name, bins, tunits = "day", funits = "cph", log = False, grid = False, type = type, withci = True, sel_dates = False)




if __name__ == '__main__':

    date = ['13/06/25 00:00:00', '13/10/27 00:00:00']
    # zoom in
    # date = ['13/07/18 00:00:00', '13/07/27 00:00:00']

    # date = None

    v = 'vel_fft_pairs'
    # map the inputs to the function blocks
    for case in switch(v):
        if case('wl_fft_all'):
            WL_FFT_analysis_all()
            break
        if case('wl_fft_pairs'):
            WL_FFT_pairs()
            break
        if  case('vel_fft_pairs'):
            Vel_FFT_pairs(date)
            break
        if case():  # default, could also just omit condition or 'if True'
            print ("something else!")
            # No need to break here, it'll stop anyway

    print "Done!"
