import numpy
import scipy


def resample(time, data, dt):
    """
    :param time: time array  time[1]-tim[0] gives the current interval
    :param data: times series data to be resampled
    :param dt:   new time interval
    :return: resampled time, resampled data
    """
    dt0 = (time[1] - time[0])
    print("dt0=%f" % dt0)
    winlen = float(dt0 / dt)
    y = data
    rdata = scipy.ndimage.interpolation.zoom(numpy.real(y), winlen)
    rtime = scipy.ndimage.interpolation.zoom(time, winlen)

    return rtime, rdata
