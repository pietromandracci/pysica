# COPYRIGHT 2020 by Pietro Mandracci

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

""" PYthon tools for SImulation and CAlculus: mathematical functions.

    This module contains some mathematical functions used to compute averages, derivatives and so on.

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import numpy
from scipy.stats import norm


# +------------------------+
# | Mathematical functions |
# +------------------------+

def interpolate(x, x0, x1, y0, y1):
    """ Returns the linear interpolation of the given values

        Parameters
        ----------

        x0: lower value of independent variable
        x1: upper value of indepenfent variable
        y0: lower value of dependent variable
        y1: upper value of dependent variable

        Returns
        -------

        y = y0 + (y1-y0) * (x-x0) / (x1-x0)
    """

    return y0 + (y1-y0) * (x-x0) / (x1-x0)


def interpolate_3p(x, x0, x1, x2, y0, y1, y2):
    """ Returns the interpolation of the given values using 3 points

        Parameters
        ----------

        x0: lower  value of independent variable
        x1: medium value of independent variable
        x2: upper  value of independent variable
        y0: lower  value of dependent variable
        y1: medium value of dependent variable
        y2: upper  value of dependent variable

        Returns
        -------

        y = interpolated value, using 3 points
    """

    d0  = x  - x0
    d1  = x  - x1
    d2  = x  - x2
    d01 = x0 - x1
    d02 = x0 - x2
    d12 = x1 - x2
    d10 = -d01
    d20 = -d02
    d21 = -d12
    L0  =  (d1 * d2) / (d01 * d02)
    L1  =  (d0 * d2) / (d10 * d12)
    L2  =  (d0 * d1) / (d20 * d21)
    y   = L0 * y0 + L1 * y1 + L2 * y2

    return y


def _derivative_5p(dx, y0minus2, y0minus1, y0, y0plus1, y0plus2):
    """ Calculates the 5-points derivative of function f(x) at point x0.

        Parameters
        ----------

        dx         finite difference of x-value
        y0minus2   f(x0-2*dx)
        y0minus1   f(x0-dx) 
        y0         f(x0)
        y0plus1    f(x0+dx)
        y0plus2    f(x0+2*dx)

        Returns
        -------

        y1 = f'(x0)
    """

    return (y0minus2 - 8 * y0minus1 + 8 * y0plus1 - y0plus2 ) / (12 * dx)


def derivative_5p(dx, y):
    """ Calculates the derivative of a function y=f(x), the y values are given as a numpy array.

        Parameters
        ----------
        dx: distance between the x-points for which the y values were given
        y: numpy array contaning

        Returns
        -------

        y1: a numpy array containing the values of f'(x) 
            the first 2 and last 2 values are nan
    """

    n = len(y)
    y1 = numpy.ones(n) * numpy.nan
    for i in range(2,n-2):
        y1[i] = _derivative_5p(dx, y[i-2], y[i-1], y[i], y[i+1], y[i+2])
                
    return y1


def median_filter(a, window, pad=False, spikes=False, k=1):
    """ Returns the moving average of an array, using a window of 2w+1 points

        Parameters
        ----------

        a:       the array to be filtered

        window:  width of the mediana window, in points
                 the i-th point of the averaged array will be the median
                 of points from i-window to i+window (included) of a
                 this value must be lower than len(a)/2
        pad:     if False, the first and last w points of the new spectrum will be padded with nan
                 if True, the window used to calculate median is padded with the first/last point
        spikes:  if True a procedure to remove spike is applied,
                 otherwise, the standard median filter is applied 
        k:       used only if spikes is True: points will be cosidered spikes if their y value 
                 differs (in modulus) from the median of the 2w+1 windows 
                 of more than k times the maximum 


        Returns
        -------

        a_median:  the medinan filtered array, since the filter uses a number of points equal 
                   to 'window' before and after each one for which it is calculated,
                   it is not possibile to calculate it for a number of values equal to 'window'
                   at the beginning and end of the array: these will be padded with nan 
    """

    n = len(a)

    if (window < 0): window = - window
    if (window >= n / 2): window = n / 2 - 1
    window = int(window)
    if (k < 0): k = - k
    a_median = numpy.ones(n) * numpy.nan

    # This array will be used in the unspike mode to remove the central point of the window
    # in order to calculate the maximum deviation from the median
    if spikes:
        window_indexes         = numpy.ones(2*window+1,'bool')
        window_indexes[window] = False        

    for i in range(n):
        # Define the window of points from which the median
        # will be calculated, from i-window to i+window, both included
        # if the point is too near the the edges, and pad is True,
        # the first or last point is used to fill the window
        if (i < window):
            if pad:
                a_window = numpy.zeros(2*window+1,'d')
                for j in range(window-i):            a_window[j] = a[0]
                for j in range(window-i,2*window+1): a_window[j] = a[j-window+i]
            else:
                continue
        elif (i >= n - window):
            if pad:
                a_window = numpy.zeros(2*window+1,'d')
                for j in range(window+n-i):            a_window[j] = a[j-window+i]
                for j in range(window+n-i,2*window+1): a_window[j] = a[n-1]
            else:
                continue
        else:
            a_window = a[i-window:i+window+1]
            
        median = a_window[a_window.argsort()][window]
                
        if spikes:
            # Calculate the maximum deviation from the median 
            max_dev = abs(a_window[window_indexes] - median ).max()
            # Substitute the point with the median,
            # if it differs from the median more than k times
            # the maximum deviation from the median of the other points
            if (abs(a[i]-median) > k * max_dev):
                a_median[i] = median
            else:
                a_median[i] = a[i]
        else:
            a_median[i] = median

    return a_median


def moving_average(a, window):
    """ Returns the moving average of an array, using a window of 2w+1 points

        Parameters
        ----------

        a: the array to be averaged

        window: width of the average window, in points
                the i-th point of the averaged array will be the average 
                of points from i-window to i+window (included) of a
                this value must be lower than len(a)/2

        Returns
        -------

        a_averaged: the averaged array, since the average uses a number of points equal 
                    to 'window' before and after each one for which it is calculated,
                    it is not possibile to calculate it for a number of values equal to 'window'
                    at the beginning and end of the array: these will be padded with nan 
    """

    n = len(a)

    if (window < 0): window = - window
    if (window >= n / 2): window = n / 2 - 1
    window = int(window)

    a_smooth = numpy.ones(n) * numpy.nan

    # First point (i = window)
    mov_average = 0
    for j in range(0, 2 * window + 1):
        mov_average = mov_average + a[j]
    mov_average = mov_average / (2 * window + 1)
    a_smooth[window] = mov_average

    # Following points (i = window + 1 .. n - window - 1)
    for i in range(window + 1, n - window):
        add_value = (a[i+window] - a[i-window-1]) / (2 * window + 1)
        mov_average = mov_average + add_value
        a_smooth[i] = mov_average

    return a_smooth


def moving_average_skip(a, window, skip=0):
    """ Returns the moving average of an array, using a window of 2w+1 points, 
        allowing to skip some points.

        Parameters
        ----------

        a: the array to be averaged

        window: width of the average window, in points
                the i-th point of the averaged array will be the average 
                of points from i-window to i+window (included) of a
                this value must be lower than len(a)/2
        skip:   number of values to skip at the beginning and end of the array

        Returns
        -------

        a_averaged: the averaged array, since the average uses a number of points equal 
                    to 'window' before and after each one for which it is calculated,
                    it is not possibile to calculate it for a number of values equal to 'window'
                    at the beginning and end of the array: these will be padded with nan 
    """

    n = len(a)

    if (skip < 0):      skip = - skip
    if (skip >= n / 2): skip = n / 2 - 1
    skip = int(skip)

    m = n - 2 * skip
    if (window < 0):      window = - window
    if (window >= m / 2): window = m / 2 - 1
    window = int(window)

    a_smooth = numpy.ones(n) * numpy.nan

    # First point (i = window + skip)
    mov_average = 0
    for j in range(skip, 2 * window + 1 + skip):
        mov_average = mov_average + a[j]
    mov_average = mov_average / (2 * window + 1)
    a_smooth[window+skip] = mov_average

    # Following points (i = window + 1 + skip .. n - window - 1 - skip)
    for i in range(window + 1 + skip, n - window - skip):
        add_value = (a[i+window] - a[i-window-1]) / (2 * window + 1)
        mov_average = mov_average + add_value
        a_smooth[i] = mov_average

    return a_smooth


def calculate_weights(n, method='linear', sigma=None, k=3):
    """ Returns an arrays of weights to be used for the calculation of a weighted moving average.

        Parameters
        ----------

        n:      number of weight
        method: method used to calculate weights: 'uniform', 'linear' or 'gaussian'
                'uniform'   -> uniform weights
                'linear'    -> linear weights
                'gaussian'  -> gaussian weights, with stadard deviation equal to (n-1)/k
        sigma:  used in the 'gaussian' mode only
                it is the standard deviation of the normal distribution used for weights calculation
                if not given, it will be calculated
        k:      used in the 'gaussian' mode, if the sigma parameter is not given
                the standard deviation of the normal distribution used 
                for wieghts calculation will be 
                    sigma = (n-1)/k
                    default is k=3, so that the the sum of the weights is ~> 0.99

        Returns
        -------

        weights: an array containing the weights 

    """

    n = int(n)
    weights = numpy.ones(n,'d') * numpy.nan

    w = n - 1
    if (method == 'uniform'):                
        for i in range(n):
            weights[i] = 1.0 / (2*n - 1)               
    elif (method == 'linear'):
        for i in range(n):
            weights[i] = 1.0 * (w - i + 1) / ( (w + 1) * (w + 1) )
    elif (method =='gaussian'):
        # If sigma was not given, we calculate it in such a way that k*sigma = w+0.5
        # since each weight is calculated as the gaussian probability from (i-0.5) to (i+0.5)
        if (sigma is None): sigma = 1.0 * (w+0.5) / k 
        for i in range(n):
            inf = i - 0.5
            sup = i + 0.5
            P_inf = norm.cdf(inf, loc=0, scale=sigma)
            P_sup = norm.cdf(sup, loc=0, scale=sigma)
            weights[i] = P_sup - P_inf
    return weights
        
def moving_average_weighted(a, w, normalize=False, skip=0):
    """ Returns the weighted moving average of an array, using an array of weights

        Each point of the averaged array will be given by
        average[i] = sum( j, from i-l to i+l of array[j] * weights[abs(j-i)] )
               
        The weights w[i] must satisfy these properties
        w[j] = w[-j]
        sum(i, from -l to l, weights[i]) = 1

        The array given contains the values weights[0] to weights[l]


        Parameters
        ----------

        a:         the array to be averaged
        w:         the array of weights, these must be real numers belonging to the interval [0,1]
        normalize: if True, each weight is divided by the sum of the weights, so that their 
                   sum becomes 1
        skip:      numebr of points to skip at the beginning and end of the signal

        Returns
        -------

        a_averaged: the averaged array, since the average uses a number of points equal 
                    to len(w)-1 before and after each one for which it is calculated,
                    it is not possibile to calculate it for a number of values equal to 'window'
                    at the beginning and end of the array: these will be padded with nan 
    """

    n =      len(a)
    window = len(w) - 1

    # If required, normalize the weights.
    if (normalize):
        N = 2 * w.sum() - w[0]
        w = w / N

    a_smooth = numpy.ones(n,'d') * numpy.nan

    # If the array length is not enough to calculate the average
    # for at least one point, return an array of NaN
    if (n < 2 * window + 1): return a_smooth

    # From a[window] to a[n-window-1]
    for i in range(window + skip, n - window - skip):
        mov_average = 0.0
        # weighted average of a[i-window] to a[i+window]
        for j in range(i-window, i+window+1):
            mov_average = mov_average + a[j] * w[abs(j-i)]
        a_smooth[i] = mov_average

    return a_smooth


def roughness(z):
    """ Returns the roughness of a surface profile.

        Parameters
        ----------

        z: surface profile

        Returns
        -------

        z_mean, Ra, Rq, Rsk: mean height, Ra value, Rq value, Rsk (skewness), Rku (kurtosis)
        in case of error, None is returned
    """

    # Remove NaNs and infs from the signal
    z_real      = z[numpy.isfinite(z)]
    n_real      = len(z_real)
    z_mean      = z_real.mean()
    z_delta     = z_real - z_mean
    z_delta2    = z_delta  * z_delta
    z_delta3    = z_delta2 * z_delta
    z_delta4    = z_delta3 * z_delta
    Ra          = numpy.abs(z_delta).sum() / n_real
    Rq          = numpy.sqrt( z_delta2.sum() / n_real )
    Rsk         = z_delta3.sum() / (n_real * Rq * Rq * Rq)
    Rku         = z_delta4.sum() / (n_real * Rq * Rq * Rq * Rq)         

    return (z_mean, Ra, Rq, Rsk, Rku)
