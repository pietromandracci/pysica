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

""" Tools for statistical analysis on univariate experimental data.

This module contains some classes, methods and functions, which 
may be useful to perform some basic statistical analysis on a 
univariate set of experimental data, such as some kind of data 
plots and comparison of datasets to theoretical distributions.        

Documentation is also available in the docstrings.
"""

# +-------------------------+
# | Import required modules |
# +-------------------------+

# Mudules from the standard Python library
import math

# Modules from the Python community
import numpy
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import chi2

# Mudules from pysica package
from ..managers import data_manager
from ..managers import gnuplot_manager as gm
from ..parameters import *
from ..functions.pdf import *
from ..functions.statistics import *



# +------------+ 
# | Histograms |
# +------------+

class DataSet:
    """ Class defining a univariate set of experimental data. """

    def __init__(self, data_array=None, filename=""):
        """ Defines an instance of the class DataSet and calculates 
            some basic statistical parameters.

            Parameters
            ----------

            data_array:  a numpy 1d array, if not given, will be read 
                         from a file
            filename:    name of the file from which the data 
                         must be read. 
                         If data_array is None, this must not be None

            Initialized data attributes
            ---------------------------

            self.data:           Raw data
            self.n_data:         Number of data values
            self.max:            Max data value
            self.min:            Min data value
            self.width:          Width of data interval (max-min)
            self.median:         Median of data values
            self.percentiles:    Array of percentile values (5% step)
            self.mean:           Average of data
            self.residuals:      Data residuals (difference of each datum 
                                 from the mean)
            self.residuals_2:    Squares of residuals
            self.residuals_3:    Third power of residuals
            self.var:            Variance of data
            self.stdv:           Standard deviation of data
            self.mu3:            Third moment of data 
            self.mu4:            Fourth moment of data 
            self.skewness:       Skewness of data
            self.kurtosis:       Excess curtosis of data
            self.var_bias:       Variance of data (with bias)
            self.stdv_bias:      Standard deviation of data (with bias)
            self.mu3_bias:       Third moment of data (with bias) 
            self.mu4_bias:       Fourth moment of data (with bias) 
            self.skewness_bias:  Skewness of data (with bias)
            self.kurtosis_bias:  Excess curtosis of data (with bias)
            self.n_intervals:    Number of intervals (initialized to zero)
            self.observed_ready: True when observed freq have been calculated
            self.expected_ready: True when expected freq have been calculated
            self.read_error:     used for error checking
                                 (status, message)
                                 status: 0 = no error
                                         1 = file not found
                                 message: an error message or 'OK'
        """

        self.read_error = (0,OK)
        
        # Raw experimental data
        if (data_array is None):
            d = data_manager.DataSequence()
            self.read_error = d.read_file(filename)
            if (self.read_error[0] != 0):
                del d
                return
            self.data = d.data_array
            del d
        else:
            self.data = data_array
                        
        # Number of experimental data values
        self.n_data = len(self.data)
        
        # Parameter used to calculate unbiased variance
        if (self.n_data > 1):
            self._var_bias_corr = 1.0 * self.n_data / (self.n_data-1)     
        else:                                                                   
            self._var_bias_corr = 0.0
            
        # Parameter used to calculate unbiased skewness     
        if (self.n_data > 2):        
            self.sk_bias_corr = (math.sqrt(1.0*self.n_data*(self.n_data-1))
                                 / (self.n_data-2))
        else:
            self.sk_bias_corr = 0.0                    

        # Parameter used to calculate unbiased kurtosis
        if (self.n_data > 3):
            self.ku_bias_corr = ((1.0 * self.n_data - 1)
                                 / ((self.n_data - 2)
                                    * (self.n_data - 3)))
        else:
            self.ku_bias_corr = 0.0
                        
        # Maximum of data values
        self.max = self.data.max()

        # Minimum of data values
        self.min = self.data.min()

        # Width of data interval
        self.width = self.max - self.min

        # Median and percentiles
        # numpy.linspace(5,100,20) = array of percentile values (5% step)
        self.percentiles = numpy.percentile(self.data,
                                            numpy.linspace(5,100,20))
        # Median of data values        
        self.median  = self.percentiles[9]
#       self.median  = numpy.median(self.data)

        # Mean
        self.mean = numpy.mean(self.data)
                
        # Residuals and their powers
        self.residuals   = self.data - self.mean
        self.residuals_2 = self.residuals * self.residuals
        self.residuals_3 = self.residuals_2 * self.residuals
        self.residuals_4 = self.residuals_3 * self.residuals
                
        # Biased parameters
       
        # Variance of data, with bias
        self.var_bias = self.residuals_2.sum() / self.n_data
        # Third moment of data, with bias
        self.mu3_bias = self.residuals_3.sum() / self.n_data
        # Fourth moment of data, with bias       
        self.mu4_bias = self.residuals_4.sum() / self.n_data
        # Standard deviation of data, with bias
        self.stdv_bias = math.sqrt(self.var_bias)
        # Skewness of data, with bias
        self.skewness_bias = self.mu3_bias / self.stdv_bias**3
        # Excess kurtosis of data, biased
        self.kurtosis_bias = self.mu4_bias / self.var_bias**2 - 3          

        # Bias-corrected parameters

        # Variance of data, bias corrected
        self.var = self._var_bias_corr * self.var_bias
        # Standard deviation of data, bias corrected
        self.stdv = math.sqrt(self.var)
        # Skewness of data, bias corrected
        self.skewness = self.skewness_bias * self.sk_bias_corr
        # Kurtosis, bias corrected
        self.kurtosis = (self.ku_bias_corr
                         * ((self.n_data + 1) * self.kurtosis_bias + 6))

        # Number of intervals (initialized to None)                
        self.n_intervals = None

        # Flags activated when observed/expected frequencies
        # have been calculated
        self.observed_ready = False
        self.expected_ready = False       
        

    def print_stats(self):
        """ Prints some basic statistical information about data.

            Returns
            -------

            (status, message)
            status:  0 = no error
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        print("Number of values:          ", self.n_data)
        print("Minimum value:             ", self.min)
        print("Maximum value:             ", self.max)
        print("Data range width:          ", self.width)
        print("Median:                    ", self.median)                
        print("5%  percentile:            ", self.percentiles[0])
        print("10% percentile:            ", self.percentiles[1])
        print("25% percentile:            ", self.percentiles[4])
        print("75% percentile:            ", self.percentiles[14])
        print("90% percentile:            ", self.percentiles[17])                
        print("95% percentile:            ", self.percentiles[18])
        print("Mean:                      ", self.mean)
        print("Variance (bias):           ", self.var_bias)
        print("Variance:                  ", self.var)
        print("Standard deviation (bias): ", self.stdv_bias)
        print("Standard deviation:        ", self.stdv)          
        print("Skewness (bias):           ", self.skewness_bias)
        print("Skewness:                  ", self.skewness)
        print("Excess kurtosis (bias):    ", self.kurtosis_bias)
        print("Excess kurtosis:           ", self.kurtosis)
                
        if ((self.n_intervals is not None) and (self.n_intervals > 0)):
           print("Number of bins:            ", self.n_intervals)

        return (status, message)

        
    def expected_frequencies(self, pdf, n_parameters=1, intervals=0,
                             int_method='sqrt', method='iterate',
                             debug=False):
        """ Calculates the expected frequencies according to a given model.

            The intervals are calculated in order to satisfy 
            some conditions, depending on the method chosen
            * method = 'fixed':    the given number of bins is used, 
                                   and they are equiparted
            * method = 'tails':    width of each interval is a multiple 
                                   of self.width/sqrt(self.n_data)
                                   expected frequency for each interval 
                                   (except the last one) 
                                   is not less than 5                      
            * method = 'iterate':  all intervals have the same width
                                   none of the intervals has 
                                       zero frequency 
                                   not more than 20% of the intervals 
                                       has a frequency lower than 5

            Parameters
            ----------

            pdf:            probability distribution function describing 
                            the model to which the experimental data must 
                            be compared
            n_parameters:   number of parameters of the pdf which were 
                            estimated (used to calculate the degrees 
                            of freedom)
            intervals:      maximum allowed number of intervals 
                            (0 means no limit given)
            int_method:     method used to calculate the maximum number 
                            of intervals, used only if intervals=0
                            'sqrt':  n_bins = sqrt(n_points)
                            'log2':  n_bins = 1 + log2(n_points)
                            'root3': n_bins = n_points**(1/3)
                            'norm':  n_bins = (n_data 
                                               / (3.5 
                                                  * stdv 
                                                  / n_points**(1/3) )
            method:         method used to calculate the intervals
                            'fixed':    a fixed number of equiparted 
                                        intervals is used
                            'tails':    width of each interval is a 
                                        multiple of 
                                        self.width/sqrt(self.n_data)
                                        expected frequency for each 
                                        interval (except the last one)
                                        is not less than 5                      
                            'iterate':  all intervals have the same width
                                        none of the intervals has 
                                        zero frequency 
                                        no more than 20% of the intervals 
                                        has a frequency lower than 5
            debug:          if set to True, enables debug mode


            Initialized data attributes
            ---------------------------

            self.histogram:         histograms of observed and 
                                    expected frequencies
                self.histogram[0]:  upper limits of the intervals (bins)
                self.histogram[1]:  expected frequencies
                self.histogram[2]:  observed frequencies, 
                                    initialized to zero
            self.n_intervals:       number of intervals of the histogram 
                                    (will be set to None on error)
            self.freedom_degrees:   number of degrees of freedom 
                                    = self.n_intervals - 1 - n_parameters

            Returns
            -------

            status: 0 = no error
                    1 = number of intervals is negative
                    2 = number of intervals is greater than the number of points
                    3 = unknown max bin number calculation method
                    4 = unknown bin calculation method
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK
      
        # If not all data are equal and the number of data is
        # at least 10, calculate the number of classes
        if ( (self.width > 0.0) & (self.n_data >= 10) ): 
            # If the number of intervals is negative, raise an error
            if (intervals < 0):
                status = 1
                message = 'Number of intervals must be positive'
                return (status, message)
            # If the number of intervals is greater than
            # the number of points, raise an error
            elif (intervals >= self.n_data):
                status  = 2
                message = ('Number of intervals must be lower '
                           + 'than the number of points')
                return (status, message)
            # If the maximum number of intervals was not given,
            # calculate a convenient value
            elif (intervals == 0):
                 self.n_intervals = calculate_bins_number(
                                        n_data=self.n_data,
                                        method=int_method,
                                        width=self.width,
                                        stdv=self.stdv)
                 if (self.n_intervals is None):
                    status  = 3
                    message = ('Unknown max bin number '
                               + 'calculation method \"'
                               + int_method+'\"')
                    return (status,message)
            else:
                self.n_intervals = intervals

            # Calculate intervals and frequencies
            if (method == 'fixed'):
                (status, message) = self._expected_frequencies_fixed(
                                        pdf, debug)
            elif (method == 'tails'):
                (status, message) = self._expected_frequencies_tails(
                                        pdf, debug)
            elif (method == 'iterate'):
                (status, message) = self._expected_frequencies_iterate(
                                        pdf, debug)
            else:
                self.n_intervals = None
                status  = 4
                message = ('Unknown bin calculation method \"'
                           + method + '\"')
                return (status,message)

            # Calculate the number of degrees of freedom
            self.freedom_degrees = self.n_intervals - 1 - n_parameters 

        # If alla data are equal, or the number of data is less
        # than 10, than make only one class            
        else:   
            self.n_intervals = 1
            self.histogram = numpy.array( [
                                           [0.0],
                                           [self.n_data],
                                           [self.n_data] ] )
            self.freedom_degrees = 0

        self.expected_ready = True
                               
        return (status, message)

        
    def _expected_frequencies_fixed(self, pdf, debug=False):
        """ This method should not be used directly, 
            but by calling expected_frequencies(method='fixed').

            Calculates expected frequencies for the histogram with a
            fixed number of equiparted bins

            Parameters
            ----------

            pdf:    probability distribution function describing 
                    the model to which experimental data should 
                    be compared
            debug:  if set to True, enables debug mode

            Returns
            -------

            (status, message)
            status:  0 = no error
            message: a string containing an error message or 'Ok'
        """

        status , message = 0, OK

        # Initialize the histogram array
        self.histogram = numpy.zeros((3, self.n_intervals))
                
        # Calculate the width of histogram intervals
        delta = self.width / self.n_intervals
        
        if debug: frequency_sum = 0
        for i in range(self.n_intervals):
            left_margin      = self.min + i * delta
            right_margin     = left_margin + delta
            (frequency, err) = quad(pdf, left_margin, right_margin)
            frequency        = self.n_data * frequency
            self.histogram[0,i] = right_margin
            self.histogram[1,i] = frequency
                        
            if debug:
                frequency_sum = frequency_sum + frequency
                print('')
                print('Left  boundary:       ', left_margin)
                print('Right boundary:       ', right_margin)
                print('Frequency:            ', frequency)
                print('Error:                ', err)
                print('Cumulative frequency: ', frequency_sum)
                print('Available particles:  ', self.n_data-frequency_sum)
                input('\nPress RETURN to continue...\n')
                                
        return (status, message)
        

    def _expected_frequencies_tails(self, pdf, debug=False):
        """ This method should not be used directly, 
            but by calling expected_frequencies(method='tails').

            WARNING: this method works only for *positive* data: 
            no negative values are allowed

            Calculates expected frequencies for the histogram 
            according to these rules
            * width of each interval is a multiple of 
              self.width/sqrt(self.n_data)
            * the expected frequency is at least 5 for each interval 
              (excluding the last one)

            Parameters
            ----------

            pdf:    probability distribution function describing 
                    the model to which experimental data must be compared 
            debug:  if set to True, enables debug mode

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = negative data found
            message: a string containing an error message or 'Ok'
        """

        status , message = 0, OK

        # Check that all data values are positive
        if (self.min < 0):
            status  = 1
            message = 'All data values must be positive'
            return (status, message)

        # We use lists instead of numpy arrays, since we don't know
        # the number of intervals (yet)
        interval_list   = []
        frequency_list  = []

        # Calculate intervals and expected frequencies

        # Calculate the width of histogram intervals
        delta = self.width / self.n_intervals   

        # Calculate the sum of the absolute expected frequencies
        # of all the intervals
                               
        # Initialize to the expected frequency between 0 and data_min
        frequency_sum = self.n_data * quad(pdf, 0.0, self.min)[0]

        # Calculate the boundaries of the intervals
        i = 0
        while (self.min+i*delta < self.max) & (frequency_sum < self.n_data - 5):

            # Define the first try for the interval boundaries

            if debug: print('\nInterval # ' + str(i) + '\n')

            if (i==0):
                # Left boundary of first interval is the min data value            
                left_margin = self.min  
            else:
                left_margin = interval_list[i-1]
            # Set the right boundary equal to the left boundary                   
            right_margin = left_margin

            # Increase the right boundary by delta and calculate
            # the expected frequency based on the model
            # repeat until the expected frequency is at least 5
            frequency = 0                               
            while (frequency < 5.0):
                right_margin        = right_margin + delta
                (frequency, err)    = quad(pdf, left_margin, right_margin)
                # convert relative to absolute frequency
                frequency           = self.n_data * frequency
                err                 = self.n_data * err
                if debug: print( '[ ' + str(left_margin) + ', '
                                 + str(right_margin)+' ]', frequency)

            # Now we are sure that the expected frequency
            # for this interval is at least 5
            interval_list.append(right_margin)
            frequency_list.append(frequency)
            frequency_sum = frequency_sum + frequency
            if debug: 
                print('')
                print('Left  boundary:       ', left_margin)
                print('Right boundary:       ', right_margin)
                print('Frequency:            ', frequency_list[i])
                print('Error:                ', err)
                print('Cumulative frequency: ', frequency_sum)
                print('Available particles:  ', self.n_data-frequency_sum)
                input('\nPress RETURN to continue...\n')
            i = i + 1

        # Now we move the upper (right) boundary of the last interval
        # until the total frequency is (nearly) equal to the number
        # of particles
        frequency_sum = frequency_sum - frequency
        tolerance = self.n_data / 100.0
        while (frequency_sum + frequency < self.n_data - tolerance):
            right_margin = right_margin + delta
            frequency    = (self.n_data
                            * quad(pdf, left_margin, right_margin)[0])
            if debug: 
                print( '[ ' + str(left_margin) + ', '
                       + str(right_margin) + ' ]', frequency)
        interval_list[i-1]  = right_margin
        frequency_list[i-1] = frequency
        frequency_sum = frequency_sum + frequency

        # Now that we know the number of intervals,
        # we can define a numpy array with 3 rows
        # row 0: intervals upper boundaries
        # row 1: expected frequencies
        # row 2: observed frequencies (to be calculated later,
        #         by calling method observed_frequencies())
        self.n_intervals = len(frequency_list)
        self.histogram   = numpy.zeros((3, self.n_intervals))
        for i in range(self.n_intervals):
            self.histogram[0,i] = interval_list[i]
            self.histogram[1,i] = frequency_list[i]

        return (status, message)

    def _expected_frequencies_iterate(self, pdf, debug=False):
        """ This method should not be used directly, 
            but by calling expected_frequencies(method='iterate').

            Calculates the expected frequencies for the histogram 
            according to the following rules
            * not more than 20% of the intervals has a frequency 
              lower than 5
            * none of the intervals has zero frequency

            Parameters
            ----------

            pdf:   probability distribution function describing the model 
                   to which experimental data have to be compared                   

            debug: if set to True, enables debug mode

            Returns
            -------

            (status, message)
            status:  0 = no error
            message: a string containing an error message or 'Ok'
        """
                      
        status, message = 0, OK

        # Calculate intervals and expected frequencies
        while True:

            frequencies = numpy.zeros(self.n_intervals)
            intervals   = numpy.zeros(self.n_intervals)
            delta       = self.width / self.n_intervals

            if debug: 
                print('\nNumber of intervals:'+str(self.n_intervals))
                print('Min value:          '+str(self.min))
                print('Max value:          '+str(self.max))
                print('Width:              '+str(self.width))
                print('Delta:              '+str(delta))
                print('\n') 

            less_than_5     = 0
            zeros           = 0
            # Calculate expected frequencies
            for i in range(self.n_intervals):
                intervals[i] = self.min + (i+1)*delta
                if (i==0):
                    left_margin = self.min
                else:
                    left_margin = intervals[i-1]
                right_margin = intervals[i]
                tol = 1.0E-1
                (frequency_rel, err) = quad(pdf, left_margin, right_margin)
                frequency = self.n_data * frequency_rel
                err = self.n_data * err
                if (frequency < 5):
                    less_than_5 = less_than_5 + 1
                if (frequency < tol):
                    zeros = zeros + 1
                frequencies[i] = frequency
                if debug: 
                    print('[ '   +str(left_margin)+', '+str(right_margin)+' ]',)
                    print('f: '  +str(frequency),) 
                    print('f<5: '+str(less_than_5),)
                    print('f=0: '+str(zeros),)
                    print('err: '+str(err))
            # end of for i in range(self.n_intervals):
                      
            if debug: 
                print( 'n: ' + str(self.n_intervals) + '; f<5: '
                       + str(less_than_5) + '; f=0: ' + str(zeros)
                       + '; sum: ' + str(frequencies.sum()) + '\n')
                input('\nPress RETURN to continue...\n')

            # Control if frequencies satisfy conditions, otherwise
            # reduce intervals number and recalculate
            if ( (self.n_intervals <= 1)
                 | ( (float(less_than_5)/float(self.n_intervals) < 0.2)
                     & (zeros==0) ) ): 
                break
            else:
                self.n_intervals = self.n_intervals - 1                 

        # end of While True
                      
        # Now that we now the number of intervals, we can define
        # a numpy array with 3 rows
        # row 0: intervals upper boundaries
        # row 1: expected frequencies
        # row 2: observed frequencies (to be calculated by calling
        #        observed_frequencies() method)
        self.histogram = numpy.zeros((3, self.n_intervals))
        for i in range(self.n_intervals):
            self.histogram[0,i] = intervals[i]
            self.histogram[1,i] = frequencies[i]

        return (status, message)
                

    def observed_frequencies(self, int_method='sqrt', debug=False):
        """ Calculates the observed frequencies based on the raw data.

            If the number of intervals was not defined yet 
            (e.g. by calling expected_frequencies() ) it is calculated

            Parameters
            ----------

            int_method: method used to calculate the maximum number 
                        of intervals
                        'sqrt':  n_bins = sqrt(n_points)
                        'log2':  n_bins = 1 + log2(n_points)
                        'root3': n_bins = n_points**(1/3)
                        'norm':  n_bins = (n_data 
                                           / (3.5 * stdv 
                                              / n_points**(1/3) ))
            debug:      if set to True, enables debug mode

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = error during interval calculation
                     message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        # If the intervals were not defined, than calculate them now
        if (self.n_intervals is None): 
            self.n_intervals = calculate_bins_number(n_data=self.n_data,
                                                     method=int_method,
                                                     width=self.width,
                                                     stdv=self.stdv)
            if (self.n_intervals is None):
                status = 1
                message = 'unable to compute histogram intervals'
                return (status, message)
            self.histogram = numpy.zeros((3, self.n_intervals))
            delta = self.width / self.n_intervals
            for i in range(self.n_intervals):
                 self.histogram[0,i] = self.min + (i+1) * delta
        # end of if (self.n_intervals is None):
                      
        # Fill the data frequencies array
        for i in range(self.n_data):
            # Find the interval to which the data belongs and increase
            # the frequency of the correct interval
            for j in range(self.n_intervals):
                if (j==0):      inf = self.min
                else:           inf = self.histogram[0,j-1]
                sup = self.histogram[0,j]                                       
                if (self.data[i] >= inf) & (self.data[i] < sup):
                   self.histogram[2,j] = self.histogram[2,j] + 1
                   break
               
        self.observed_ready = True

        return (status, message)
    

    def chisquare_estimation(self):
        """ Calculates the chisquare estimation based on the expected 
            and observed frequencies.

            WARNING: you *must* call the methods expected_frequencies 
            and observed_frequencies (the order matters) *before* 
            calling this method. 

            Initialized data attributes
            ---------------------------

            self.chisquare: chisquare estimation
            self.p_value:   p-value associated to the chisquare estimation
                            This is the integral of the chisquare pdf 
                            calculated between the value of the 
                            chisquare estimation and infinite
                            You would say that the differences between 
                            expected and observed frequencies are 
                            statistically significant if this number 
                            is lower than a fixed significance level 
                            (usually 0.05)

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = the self.histogram attribute was not defined yet
                     message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        error_message   = ('You must calculate expected and observed '
                           +' frequencies before calling this method')
        try:
            D = (self.histogram[2]-self.histogram[1])**2 / self.histogram[1]
        except AttributeError:
            status  = 1
            message = error_message
            return (status, message)

        self.chisquare = D.sum()
        try:
            self.p_value   = chi2.sf(self.chisquare, self.freedom_degrees)
        except AttributeError:
            status  = 1
            message = error_message
            return (status, message)                       

        return (status, message)

                      
    def _autocorrelation(self, lag, alternate_method=False):
        """ Calculates the autocorrelation coefficient for a given lag value.

            Parameters
            ----------

            lag: lag value

            Returns
            -------

            autocorrelation: autocorrelation coefficient
        """

        # Calculate autocovariance of data for the given lag
        C_lag = (self.residuals[0:self.n_data-1-lag]
                 * self.residuals[lag:self.n_data-1]).sum()
        if alternate_method:
            C_lag = C_lag / (self.n_data - lag)
        else:
            C_lag = C_lag / self.n_data
                
        # Calculate autocorrelation coefficient
        return C_lag / self.var_bias        


    def plot_histogram(self, plot_expected=True, plot_observed=True,
                       interface=PLOT_INTERFACE, close=False):
        """ Plots the histograms of expected and/or observed frequencies.

            Parameters
            ----------
                        
            interface:     'gnuplot' or 'pylab'
            plot_expected: if True, plot expected frequencies histogram
            plot_observed: if True, plot observed frequencies
            Either plot_observed, plot_expected or both must be True
            otherwise an error is raised

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = no frequency to plot was selected
                     2 = expected frequencies not available yet
                     3 = observed frequencies not available yet
                     4 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        if not (plot_expected or plot_observed):
            status  = 1
            message = "no frequency to plot was selected"
            return (status, message)
        if (plot_expected and not self.expected_ready):
            status = 2
            message = "expected frequencies not available yet"
            return (status, message)
        if (plot_observed and not self.observed_ready):
            status = 3
            message = "observed frequencies not available yet"
            return (status, message)        

        if (interface=='gnuplot'):
            # Define a new histogram, adding the left margin of the
            # first bin (nicer plot)
            left_margin = numpy.array( [ [self.min],
                                         [self.histogram[1,0]],
                                         [self.histogram[2,0]] ] )
            histogram_ext = numpy.hstack((left_margin, self.histogram))
            histogram_plot = gm.new_plot(persistence=True,
                                         style='fsteps',
                                         title='Data histogram',
                                         xlabel='Data intervals',
                                         ylabel='Absolute frequencies')
            if (plot_observed):
                gm.plot2d(histogram_plot,
                          histogram_ext[0],
                          histogram_ext[2], 
                          label='Observed frequencies')
                if (plot_expected):
                    gm.plot2d(histogram_plot,
                              histogram_ext[0],
                              histogram_ext[1],
                              label='Expected frequencies',
                              replot=True)
            elif (plot_expected):
                gm.plot2d(histogram_plot,
                          histogram_ext[0],
                          histogram_ext[1],
                          label='Expected frequencies')
            if close: gm.plot_close(histogram_plot)
                      
        elif (interface=='pylab'):
            # Define a new histogram
            histogram_plot = numpy.zeros((3,self.n_intervals+1))
            for i in range(self.n_intervals+1):
                if (i==0): histogram_plot[0,i] = self.min
                else:      histogram_plot[0,i] = self.histogram[0,i-1]
            for i in range(self.n_intervals+1):
                if (i==self.n_intervals): 
                    histogram_plot[1,i] = self.histogram[1,i-1]
                    histogram_plot[2,i] = self.histogram[2,i-1]
                else:                     
                    histogram_plot[1,i] = self.histogram[1,i]
                    histogram_plot[2,i] = self.histogram[2,i]
            plt.title('Data histogram')
            plt.xlabel('Data intervals')
            plt.ylabel('Absolute frequencies')
            if (plot_expected):
                plt.plot(histogram_plot[0],
                           histogram_plot[1],
                           #linestyle='steps',
                           drawstyle='steps',
                           color='blue',
                           label='Expected frequencies')
            if (plot_observed):
                plt.plot(histogram_plot[0],
                           histogram_plot[2],
                           #linestyle='steps',
                           drawstyle='steps',
                           color='red', 
                           label='Observed frequencies')

            plt.legend()
            plt.show()
                      
        else:
            status = 4
            message = 'Unknown plot interface \"' + interface + '\"'
            return (status, message)

        return (status, message)


    def plot_lag(self, lag=1, interface=PLOT_INTERFACE, symbol=PLOT_SYMBOL,
                 linestyle=PLOT_LINE, color=PLOT_COLOR, close=False):
        """ Plots a lag plot of the data.

            Parameters
            ----------

            lag:       lag
            interface: 'gnuplot' or 'pylab'

            Returns
            -------

            status:  0 = no error
                     1 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        lag = int(lag)
        if (lag > self.n_data-2):
            status  = 1
            message = 'Lag exceeds maximum allowed value'
            return (status, message)
                                
        x_values = self.data[0:self.n_data-1-lag]
        y_values = self.data[lag:self.n_data-1]
        title = 'Lag plot (lag='+str(lag)+')'

        if (interface=='gnuplot'):      
            lag_plot = gm.new_plot(persistence=True,
                                   title=title,
                                   xlabel='X(i-' + str(lag) + ')',
                                   ylabel='X(i)')
            gm.plot2d(lag_plot, x_values, y_values)
            if close: gm.plot_close(lag_plot)
        elif (interface=='pylab'):
            plt.title(title)
            plt.xlabel('X(i-lag)')
            plt.ylabel('X(i)')
            #plt.scatter(x_values, y_values,
            #          linestyle=linestyle,, marker=symbol, color=color)
            plt.plot(x_values, y_values, marker=symbol,
                       linestyle=linestyle, color=color)
            plt.show()
        else:
            status = 1
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)
                      
        return (status, message)

    def plot_autocorrelation(self, start=0, step=1, n_points=None,
                             alternate_method=False,
                             interface=PLOT_INTERFACE, symbol=PLOT_SYMBOL,
                             linestyle=PLOT_LINE, color=PLOT_COLOR, close=False):
                
        """ Plots an autocorrelation plot of the data.

            Parameters
            ----------
                        
            interface:      'gnuplot' or 'pylab'
            start:          start the plot from the data value with this index
            step:           plot only every step values
            n_points:       number of points to plot
                        
            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = start is negative or greater than (or equal to) 
                         the number of data points
                     2 = step is zero, negative or greater than the last-first
                     3 = n_points is lower than first point or greater 
                         than the number of data points
                     4 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """
        status, message = 0, OK

        start = int(start)
        if ((start < 0) or (start >= self.n_data)):
            status  = 1
            message = 'start must be in range [0..' + str(self.n_data-1) + ']'
            return (status, message)
        step = int(step)
        if ((step < 1) or (step > (self.n_data-start))):
            status = 2
            message = 'step must be in range [1..'+str(self.n_data-start)+']'
            return (status, message)
        if (n_points is None):
            n_points = int( math.ceil( 1.0 * (self.n_data-start) / step) )
        else:
            n_points = int(n_points)
            if (n_points <= 0):
                status = 3
                message = 'n_points must be strictly positive'
                return (status, message)
            elif (start + n_points * step >= self.n_data ):
                n_points = int( math.ceil( 1.0 * (self.n_data-start) / step) )

        last_point = start + step * (n_points - 1)
        # This check should be unnecessary but...
        if (last_point > self.n_data-1): last_point = self.n_data-1
                
        lag_list = range(start, last_point+1, step)
        autocorrelation_array = numpy.zeros(len(lag_list))
        for i in range(len(lag_list)):
            autocorrelation_array[i] = self._autocorrelation(
                lag=lag_list[i], alternate_method=alternate_method)
            
        if (interface=='gnuplot'):
            autocorrelation_plot = gm.new_plot(persistence=True,
                                               title='Autocorrelation plot',
                                               xlabel='lag',
                                               ylabel='Autocorrelation')
            gm.plot2d(autocorrelation_plot, numpy.array(lag_list), autocorrelation_array)
            if close: gm.plot_close(autocorrelation_plot)

        elif (interface=='pylab'):
            plt.grid()
            plt.title('Autocorrelation plot')
            plt.xlabel('lag')
            plt.ylabel('Autocorrelation')
            plt.plot(numpy.array(lag_list),
                       autocorrelation_array,
                       marker=symbol,
                       linestyle=linestyle,
                       color=color)
            plt.show()
        else:
            status = 5
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)

        return (status, message)


    def plot_runsequence(self, mean=False, n_stdv=None, interface=PLOT_INTERFACE,
                         symbol=PLOT_SYMBOL, color=PLOT_COLOR,
                         color_mean=PLOT_COLORS[0], color_stdv=PLOT_COLORS[1], close=False):
        """ Plots a runsequence plot of the data.

            Parameters
            ----------

            interface: 'gnuplot' or 'pylab'
            mean:      if True, plot a line representing the mean value 
                       of the data
            n_stdv:    if not None, plot lines representing n_stdv times 
                       the standard deviation

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = unknown plot interface
                     2 = negative or zero n_stdv value was given
                   
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        if (mean):
            meanline = self.mean * numpy.ones(self.n_data)                
        if (n_stdv is not None):
            if (n_stdv <= 0):
                status = 2
                message = 'n_stdv must be strictly positive'
                return (status, message)
            else:
                stdvline_h = ((self.mean + n_stdv * self.stdv)
                              * numpy.ones(self.n_data))
                stdvline_l = ((self.mean - n_stdv * self.stdv)
                              * numpy.ones(self.n_data))
                str_h      = 'mean+' + (n_stdv!=1) * str(n_stdv) + 'stdv'
                str_l      = 'mean-' + (n_stdv!=1) * str(n_stdv) + 'stdv'

        if (interface=='gnuplot'):      
            # Define a new plot window and set titles
            runsequence_plot = gm.new_plot(persistence=True,
                                           title='Runsequence plot',
                                           xlabel='i',
                                           ylabel='X(i)')

            l = []
            l.append([None, self.data, None, 'with points'])
            if mean:
                l.append([None, meanline, 'mean', 'with lines'])
            if (n_stdv is not None):
                l.append([None, stdvline_h, str_h, 'with lines'])
                l.append([None, stdvline_l, str_l, 'with lines'])
            gm.plot_curves(runsequence_plot, l)
            #gm.plot1d(runsequence_plot, self.data)
            #if (mean):
            #    gm.plot1d(runsequence_plot, meanline, label='mean', replot=True)
            #if (n_stdv is not None):
            #    gm.plot1d(runsequence_plot, stdvline_h, label=str_h, replot=True)
            #    gm.plot1d(runsequence_plot, stdvline_l, label=str_l, replot=True)
            if close: gm.plot_close(runsequence_plot)

        elif (interface=='pylab'):
            # Plot the data using matplotlib
            plt.title('Runsequence plot')
            plt.xlabel('i')
            plt.ylabel('X(i)')
            plt.plot(self.data, marker=symbol, color=color, linestyle='None')
            if (mean):
                plt.plot(meanline, marker=None, color=color_mean,
                           linestyle='-', label='mean')
            if (n_stdv is not None):
                plt.plot(stdvline_h, marker=None, color=color_stdv,
                           linestyle='--', label=str_h)
                plt.plot(stdvline_l, marker=None, color=color_stdv,
                           linestyle='--', label=str_l)
            if (mean or (n_stdv is not None)):
                plt.legend()         
            plt.show()
        else:
            status = 1
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)

        return (status, message)

    def plot_residuals(self, interface=PLOT_INTERFACE, symbol=PLOT_SYMBOL,
                       color=PLOT_COLOR, close=False):
        """ Plots the residuals of data.

            Parameters
            ----------

            interface: 'gnuplot' or 'pylab'

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         
                        
        if (interface=='gnuplot'):      
            # Define a new plot window and set titles
            residuals_plot = gm.new_plot(persistence=True,
                                         title='Residuals plot',
                                         xlabel='i',
                                         ylabel='X(i)-<X>')
            gm.plot1d(residuals_plot, self.residuals)
            if close: gm.plot_close(residuals_plot)
                      
        elif (interface=='pylab'):
            # Plot the histograms using matplotlib
            plt.title('Residuals plot')
            plt.xlabel('i')
            plt.ylabel('X(i)-<X>')
            plt.plot(self.residuals, marker=symbol, color=color,
                       linestyle='None')
            plt.show()
                      
        else:
            status = 1
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)

        return (status, message)


    def plot_residuals_lag(self, lag=1, interface=PLOT_INTERFACE,
                           symbol=PLOT_SYMBOL, color=PLOT_COLOR, close=False):
        """ Plots a lag plot of the data residuals.

            Parameters
            ----------

            lag:       lag
            interface: 'gnuplot' or 'pylab'

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        if (lag > self.n_data-2):
            status  = 1
            message = 'Lag exceeds maximum allowed value'
            return (status, message)
                                
        x_values = self.residuals[0:self.n_data-1-lag]
        y_values = self.residuals[lag:self.n_data-1]
        title = 'Lag plot (lag='+str(lag)+')'

        if (interface=='gnuplot'):      
            # Define a new plot window and set titles
            lag_plot = gm.new_plot(persistence=True,
                                   title=title,
                                   xlabel='X(i-lag)',
                                   ylabel='X(i)')
            gm.plot2d(lag_plot, x_values, y_values )
            if close: gm.plot_close(lag_plot)
                      
        elif (interface=='pylab'):
            # Plot the histograms using matplotlib
            plt.title(title)
            plt.xlabel('X(i-lag)')
            plt.ylabel('X(i)')
            plt.scatter(x_values, y_values, marker=symbol, color=color)
            plt.show()
                      
        else:
            status = 1
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)

        return (status, message)
        

    def plot_box(self, whiskers='range', interface=PLOT_INTERFACE,
                 grid=True, color=PLOT_COLOR, close=False):
        """ Plots a boxplot of the data.

            Parameters
            ----------

            interface: 'gnuplot' or 'pylab'
            whiskers:  extension of the whiskers
                       - if set to 'range' (default) forces the 
                         whiskers to match the highest and lowest 
                         values of the data
                       - if set to a number
                         top whisker    is Q3 + whiskers * (Q3-Q1)
                         bottom whisker is Q1 - whiskers * (Q3-Q1)
            grid:      if set to True, plots a grid on the graph
            color:     color to use in the plot

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = non numerical parameter given
                     2 = unknown plot interface
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        title = 'Box plot'

        if (whiskers!='range'):
            try:
                whiskers = round(whiskers, 2)
            except TypeError:
                status = 1
                message = 'Parameter \"'+ whiskers + '\" is not a number'
                return (status, message)

        if (interface=='gnuplot'):
            width = int(self.n_data/2.0)          
            box_plot = gm.new_plot(persistence=True,
                                   title=title,
                                   xlabel='Data',
                                   ylabel='Min, Q1, Median, Q3, Max')
            #gm.plot_command( box_plot, 'set boxwidth ' + str(width) )
            if (whiskers=='range'):
                gm.plot_command(box_plot, 'set style boxplot fraction 1')
            else:
                gm.plot_command(box_plot, 'set style boxplot range ' + str(whiskers) )
            gm.plot_box(box_plot, self.data, width=width)
            if close: gm.plot_close(box_plot)

        elif (interface=='pylab'):
            # Plot the boxplot using matplotlib
            plt.title(title)
            plt.ylabel('Min, Q1, Median, Q3, Max')
            plt.grid()
            plt.boxplot(self.data, whis=whiskers, showmeans=True,
                        meanline=True)
            plt.show()

        else:
            status = 1
            message = 'Unknown plot interface' + '\"'+interface + '\"'
            return (status, message)

        return (status, message)

    
    def plot_clear(self):
        """ Close all the gnuplot open plots. """

        return gm.plot_close_all(purge=True)
    
