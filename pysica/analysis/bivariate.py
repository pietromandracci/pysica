# COPYRIGHT (c) 2020-2024 Pietro Mandracci

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

""" Tools for simple statistical analysis on bivariate experimental data.

    This module contains some classes, methods and functions,useful to 
    perform some basic analysis on a bivariate set of experimental data.

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
#import matplotlib.patches as patches 

#from scipy.integrate   import quad
#from scipy.stats       import chi2

# Mudules from plasma.ccpla package
from pysica.parameters import *
from pysica.managers import data_manager

# +---------------+ 
# |Bivariate data |
# +---------------+

class DataBivariate:
    """ Class defining a bivariate set of experimental data to be analyzed. """

    def __init__(self, data_matrix=None, filename='', name_x='x', name_y='y', skip=0):
        """ Defines an instance of the class DataBivariate 
            and calculates some basic statistical parameters.

            data_matrix:     a numpy 2d array: if not given, will be read from a file
            filename:        name of the file from which the data must be read, 
                             if data_matrix is not given, this must not be blank
            name_x:          a string used to identify the x values
            name_y:          a string used to identify the y values
            skip:            used only if data_matrix is not given
                             number of rows to skip while reading data from the file
            self.read_error: used for error checking
                             (status, message)
                             status: 0 = no error
                                     1 = could not open file
                                     2 = syntax error
                                     3 = non-numerical value
                                     4 = invalid number of rows/columns
                                     5 = invalid number of rows to skip
                                     6 = non-array argument was given
                                     7 = the array is not 2-d
                                     8 = wrong array structure read from file
                             message: an error message or 'OK'
        """

        status, message = 0, OK 

        if (data_matrix == None):
            matrix = data_manager.DataGrid()
            self.read_error = matrix.read_file(filename, transpose=True, skip=skip)
            if (self.read_error[0] != 0):
                status  = self.read_error[0]
                message = self.read_error[1]
                self.read_error = (status, message)                
                return
            self.data = matrix.data_array
            del matrix
            n = len(self.data)
            if (n != 2):
                status = 8
                message = 'wrong array structure in file \"'+filename+' \"'
                self.read_error = (status, message)
                return
        else:   
            try:
                n = len(data_matrix.shape)
            except AttributeError:
                status = 6
                message = 'non-array argument given'
                self.read_error = (status, message)
                return
            if (n == 2): 
                self.data = data_matrix
            else:
                status = 7
                message = 'array must be 2-dimensional'
                self.read_error = (status, message)
                return

        self.name_x =   name_x
        self.name_y =   name_y
            
        self.n_data_x = len(self.data[0])
        self.n_data_y = len(self.data[1])

        self.max_x =    self.data[0].max()        
        self.max_y =    self.data[1].max()

        self.min_x =    self.data[0].min()        
        self.min_y =    self.data[1].min()

        self.width_x =  self.max_x - self.min_x
        self.width_y =  self.max_y - self.min_y

        # Median and percentiles
        # numpy.linspace(5,100,20) = array of percentile values (5% step)
        self.percentiles_x = numpy.percentile(self.data[0],
                                              numpy.linspace(5,100,20))
        self.percentiles_y = numpy.percentile(self.data[1],
                                              numpy.linspace(5,100,20))

        
        self.median_x = self.percentiles_x[9]
        self.median_y = self.percentiles_y[9]
        
        self.mean_x =   numpy.mean(self.data[0])
        self.mean_y =   numpy.mean(self.data[1])

        # Residuals and their powers
        self.residuals_x   = self.data[0] - self.mean_x
        self.residuals_2_x = self.residuals_x * self.residuals_x
        self.residuals_3_x = self.residuals_2_x * self.residuals_x
        self.residuals_4_x = self.residuals_3_x * self.residuals_x

        # Residuals and their powers
        self.residuals_y   = self.data[1] - self.mean_y
        self.residuals_2_y = self.residuals_y * self.residuals_y
        self.residuals_3_y = self.residuals_2_y * self.residuals_y
        self.residuals_4_y = self.residuals_3_y * self.residuals_y        
                
     
        # Variance of data, with and without bias
        self.var_bias_x = self.residuals_2_x.sum() / self.n_data_x
        self.var_x      = self.var_bias_x * self.n_data_x / (self.n_data_x-1)

        self.var_bias_y = self.residuals_2_y.sum() / self.n_data_y
        self.var_y      = self.var_bias_y * self.n_data_x / (self.n_data_y-1)

        # Standard deviation, with and without bias
        self.stdv_x      = numpy.sqrt(self.var_x)
        self.stdv_y      = numpy.sqrt(self.var_y)
        self.stdv_bias_x = numpy.sqrt(self.var_bias_x)
        self.stdv_bias_y = numpy.sqrt(self.var_bias_y)

        # Standard deviation of the mean
        self.stdv_mean_x = self.stdv_x / numpy.sqrt(self.n_data_x)
        self.stdv_mean_y = self.stdv_y / numpy.sqrt(self.n_data_y)
        

    def print_stats(self):
        """ Prints some basic statistical information about data.

            Returns
            -------

            (status, message)
            status:  0 = no error
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        print("Stats of data: \"" + self.name_x + "\"")
        print("Number of values:               ", self.n_data_x)
        print("Minimum value:                  ", self.min_x)
        print("Maximum value:                  ", self.max_x)
        print("Data range width:               ", self.width_x)
        print("Median:                         ", self.median_x)                
        print("5%  percentile:                 ", self.percentiles_x[0])
        print("10% percentile:                 ", self.percentiles_x[1])
        print("25% percentile:                 ", self.percentiles_x[4])
        print("75% percentile:                 ", self.percentiles_x[14])
        print("90% percentile:                 ", self.percentiles_x[17])                
        print("95% percentile:                 ", self.percentiles_x[18])
        print("Mean:                           ", self.mean_x)
        print("Standard deviation:             ", self.stdv_x)
        print("Standard deviation of the mean: ", self.stdv_mean_x)        

        print("\n")
        
        print("Stats of data: \"" + self.name_y + "\"")
        print("Number of values:               ", self.n_data_y)
        print("Minimum value:                  ", self.min_y)
        print("Maximum value:                  ", self.max_y)
        print("Data range width:               ", self.width_y)
        print("Median:                         ", self.median_y)                
        print("5%  percentile:                 ", self.percentiles_y[0])
        print("10% percentile:                 ", self.percentiles_y[1])
        print("25% percentile:                 ", self.percentiles_y[4])
        print("75% percentile:                 ", self.percentiles_y[14])
        print("90% percentile:                 ", self.percentiles_y[17])                
        print("95% percentile:                 ", self.percentiles_y[18])
        print("Mean:                           ", self.mean_y)
        print("Standard deviation:             ", self.stdv_y)
        print("Standard deviation of the mean: ", self.stdv_mean_y)                

        return (status, message)          

    def plot_scatter(self, logx=False, logy=False, legend='', title='', 
                     color = PLOT_COLOR, line = PLOT_LINE,
                     symbol = PLOT_SYMBOL):
        """ Plots a scatter plot of the data.

            xlabel          string to be used as label for the x axis
            ylabel          string to be used as label for the y axis
            legend          string to be used as legend for the data plotted
            title           string to be used as title
            logx:           set logaritmic x axis
            logy:           set logaritmic y axis
            return:         (status, message)
                            status:  0 = no error
                                     1 = unknown plot interface
                            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        if (title != ''): plt.title(title)
        plt.xlabel(self.name_x)
        plt.ylabel(self.name_y)
        if (logx and logy):
            if (legend == ''):
                plt.loglog(self.data[0], self.data[1], 
                           marker=symbol, linestyle=line, color=color)
            else:
                plt.loglog(self.data[0], self.data[1], 
                           label=legend, marker=symbol, linestyle=line, color=color)
        elif logx:
            if (legend == ''):
                plt.semilogx(self.data[0], self.data[1], 
                             marker=symbol, linestyle=line, color=color)
            else:
                plt.semilogx(self.data[0], self.data[1], 
                             label=legend, marker=symbol, linestyle=line, color=color)
        elif logy:
            if (legend == ''):
                plt.semilogy(self.data[0], self.data[1], 
                             marker=symbol, linestyle=line, color=color)
            else:
                plt.semilogy(self.data[0], self.data[1], 
                             label=legend, marker=symbol, linestyle=line, color=color)
        else:
            if (legend == ''):
                plt.plot(self.data[0], self.data[1], 
                    marker=symbol, linestyle=line, color=color)
            else:
                plt.plot(self.data[0], self.data[1], 
                    label=legend, marker=symbol, linestyle=line, color=color)
        if (legend != ''): plt.legend()               
        plt.show()

        return (status, message)


    def plot_box(self, stdv=False, stdv_mean=False,
                 quartiles=True, deciles=False, maxmin=False,
                 grid=True):
        """ Plots a boxplot of the data.

            Parameters
            ----------

            stdv:      plot the mean with standard deviation
                       as errorbars
            stdv_mean: plot the mean with standard deviation of the mean 
                       as errorbars
            quartiles: plot the median with 25% and 75% as errorbars
            deciles:   plot the median with 10% and 90% as errorbars
            maxmin:    plot the median with min and max values as errorbars
            grid:      if set to True, plots a grid on the graph

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
        plt.title(title)
        plt.xlabel(self.name_x)
        plt.ylabel(self.name_y)
        if grid: plt.grid()
        if stdv:
            self._plot_box(plt,
                           self.mean_x, self.mean_y,
                           self.stdv_x, self.stdv_x,
                           self.stdv_y, self.stdv_y,
                           box=False,   bars=True,
                           bar_color='orange',
                           label='mean,   ' + 'stdv')
        if stdv_mean:
            self._plot_box(plt,
                           self.mean_x,      self.mean_y,
                           self.stdv_mean_x, self.stdv_mean_x,
                           self.stdv_mean_y, self.stdv_mean_y,
                           box=False,        bars=True,
                           bar_color='magenta',
                           label='mean,   ' + 'stdv_mean')     
        if maxmin:
            dx_low   = self.median_x - self.min_x
            dx_upp   = self.max_x    - self.median_x
            dy_low   = self.median_y - self.min_y
            dy_upp   = self.max_y    - self.median_y
            self._plot_box(plt,
                           self.median_x, self.median_y,
                           dx_low,        dx_upp,
                           dy_low,        dy_upp,                       
                           box=False,      bars=True,
                           bar_color='black',
                           label='median, ' + 'min-max')           
        if deciles:
            dx_low   = self.median_x          - self.percentiles_x[1]
            dx_upp   = self.percentiles_x[17] - self.median_x
            dy_low   = self.median_y          - self.percentiles_y[1]
            dy_upp   = self.percentiles_y[17] - self.median_y
            self._plot_box(plt,
                           self.median_x, self.median_y,
                           dx_low, dx_upp,
                           dy_low, dy_upp,                       
                           box=True,         bars=True,
                           line='--',
                           box_color='blue', bar_color='blue',
                           label='median, ' + '10%-90%')
        if quartiles:
            dx_low   = self.median_x          - self.percentiles_x[4]
            dx_upp   = self.percentiles_x[14] - self.median_x
            dy_low   = self.median_y          - self.percentiles_y[4]
            dy_upp   = self.percentiles_y[14] - self.median_y
            self._plot_box(plt,
                           self.median_x, self.median_y,
                           dx_low, dx_upp,
                           dy_low, dy_upp,                       
                           box=True, bars=True,
                           line='--',
                           box_color='green', bar_color='green',
                           label='median, ' + '25%-75%')
        plt.legend()
        plt.show()

        return (status, message)
    

    def _plot_box(self, plt, x, y, dx_low, dx_upp, dy_low, dy_upp,
                  box=False, bars=True,
                  box_color='black', bar_color='black', line='-',
                  label=''):
        """ Draw a data point with error bars and a rectangle

            x:         x-value of the data point
            y:         y-value of the data point
            dx_low:    lower-end error for x
            dx_upp:    upper-end error for x
            dy_low:    lower-end error for y
            dy_low:    upper-end error for y
            box:       draw a rectangle defining the error aera
            bars:      draw error bars 
            box_color: color of the rectangle
            bar_color: color of the errorbars 
            line:      style of the box line
                       {'-', '--', '-.', ':', '', (offset, on-off-seq), ...}
            label:     string to use in the legend
        """

        if box:
            box_origin = (x-dx_low, y-dy_low)
            box = plt.Rectangle(box_origin, dx_low+dx_upp, dy_low+dy_upp,
                                linewidth=1,
                                edgecolor=box_color,
                                linestyle=line,
                                facecolor='none')
            plt.gca().add_patch(box)
        if bars:
            # NOTE: errobars has a strange parameter order:
            # errorbar(x, y, yerr=None, xerr=None,...) yerr is *before* xerr       
            plt.errorbar(x, y,
                         numpy.array([[dy_low, dy_upp]]).T,
                         numpy.array([[dx_low, dx_upp]]).T,
                         marker='s', capsize=10, color=bar_color,
                         label=label)
    

    def calculate_correlation(self):
        """ Calculates the correlation coefficient of data. 

            Initialized data attributes
            ---------------------------

            self.correlation: 

            return: (status, message)
                    status:  0 = no error
                             1 = unknown plot interface
                    message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        return (status, message)        


