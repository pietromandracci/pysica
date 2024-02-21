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

""" Classes and methods to plot, compare and make some simple 
    calculations on bivariate data sets, such as spectra.

    This module contains some classes, methods and functions, 
    that may be useful to perform some basic analysis on a bivariate 
    set of experimental data, such as plotting and comparing experimental 
    spectra.

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
#import pylab
#from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from scipy.integrate   import quad

# Mudules from plasma.ccpla package
from pysica.parameters import *
from pysica.functions.mathematics import *
from pysica.functions.fortran.fmathematics import f_math
from pysica.functions.optics import *
from pysica.managers import data_manager, unit_manager


# +----------------+ 
# | Set of spectra |
# +----------------+

class SpectraSet:
    """ Class defining a set of spectra to be analyzed"""

    def __init__(self, xarray=None, xmin=None, xmax=None, npoints=None):
        """ Defines a class of spectra.

            Defines an instance of the class SpectraSet, which consists 
            of a set of spectra, each one represented by an array of 
            depending variables (y-values), defined on the same set 
            of x-values.

            The set of spectra is contained in the list self.spectra. 
            Each element of this list is itself a list of 3 elements:
            - a string representing the name of the spectrum
            - a numpy array containing the y-values
            - a string describing the type of spectrum
            The x-values are the same for each spectrum and are contained 
            in the numpy array self.xvalues.

            self.spectra[i][0] -> name of the spectrum #i
            self.spectra[i][1] -> numpy array of the y-values of
                                  the spectrum #i
            self.spectra[i][2] -> string describing the type of spectrum
            self.xvalues       -> numpy array of the x-values for all 
                                  the spectra
                
            Parameters
            ----------

            xarray:     a numpy array which is used as a model to define:
                        - the minimum and maximum of the x-values, 
                        - the number of x-values
                        If it is not given, then the x-values are 
                        calculated using the values xmin, xmax, npoints
            xmin:       minimum value of the range to be calculated 
                                    (used only if xarray is not given)
            xmax:       maximum value of the range to be calculated 
                        (used only if xarray is not given)
                        must be greater than xmin (or the range will be left empty)

            npoints:    number of points of the range to be calculated 
                        (used only if xarray is not given)
                        must be at leats 2  (or the range will be left empty)
                                    if it is not an integer, it will be truncated

            Initialized data attributes
            ---------------------------
                                
            self.npoints:      number of x-values
            self.delta:        separation between x-values
            self.xvalues:      a numpy array containing the x-values, 
                               which are equiparted
            self.spectra:      an empty list to which the y-values and 
                               the names of the spectra are added
            self.error:        (status, message) used for error checking
                               status: 0 = no error
                                       1 = missing parameter
                                       2 = max is lower than min
                                       3 = npoints is lower than 2
                               message: an error message or 'OK'
        """

        status, message = 0, OK 

        if (xarray is None):
            if ( (xmin is None) or (xmax is None) or (npoints is None) ):
                status = 1
                message = 'missing parameter'
                self.error = (status, message)
                return
            else:
                # If values of xmin, xmax and npoints were given,
                # calculate the x-values
                if (xmax <= xmin):
                    status = 2
                    message = 'maximum value must be greater than minimum'
                    self.error = (status, message)
                    return
                elif (npoints < 2):     
                    status = 3
                    message = 'number of points is lower than 2'
                    self.error = (status, message)
                    return
                else:   
                    self.npoints = int(npoints)
                    self.xvalues  = numpy.zeros(self.npoints)
                    (self.xvalues, self.delta) = numpy.linspace(
                        xmin, xmax, self.npoints, retstep=True)
        else:
            # If an array of values was given, create an
            # equiparted array, with the same min,max and number
            # of values
            self.npoints = len(xarray)
            self.xvalues = numpy.zeros(self.npoints)
            xmin = xarray.min()
            xmax = xarray.max()
            (self.xvalues, self.delta) = numpy.linspace(
                xmin, xmax, self.npoints, retstep=True)
            
        self.spectra = []
        self.error = (status, message)
               

#   +------------------------------------------------+
#   | Methods to perform basic operations on spectra |
#   +------------------------------------------------+


    def list_spectra(self):
        """ Lists the names of the spectra.

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = no spectra is present
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        if (len(self.spectra) == 0):
            status = 1
            message = "no spectra have been loaded yet"
            return (status, message)
        else:
            print( 'Index', 'Type'.ljust(20), 'Name')
            print( '-----', '----'.ljust(20), '----')
            for i in range(len(self.spectra)):
                print(str(i).rjust(5), self.spectra[i][2].ljust(20), self.spectra[i][0])

        return (status, message)

        
    def rename_spectrum(self, index, name):
        """ Renames a spectrum from the set of spectra.

            Parameters
            ----------

            index: index of the spectrum to be renamed
            name:  string contanining the new name

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = name given is not a string
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)
        if (plt.is_string_like(name) == False):
            status = 2
            message = 'name given is not a string'
            return (status, message)
        self.spectra[index][0] = name

        return (status, message)


    def delete_spectrum(self, index):
        """ Removes a spectrum from the set of spectra.

            Parameters
            ----------

            index: index of the specrum to remove

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)
        else:                   
            self.spectra.pop(index)

        return (status, message)                            

   
        
    def plot_spectrum(self, index, logx=False, logy=False,
                      xlabel='x', ylabel='y', explabel=True, 
                      title=None, grid=True, color=PLOT_COLOR, line=PLOT_LINE,
                      symbol=PLOT_SYMBOL, fill_color=FILL_COLOR):
        """ Plots a scatter plot of a spectrum.

                Parameters
                ----------

                index:          index identifying which spectrum to plot
                logx:           set logaritmic x axis
                logy:           set logaritmic y axis
                xlabel:         string to be used as label for the x axis
                ylabel:         string to be used as label for the y axis
                explabel:       if True, use scientific format for labels if they exceed 100
                title:          string to be used as title
                grid:           wether to show a grid on the plot
                color:          color to be used for plotting
                line:           type of line
                symbol:         type of symbol to use for data points
                fill_color:     color used to fill symbols area

                Returns
                -------

                status:  0 = no error
                         1 = no spectrum has been loaded yet
                         2 = wrong index value
                message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        if (len(self.spectra) == 0):
            status = 1
            message = "no spectra have been loaded yet"
            return (status, message)
        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status  = 2
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        if (title is not None): plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        label = '(' + str(index).rjust(2) + ') ' + self.spectra[index][0]

        plt.grid(grid)

        if explabel: plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))

        # If the spectrum is a histogram, make a bar plot
        if (self.spectra[index][2].startswith('histogram')):
            left   = self.xvalues - self.delta / 2.0
            height = self.spectra[index][1]
            width  = self.delta
            plt.bar(left, height, width=width, log=logy, label=label,
                    color=fill_color, edgecolor=color)
            
        # If the spectrum is not a histogram, plot normally
        else:   
            if (logx and logy):
                plt.loglog(self.xvalues, self.spectra[index][1], label=label, 
                           marker=symbol, linestyle=line, color=color)
            elif logx:
                plt.semilogx(self.xvalues, self.spectra[index][1], label=label, 
                             marker=symbol, linestyle=line, color=color)
            elif logy:
                plt.semilogy(self.xvalues, self.spectra[index][1], label=label, 
                             marker=symbol, linestyle=line, color=color)
            else:
                plt.plot(self.xvalues, self.spectra[index][1], label=label, 
                         marker=symbol, linestyle=line, color=color)
        plt.legend()            
        plt.show()

        return (status, message)

    
    def plot_3D(self, index_list, y0, dy, cstride=1, rstride=1, type='scatter',
                xlabel='x', ylabel='y', zlabel='z', explabel=True, 
                title=None, grid=False,
                 color=PLOT_COLOR, line=PLOT_LINE, symbol=PLOT_SYMBOL, fill_color=FILL_COLOR):
        """ Plots a 3D plot of a set of spectra.

            Parameters
            ----------

            index_list:     list containing the indexes identifying which spectra to plot
            xlabel:         string to be used as label for the x axis
            ylabel:         string to be used as label for the y axis
            zlabel:         string to be used as label for the z axis
            explabel:       if True, use scientific format for labels if they exceed 100
            title:          string to be used as title
            grid:           wether to show a grid on the plot
            color:          color to be used for plotting
            line:           type of line
            symbol:         type of symbol to use for data points
            fill_color:     color used to fill symbols area

            Returns
            -------

            status:  0 = no error
                     1 = no spectrum has been loaded yet
                     2 = wrong index value
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        if (len(self.spectra) == 0):
            status = 1
            message = "no spectra have been loaded yet"
            return (status, message)
        n_indexes = len(index_list)

        if (n_indexes < 2):
            status = 2
            message = "the indexes of two spectra at least must be given"
            return (status, message)
        if (type not in ('scatter', 'wireframe', 'surface')):
            status = 3
            message = "unkown plot type \'" + type + "\'"
            return (status, message)

        for i in range(n_indexes):
            if ( (index_list[i] < 0) or (index_list[i] >= len(self.spectra)) ):
                status  = 4
                message = 'spectrum ' + str(index_list[i]) + ' does not exist'
                return (status, message)
            if (self.spectra[index_list[i]][2].startswith('histogram')):
                status  = 5
                message = 'spectrum ' + str(index_list[i]) + ' is a histogram'
                return (status, message)                        

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')                        
        if (title is not None): ax.title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)                
        if explabel: ax.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))

        if   (type == 'scatter'):
            # Calculate x,y,z values arrays
            nvalues = n_indexes*self.npoints
            x = numpy.zeros(nvalues)
            y = numpy.zeros(nvalues)
            z = numpy.zeros(nvalues)
            k = 0
            for i in range(self.npoints):
                for j in range(n_indexes):
                        x[k] = self.xvalues[i]
                        y[k] = y0 + j * dy
                        z[k] = self.spectra[ index_list[j] ][1][i]
                        k = k + 1
            ax.scatter(x, y, z, marker=symbol, color=color)
            
        elif ( (type == 'wireframe') or (type == 'surface') ):
            # Calculate x,y,z values arrays
            x   = numpy.linspace(0, self.npoints-1, self.npoints)
            y   = numpy.linspace(0, n_indexes-1, n_indexes) * dy
            X,Y = numpy.meshgrid(x, y)
            Z   = numpy.zeros([n_indexes, self.npoints], dtype='d')
            for i in range(n_indexes):
                for j in range(self.npoints):
                    Z[i][j] = self.spectra[ index_list[i] ][1][j]

            if   (type=='wireframe'): ax.plot_wireframe(X, Y, Z, cstride=cstride, rstride=rstride)
            elif (type=='surface'):   ax.plot_surface(  X, Y, Z, cstride=cstride, rstride=rstride)

        plt.legend()            
        plt.show()

        return (status, message)



    def plot_spectra(self, logx=False, logy=False, xlabel='x', ylabel='y', kind='all', explabel=True, 
                     title=None, grid=True, line=PLOT_LINE, symbol=PLOT_SYMBOL):
        """ Plots a scatter plot of all the spectra together.

            Parameters
            ----------

            logx:           set logaritmic x axis
            logy:           set logaritmic y axis
            xlabel:         string to be used as label for the x axis
            ylabel:         string to be used as label for the y axis
            kind:           define which type to spectra plot
                            'all'       -> plot all
                            'spectrum'  -> only plot spectra (no histogram)
                            'histogram' -> only plot histograms
            explabel:       if True, use scientific format for labels if they execeed 100
            title:          string to be used as title
            grid:           wether to show a grid on the plot
            line:           type of line to be used for plotting
            symbol:         type of symbol to use for data point

            Returns
            -------

            status:  0 = no error
                     1 = no spectra have been loaded yet
                     2 = unkwnowkn kind
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK         

        if (len(self.spectra) == 0):
            status = 1
            message = "no spectra have been loaded yet"
            return (status, message)
        if (kind not in ('spectrum', 'histogram', 'all' )):
            status = 2
            message = 'unknow kind \"' + str(kind) + '\"'
            return (status, message)

        if (title is not None): plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(grid)
        if explabel: plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
        for i in range(len(self.spectra)):
            label = '(' + str(i).rjust(2) + ') ' + self.spectra[i][0]
            if (self.spectra[i][2].startswith('histogram') and ((kind in ('all', 'histogram')))):
                left   = self.xvalues - self.delta / 2.0
                height = self.spectra[i][1]
                width  = self.delta
                plt.bar(left=left, height=height, width=width, log=logy, label=label, fill=False) 
            elif (self.spectra[i][2].startswith('spectrum') and ((kind in ('all', 'spectrum')))):
                if (logx and logy):
                    plt.loglog(self.xvalues,   self.spectra[i][1], label=label, 
                              marker=symbol, linestyle=line, color=PLOT_COLORS[i%N_PLOT_COLORS])
                elif logx:
                    plt.semilogx(self.xvalues, self.spectra[i][1], label=label, 
                                 marker=symbol, linestyle=line, color=PLOT_COLORS[i%N_PLOT_COLORS])
                elif logy:
                    plt.semilogy(self.xvalues, self.spectra[i][1], label=label,
                                 marker=symbol, linestyle=line, color=PLOT_COLORS[i%N_PLOT_COLORS])
                else:
                    plt.plot(self.xvalues, self.spectra[i][1], label=label, 
                             marker=symbol, linestyle=line, color=PLOT_COLORS[i%N_PLOT_COLORS])
        plt.legend()            
        plt.show()

        return (status, message)


    def save_spectrum(self, index, filename, sep='\t', nan2zero=None):
        """ Saves a spectrum (range of y-values) to an ASCII file.

            Parameters
            ----------

            index:    index of the spectrum to be saved
            filename: name of the file to which the spectrum must be written
            sep:      string to be used as separator between x and y columns 
                      (default is a tab character)
            nan2zero: if not set (i.e. set to None), nan values will be removed before saving
                      if set to True, nan values will be saved as zeros
                      if set to False, nan values will be left in the file

            Returns
            -------

            status:  0 = no error
                     1 = could not open file
                     2 = table inconsistency
                     3 = error during file writing (this *should* never happen...)
                     4 = spectrum does not exist
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status  = 4
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        if (nan2zero is None):
            not_nan_indexes = numpy.isfinite(self.spectra[index][1])
            save_x = self.xvalues[not_nan_indexes]
            save_y = self.spectra[index][1][not_nan_indexes]
        elif (nan2zero == True):
            save_x = numpy.nan_to_num(self.xvalues)
            save_y = numpy.nan_to_num(self.spectra[index][1])
        else:
            save_x = self.xvalues
            save_y = self.spectra[index][1]
        save_data = data_manager.DataGrid()
        save_data.data_array = numpy.array([save_x, save_y])

        (status, message) = save_data.write_file(filename, sep=sep, transpose= True)

        return (status, message)
           
        
 #  +---------------------------------------------------------+
 #  ! Methods to create new spectra from data or calculations |
 #  +---------------------------------------------------------+

    def calculate_spectrum_x(self, name=None, function=None):
        """ Adds a spectrum calculating y-values as a given function of x-values.

            Parameters
            ----------

            name       string used to identify the spectrum
            function   function to be used to calculate y-values

            Returns
            -------

            status:  0 = no error
                     1 = no function was given
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        # Check that a function was given
        if (function is None):
            status = 1
            message = 'no function was given'
            return (status, message)
        # Check that a name was given and that it is a string
        if  (name is None):
            name = 'f(x)'
        else:
            name = str(name)
        self.spectra.append([name, numpy.zeros(self.npoints), 'spectrum'])
        spectrum_index = len(self.spectra)-1
        for i in range(self.npoints):
            self.spectra[spectrum_index][1][i] = function(self.xvalues[i])

        return (status, message)


    def calculate_spectrum_y(self, index, name=None, function=None):
        """ Adds a spectrum calculating y-values as a given function of y-values of another spectrum.

            Parameters
            ----------

            index:     index of the spectrum of which the y-values must be used as y variable
            name:      string used to identify the new spectrum
            function:  function(y) to be used to calculate the new y-values
                       y -> y-values of the spectrum identified by index                    

            Returns
            -------

            status:  0 = no error
                     1 = incorrect spectrum index
                     2 = no function was given
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)
        # Check that a function was given
        if (function is None):
            status = 2
            message = 'no function was given'
            return (status, message)
        # Check that a name was given and that it is a string
        if (name is None):
            name = 'f(' + self.spectra[index][0] + ')'
        else:
            name = str(name)

        self.spectra.append([name, numpy.zeros(self.npoints), 'spectrum'])
        new_spectrum_index = len(self.spectra)-1

        for i in range(self.npoints):
            self.spectra[new_spectrum_index][1][i] = function(self.spectra[index][1][i])

        return (status, message)


    def calculate_spectrum_xy(self, index, name=None, function=None):
        """ Adds a spectrum calculating new y-values as a given function of the x-values 
            and of the y-values of another spectrum.

            Parameters
            ----------

            index:     index of the spectrum of which the y-values must be used as y variable
            name:      string used to identify the new spectrum
            function:  function(x,y) to be used to calculate the new y-values
                       x -> x-values 
                       y -> y-values of the spectrum identified by index                    

            Returns
            -------

            status:  0 = no error
                     1 = incorrect spectrum index
                     2 = no function was given
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        # Check that a function was given
        if (function is None):
            status = 2
            message = 'no function was given'
            return (status, message)

        # Check that a name was given and that it is a string
        if (name is None):
            name = 'f(x,' + self.spectra[index][0] + ')'
        else:
            name = str(name)

        self.spectra.append([name, numpy.zeros(self.npoints), 'spectrum'])
        new_spectrum_index = len(self.spectra)-1

        for i in range(self.npoints):
            self.spectra[new_spectrum_index][1][i] = function(self.xvalues[i],
                                                              self.spectra[index][1][i])

        return (status, message)


    def calculate_spectrum_yy(self, index1, index2, name=None, function=None):
        """ Adds a spectrum calculating y-values as a given function of y-values of two other spectra.

            Parameters
            ----------

            index1,     
            index2:    indexes of the spectra used for calculation of the new y-values
            name:      string used to identify the new spectrum
            function:  function(y1,y2) to be used to calculate the new y-values
                       y1 -> y-values of the spectrum identified by index1  
                       y2 -> y-values of the spectrum identified by index2

            Returns
            -------

            status:  0 = no error
                     1 = incorrect spectrum index
                     2 = no function was given
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index1 = int(index1)
        index2 = int(index2)
        if ( (index1 < 0) or (index1 >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index1) + ' does not exist'
            return (status, message)
        if ( (index2 < 0) or (index2 >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index2) + ' does not exist'
            return (status, message)
        # Check that a function was given
        if (function is None):
            status = 2
            message = 'no function was given'
            return (status, message)
        # Check that a name was given and that it is a string
        if (name is None):
            name = 'f(' + self.spectra[index1][0] + ', ' + self.spectra[index2][0] + ')'
        else:
            name = str(name)

        self.spectra.append([name, numpy.zeros(self.npoints), 'spectrum'])
        new_spectrum_index = len(self.spectra)-1
        for i in range(self.npoints):
            self.spectra[new_spectrum_index][1][i] = function(self.spectra[index1][1][i],
                                                              self.spectra[index2][1][i])

        return (status, message)

            
    def add_spectrum(self, name=None, yvalues=None, left=None, right=None,
                     filename='', skip=0, sep='\t', plot=False):
        """ Adds a spectrum from a given range of y-values or reading data from a ASCII file.

            Parameters
            ----------

            name:     string used to identify the spectrum
            yvalues:  numpy array defining the y values
                      the number of values must be the same of self.xrange
                      if it is not given, it will be read from a file
            left:     value to be used as y-value for x-values which are 
                      lower then the minimum x value
            right:    value to be used as y-value for x-values which are 
                      greater then the maximum x value
            filename: name of the file from which the spectrum must be read
                      the file must contain a table of (x,y) values
                      the y values will be interpolated at the values in xrange
            skip:     number of lines to skip at the beginning of the file 
            sep:      character used to separate the x and y value in each line of the file
            plot:     if set to True, the spectrum will be plotted after being read from file

            Returns
            -------

            status:  0 = no error
                     1 = could not open file
                     2 = syntax error
                     3 = non-numerical value
                     4 = invalid number of rows/columns
                     5 = invalid number of rows to skip
                     6 = lenght of array yrange is wrong
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        # Check that if name was given
        if (name is None):
            name = 'Unnamed'
        else:
            name = str(name)
        # If an array is not given, try to read the spectrum from a file 
        if (yvalues is None):
            # Read the x and y values as a matrix from a file (two colums expected in the file)
            matrix = data_manager.DataGrid()
            (status, message) = matrix.read_file(filename, n_columns=2,
                                                 transpose=True, sep=sep, skip=skip)
            if (status != 0): return (status, message)
            # Order the (x,y) matrix using the x values as key (necessary for interpolation)
            order = numpy.argsort(matrix.data_array[0])
            # Calculate by interpolation the y-values corresponding to
            # the x-values stored in self.xvalues
            yvalues = numpy.interp(self.xvalues, matrix.data_array[0][order],
                                   matrix.data_array[1][order],
                                   left=left, right=right)
            if plot:
                plt.title(filename)
                plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
                l1 = 'read values'
                l2 = 'interpolated values'
                plt.plot(matrix.data_array[0], matrix.data_array[1], label = l1, 
                           marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'blue')
                plt.plot(self.xvalues, yvalues, label = l2, 
                           marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'red')
                plt.legend()
                plt.show()
            del matrix
            
        # If an array of values was given, check that it has the correct shape    
        else:                   
            if ( len(yvalues) != self.npoints ):
                status  = 6
                message = 'lenght of array yrange is wrong'
                return (status, message)
            if plot:
                plt.plot(self.xvalues, yvalues, PLOT_SYMBOL)
                plt.show()

        self.spectra.append([name, yvalues, 'spectrum'])

        del yvalues

        return (status, message)


    def add_histogram(self, name=None, values=None, filename='', htype='absolute', function=None):
        """ Adds a spectrum as a histogram.

            Adds a new histogram, that is a spectrum calculated as the frequencies
            of an array of values. These can be given as a parameter or loaded from a file. 
            The bins are defined in such a way that each x-value is in the center of a bin.

            self.xvalues[i]              -> i-th xvalue
            self.delta                   -> separation between xvalues
            self.xvalues[i]-self.delta/2 -> lower boundary of i-th bin
            self.xvalues[i]+self.delta/2 -> upper boundary of i-th bin

            Parameters
            ----------

            name:     string used to identify the new spectrum
            values:   array of values from which the frequencies will be calculated
            filename: name of the file from which the frequencies will be calcualted
                      if no array of values was given
                      the file must contain a sequence of values
            htype:    'absolute' -> stores absolute frequencies
                      'relative' -> stores relative frequencies (sum of all frequencies is 1)
                      'density'  -> stores probability densities (integral of the histogram is 1)
            function: if given, this function is applied to the values, before creating the histogram

            Returns
            -------

            status:  0 = no error
                     1 = could not open file
                     2 = syntax error
                     3 = non-numerical value
                     4 = invalid number of rows/columns
                     5 = invalid number of rows to skip
                     6 = lenght of array yrange is wrong
                     7 = unknown histogram type
            message: a string containing an error message or 'Ok'
        """

        status,message = 0, OK

        # Check if a name was given
        if (name is None):
            name = 'Unnamed'
        else:
            name = str(name)

        if (htype not in ('absolute', 'relative', 'density')):
            status = 7
            message = 'unknown histogram type'
            return (status, message)

        # If an array of values was not given, try to read values from a file   
        if (values is None):
            data = data_manager.DataSequence()
            (status, message) = data.read_file(filename)
            if (status != 0): return (status, message)
            values = data.data_array
            del data

        # If a function was given, apply it to the values before filling the histogram
        if (function is not None):
            for i in range(len(values)):
                values[i] = function(values[i])

        xmin = self.xvalues.min() - self.delta / 2.0
        xmax = self.xvalues.max() + self.delta / 2.0
        if (htype == 'density'):
            hist = numpy.histogram(values, bins=self.npoints, range=(xmin,xmax), density=True)      
            yvalues = hist[0]
        else:
            hist = numpy.histogram(values, bins=self.npoints, range=(xmin,xmax), density=False)
            if (htype == 'relative'):
                yvalues = hist[0] / float(numpy.sum(hist[0]))
            else:
                yvalues = hist[0]

        spectrum_type = 'histogram ' + htype            
        self.spectra.append([name, yvalues, spectrum_type])     

        return (status, message)

    
    def add_histogram_pdf(self, name=None,  htype='absolute', pdf=None, N=None, normalize=False):
        """Adds a spectrum as a histogram from a pdf (probability density function).

            Adds a new spectrum, which is an histogram calculated from a given pdf. 
            The bins are defined in such a way that each x-value is in the center of a bin.

            self.xvalues[i]              -> i-th xvalue
            self.delta                   -> separation between xvalues
            self.xvalues[i]-self.delta/2 -> lower boundary of i-th bin
            self.xvalues[i]+self.delta/2 -> upper boundary of i-th bin

            Parameters
            ----------

            name:      string used to identify the new spectrum
            htype:     'absolute' -> stores absolute frequencies
                       'relative' -> stores relative frequencies (sum of all frequencies is 1)
                       'density'  -> stores probability densities (integral of the histogram is 1)
            pdf:       pdf from which the histogram is calculated
            N:         number of values, needed to calculate absolute frequencies
            normalize: if set to True, force the sum (or integral for htype ='density') to be 1 
                       (or N, for htype='absolute') even if the range of xvalues does not cover the 
                       whole range over which the pdf is normalized

            Returns
            -------

            status:  0 = no error
                     1 = no pdf was given
                     2 = unknown histogram type
                     3 = number of values missing
            message: a string containing an error message or 'Ok'
        """

        status,message = 0, OK

        # Check if a name was given
        if (name is None):
            name = 'Unnamed'
        else:
            name = str(name)

        if (pdf is None):
            status = 1
            message = 'no pdf was given'
            return (status, message)

        if (htype not in ('absolute', 'relative', 'density')):
            status = 2
            message = 'unknown histogram type'
            return (status, message)

        if ( (htype == 'absolute') and (N is None) ):
            status = 3
            message = 'number of values missing'
            return (status, message)

        # Initialize to zero the frequencies
        yvalues = numpy.zeros(self.npoints)

        # Calculate frequency for each bin, integrating the pdf over the bin width, 
        # and then multiplying by the normalization factor, depending on the histogram type
        for i in range(self.npoints-1):
            left_margin  = self.xvalues[i] - self.delta / 2.0
            right_margin = self.xvalues[i] + self.delta / 2.0 
            (frequency_rel, err) = quad(pdf, left_margin, right_margin)
            if   (htype == 'absolute'): yvalues[i] = N * frequency_rel
            elif (htype == 'density'):  yvalues[i] = frequency_rel / self.delta
            else:                       yvalues[i] = frequency_rel

        # If required, force normalization
        if normalize:
            if   (htype == 'absolute'): yvalues = yvalues * N / yvalues.sum()
            elif (htype == 'density'):  yvalues = yvalues / (yvalues.sum() * self.delta)
            else:                       yvalues = yvalues / yvalues.sum()

        spectrum_type = 'histogram ' + htype    

        self.spectra.append([name, yvalues, spectrum_type])     

        return (status, message)


    def change_xvalues(self, function, plot=False):
        """ Calculate new x-values and changes y-values accordingly.

            Calculates new x-values on the basis of a given function
            the new x-values are then ordered and equiparted
            and y-values of all spectra are recalculated by linear interpolation

            WARNING. This operation cannot be applied on the histograms,
            since changing the x-values, also the bins are changed, and the frequencies
            lose their meaning. So histograms are removed from the DataSet.

            Parameters
            ----------

            function:     the function which will be used to calculate the new x-values
            plot:         if True, plots the new spectra 

            Returns
            -------

            status:  0 = no error
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        # Calculate the new x-values based on the given function
        # the new array will not not always be ordered nor equiparted, 
        # while self.xvalues must be ordered and equiparted,
        # so it will be necessary to modify it later
        new_xvalues = numpy.ones(self.npoints) * numpy.nan
        for i in range(self.npoints):
            new_xvalues[i] = function(self.xvalues[i])

        # Define the new x-values array, ordered and equiparted
        xmin = new_xvalues.min()
        xmax = new_xvalues.max()
        (self.xvalues, self.delta) = numpy.linspace(xmin, xmax, self.npoints, retstep = True)

        # For each spectrum, calculate by interpolation the y-values 
        # corresponding to the ordered and equiparted x-values
        # histograms will be removed, since changing the x-values we are changing the
        # bins and frequencies become meaningless
        order = numpy.argsort(new_xvalues)
        remove_list = []
        for i in range(len(self.spectra)):
            # If it is a histogram, add to the remove list
            if (self.spectra[i][2] != 'spectrum'):
                status = 1
                message = 'warning: some spectra were removed'
                remove_list.append(i)
                continue
            # If it is a spectrum, calculate the y-values at the new x-values
            yvalues = numpy.interp(self.xvalues, new_xvalues[order], self.spectra[i][1][order])
            if plot:
                if (i == 0):
                    l1 = 'original values'
                    l2 = 'interpolated values'
                    plt.plot(new_xvalues, self.spectra[i][1], label = l1, 
                             marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'blue')
                    plt.plot(self.xvalues, yvalues, label = l2,
                             marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'red')
                    plt.legend()
                else:   
                    plt.plot(new_xvalues, self.spectra[i][1], 
                             marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'blue')
                    plt.plot(self.xvalues, yvalues,
                             marker = PLOT_SYMBOL, linestyle = PLOT_LINE, color = 'red')

                plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2)) 
                plt.show()
            self.spectra[i][1] = yvalues
            
        # if the remove list is not empty, delete the histograms
        if (len(remove_list) > 0):
            # Delete the spectra in remove list, starting from the last one added, since
            # otherwise the indexes change after every deletion
            for i in range(len(remove_list), 0, -1):
                self.delete_spectrum(remove_list[i-1])

        return (status, message)
            


#   +-----------------------------------------------------------+
#   | Methods to calculate smoothed, averaged, unspiked spectra |
#   +-----------------------------------------------------------+
        
    def add_mean_line(self, index, name=None):
        """ Adds a spectrum as the mean of all the points of another spectrum.

            Calculates a new spectrum, in which all y-values are equal to the 
            mean of the original spectrum

            Parameters
            ----------

            index:   index of the spectrum of which the mean line is requested
            name:    string to be used as name of the mean line spectrum
                     if not given, a defalut name will be used

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        if (name is None):
            name = self.spectra[index][0] + ' mean'
        else:
            name = str(name)

        # Calculate the mean of the spectrum, excluding inf and NaNs from the calculation
        mean = self.spectra[index][1][ numpy.isfinite( self.spectra[index][1] ) ].mean()

        self.spectra.append([name, numpy.ones(self.npoints) * mean, 'spectrum mean line'])

        return (status, message)

            

    def add_smooth(self, index, name=None, method='uniform',
                   weights=None, window=1, times=1, pad=False,
                   k=1, dx=1.0, fortran=False):
        """ Adds the smoothed version of a spectrum.

            Parameters
            ----------

            index:   index of the spectrum of which the smoothed version is requested
            name:    string to be used as name of the smoothed spectrum
                     if not given, a defalut name will be used
            method:  method used to calculate the moving average: 
                     'uniform':  moving average without weights (i.e. uniform weights)
                     'weights':  moving average with weights given as a parameter
                     'linear':   moving average with linear weights
                     'gaussian': moving average with gaussian weights
                     'ISO':      moving average with gaussian weights, 
                                 with standard deviation calculated following 
                                 ISO11562-1996 for roughness measurament
                     'median':   median filter
                     'spikes':   median filter, applied only to points that are recognised as spikes
            weights: if method='weights' this parameter must be an array with the weights
            window:  the spectrum is smoothed by a moving average, or a median filter,
                     calculated over 2*window+1 points
            times:   used only with the 'uniform', 'linear', 'weights' and 'gaussian' methods:
                     number of times smoothing must be repeated
            pad:     used only with the 'median' and spikes' method: 
                     pad the window with the first/last point 
            k:       used only with the 'spikes' method: increasing k will become 
                     more difficult for a point to be considered a spike
            dx:      used only with the 'ISO' method: distance between points
            fortran: use fortran compiled module for calculation, if available

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = unknown method
                     3 = window value must be strictly positive
                     4 = window value is more or equal to half the number of points 
                         (except for 'median' and 'spikes' method with pad option active)
                     5 = number of times must be strictly positive
                     6 = k must be strictly positive
                     7 = wights method was selected, but weights wer noto given
                     8 = fortran option requested, but not available

            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)
        if (method not in ('weights', 'uniform', 'linear', 'gaussian', 'ISO', 'median', 'spikes')):
            status = 2
            message = 'unknown method'
            return (status, message)
        if (method == 'weights'):
            if (weights is None):
                status = 7
                message = 'weights were not given'
                return (status, message)
            else:
                window = len(weights) - 1
        if ( fortran and (method not in ('weights')) ):
            status = 8
            message = 'fortran option is not available for method \"' + method + '\"'
            return (status, message)
        window = int(window)
        if (window <= 0):
            status = 3
            message = 'window value must be strictly positive'
            return (status, message)
        elif (window*times >= self.npoints/2):
            if ( (method not in ('median', 'spikes')) or (pad==False) ):
                status = 4
                message = 'window value is greater or equal to half the number of points'
                return (status, message)                
        times = int(times)
        if (times <= 0):
            status = 5
            message = 'number of times must be strictly positive'
            return (status, message)
        if (k<=0):
            status = 6
            message = 'k must be greater than zero'
            return (status, message)
        if (name is None):
            name = self.spectra[index][0] + ' s' + method[0] + ' w' + str(int(window))
            if (times > 1): name = name + 'x' + str(times)
        else:
            name = str(name)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum smooth'])
        i_smooth = len(self.spectra)-1

        if (method == 'uniform'):
            self.spectra[i_smooth][1] = moving_average(self.spectra[index][1], window)
        elif (method == 'weights'):
            if (fortran):
                self.spectra[i_smooth][1] = f_math.moving_average_weighted(
                    self.spectra[index][1], weights, fill=numpy.nan)
            else:
                self.spectra[i_smooth][1] = moving_average_weighted(self.spectra[index][1], weights)
        elif (method in ('linear', 'gaussian', 'ISO')):
            if (method == 'ISO'):
                method  = 'gaussian'
                Lc      = 2 * window * dx  # This is the cutoff length
                sigma   = numpy.sqrt( numpy.log(2.0) / 2.0 ) * Lc / numpy.pi
                weights = calculate_weights(window+1, method=method, sigma=sigma)
            else:
                weights = calculate_weights(window+1, method=method)
            if (fortran):
                self.spectra[i_smooth][1] = f_math.moving_average_weighted(
                    self.spectra[index][1], weights, fill=numpy.nan)
            else:
                self.spectra[i_smooth][1] = moving_average_weighted(self.spectra[index][1], weights)
        else:                        
            self.spectra[i_smooth][1] = median_filter(self.spectra[index][1], window,
                                                      spikes=(method=='spikes'), pad=pad, k=k)
        if ((times > 1) and (method in ('uniform', 'linear', 'weights', 'gaussian'))):
            skip = 0
            for i in range(times-1):
                skip = skip + window
                if (method=='uniform'):
                    self.spectra[i_smooth][1] = moving_average_skip(
                        self.spectra[i_smooth][1], window, skip)
                elif (method in ('linear', 'weights', 'gaussian')):
                    if (fortan):
                        self.spectra[i_smooth][1] = f_math.moving_average_weighted(
                            self.spectra[i_smooth][1], weights, skip, fill=numpy.nan)
                    else:
                        self.spectra[i_smooth][1] = moving_average_weighted(
                            self.spectra[i_smooth][1], weights, skip)

        return (status, message)


    def add_derivative(self, index, name=None, window=0, times=1):
        """ Adds the derivative of a spectrum.

            Parameters
            ----------

            index:   index of the spectrum of which the derivative is requested
            name:    string to use as name of the derivative spectrum
                     if not given, a default value will be used
            window:  if non-zero, the derivative is smoothed by a moving average
                     calculated over 2*window+1 points
            times:   if smoothing is requested, it is repeated this number of times

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = window value must be strictly positive
                     3 = number of times must be strictly positive
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ((index < 0) or (index >= len(self.spectra))):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)        

        window = int(window)
        if (window < 0):
            status = 2
            message = 'window value must be positive'
            return (status, message)

        times = int(times)
        if (times <= 0):
            status = 3
            message = 'number of times must be strictly positive'
            return (status, message)

        if (name is None):
            name = self.spectra[index][0] + ' der'
            if (window > 0):
                name = name + ' w' + str(int(window))
            if (times > 1):
                name = name + 'x' + str(times)
        else:
            name = str(name)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum der'])
        i_der = len(self.spectra)-1

        yvalues = derivative_5p(self.delta, self.spectra[index][1])
        if (window > 0):
            self.spectra[i_der][1] = moving_average_skip(yvalues, window, 2)
            if (times > 1):
                skip = 2
                for i in range(times-1):
                    skip = skip + window
                    self.spectra[i_der][1] = moving_average_skip(
                        self.spectra[i_der][1], window, skip)      
        else:
            self.spectra[i_der][1] = yvalues                

        return (status, message)

    
    def add_min_max(self, index, name_min=None, name_max=None, window=0,
                    times=1, xmin=None, xmax=None, purge=True):
        """ Find minima and maxima of a spectrum.

            To find minima and maxima, a derivative is calculated and
            the intersections with the y=0 axis are searched.

            Values of minima and maxima are stored in numpy arrays
            as new spectra (items of self.spectra)

            Parameters
            ----------

            index:    index of the spectrum of which the minima and maxima are requested
            name_min, strings to be used for the names of the min and max spectra
            name_max: if not given, default values will be used
            window:   if given and non-zero, the derivative is smoothed by a moving average
                      calculated over 2*window+1 points
            times:    number of times which smoothing must be repeated (used only if window > 0)
            xmin,
            xmax:     if given, searches minima and maxima only in the window [xmin, xmax]
            purge:    if True, the spectrum of the derivative will be removed after being used
                      if False, will be left in the spectra set

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = window value must be positive
                     3 = number of times must be strictly positive
                     4 = the given window is outside the spectra
                     5 = xmin is greater or equal to xmax
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        window = int(window)
        if (window < 0):
            status = 2
            message = 'window value must be positive'
            return (status, message)
        
        times = int(times)
        if (times <= 0):
            status = 3
            message = 'number of times must be strictly positive'
            return (status, message)

        if (xmin is not None):
            if ((xmin < self.xvalues[0]) or (xmin > self.xvalues[self.npoints-1])):
                status = 4
                message = 'xmin is outside the x-values extremes'
                return (status, message)
            
        if (xmax is not None):
            if ((xmax < self.xvalues[0]) or (xmax > self.xvalues[self.npoints-1])):
                status = 4
                message = 'xmax is outside the x-values extremes'
                return (status, message)
            if ((xmin is not None) and (xmax <= xmin)):
                status = 5
                message = 'xmax must be greater than xmin'
                return (status, message)

        if ((name_min is None) or (name_max is None)):
            add_string = ''
            if (window > 0):
                add_string = add_string + ' w' + str(int(window))
                if (times > 1): 
                    add_string = add_string + 'x' + str(times)
            if ( (xmin is not None) or (xmax is not None) ):                
                if (xmin is None): str_min = str(self.xvalues[0])
                else:              str_min = str(xmin)
                if (xmax is None): str_max = str(self.xvalues[self.npoints-1])
                else:              str_max = str(xmax) 
                add_string = add_string + ' [' + str_min + ', ' + str_max + ']'

        if (name_min is None):  name_min = self.spectra[index][0] + ' min' + add_string
        else:                   name_min = str(name_min)
        if (name_max is None):  name_max = self.spectra[index][0] + ' Max' + add_string
        else:                   name_max = str(name_max)

        # Calculate the derivative of the spectrum      
        (status, message) = self.add_derivative(index, window = window, times = times)
        if (status != 0): return (status, message)
        i_der = len(self.spectra)-1

        # Create the numpy arrays into which minima and maxima will be stored
        self.spectra.append([name_min, numpy.ones(self.npoints) * numpy.nan, 'spectrum min'])
        self.spectra.append([name_max, numpy.ones(self.npoints) * numpy.nan, 'spectrum Max'])
        i_min = len(self.spectra)-2
        i_max = len(self.spectra)-1

        # Search for minima and maxima
        for i in range(1, self.npoints):
            # Check that the point is inside the given window, if a window was given
            if ( (xmin is not None) and (self.xvalues[i] < xmin) ): continue
            if ( (xmax is not None) and (self.xvalues[i] > xmax) ): break
            y0 = self.spectra[i_der][1][i-1]
            y1 = self.spectra[i_der][1][i]
            # If derivative passes from negative to positive, we have a minimum
            if ( (y0 < 0) and (y1 > 0) ):
                # Choose as index of the minimum the one for which
                # the derivative y-value is closer to zero
                if (abs(y0) < abs(y1)):
                    self.spectra[i_min][1][i-1] = self.spectra[index][1][i-1]
                else:
                    self.spectra[i_min][1][i] = self.spectra[index][1][i]
            # If derivative passes from positive to negative, we have a maximum
            if ( (y0 > 0) and (y1< 0) ):
                # Choose as index of the maximum the one for which
                # the derivative y-value is closer to zero
                if (abs(y0) < abs(y1)): 
                    self.spectra[i_max][1][i-1] = self.spectra[index][1][i-1]
                else:
                    self.spectra[i_max][1][i] = self.spectra[index][1][i]

        if purge: self.delete_spectrum(i_der)

        return (status, message)
    

    def add_interpolate_3p(self, index, name=None, xmin=None, xmax=None, mode='left'):
        """ Adds a new spectrum by interpolation from another one, using linear interpolation.

            Creates a new spectrum by interpolation from another one, 
            which has some (usually many) points missing, which are represented with nans.

            Parameters
            ----------

            index:   index of the spectrum of which the interpolated version is requested
            name:    string to be used as name of the new spectrum
                     if not given, a default value will be used
            xmin,
            xmax:    if given, the interpolated points are calculated only in the window [xmin, xmax]
                     other points will be stored as NaNs
            mode:    sets the mode by which the three interpolation points are changed
                     'left' - > shift when the middle point is passed
                     'right' -> shift when the right  point is passed

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = the given window is outside the spectra
                     3 = xmin is greater or equal to xmax
                     4 = less then 3 points are available for interpolation
                     5 = unknown mode
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ((index < 0) or (index >= len(self.spectra))):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)
        
        if ((xmin is not None) or (xmax is not None)):
            if ((xmin < self.xvalues[0]) or (xmax > self.xvalues[self.npoints-1])):
                status = 2
                message = 'xmin or xmax is outside the x-values extremes'
                return (status, message)
            if (xmin >= xmax):
                status = 3
                message = 'xmax must be greater than xmin'
                return (status, message)

        if (mode not in ('left', 'right')):
                status = 5
                message = 'unknown mode'
                return (status, message)

        if (name is None):
                name = self.spectra[index][0] + ' i3p'
                if ( (xmin is not None) or (xmax is not None) ):                
                        if (xmin is None): str_min = str(self.xvalues[0])
                        else:              str_min = str(xmin)
                        if (xmax is None): str_max = str(self.xvalues[self.npoints-1])
                        else:              str_max = str(xmax) 
                        name = name + ' [' + str_min + ', ' + str_max + ']'
        else:
                name = str(name)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum'])
        i_interp = len(self.spectra)-1

        # Search for the first three non-NaN points
        i0 = None
        i1 = None
        i2 = None
        for i in range(self.npoints):
            if not math.isnan(self.spectra[index][1][i]):
                if   (i0 is None): i0 = i
                elif (i1 is None): i1 = i
                elif (i2 is None):
                    i2 = i
                    break
        # If there are less than three non-NaN values, it is not possible to interpolate
        if (i2 is None):
            error = 4
            message = 'not enough points'
            self.spectra.pop(i_interp)
            return (error, message)
        x0 = self.xvalues[i0]
        x1 = self.xvalues[i1]
        x2 = self.xvalues[i2]
        y0 = self.spectra[index][1][i0]
        y1 = self.spectra[index][1][i1]
        y2 = self.spectra[index][1][i2]

        only_nan_remaining = False
        for k in range(self.npoints):
            # Decide when the three interpolation points must be shifted
            # 'left'  -> shift when the middle point is passed
            #            interpolation between first and middle point
            # 'right' -> shift when the right point is passed
            #            interpolation between middle and right point
            if (mode=='left'): ichange = i1
            else:              ichange = i2 

            # If the middle (or right) point is passed, shift the three points
            # if all the points on the right of i2 are Nan there shift is not possible
            # so we need to continue to the end with these points
            if ((k > ichange) and not only_nan_remaining):
                only_nan_remaining = True
                for i in range(i2+1, self.npoints):
                    if not math.isnan(self.spectra[index][1][i]):
                        i0 = i1
                        x0 = x1
                        y0 = y1
                        i1 = i2
                        x1 = x2
                        y1 = y2
                        i2 = i
                        x2 = self.xvalues[i2]
                        y2 = self.spectra[index][1][i2]
                        only_nan_remaining = False
                        break

            # Check that the point is inside the given xvalues window, if a window was given
            if ( (xmin is not None) and (self.xvalues[k] < xmin) ): continue
            if ( (xmax is not None) and (self.xvalues[k] > xmax) ): break

            # Calculate the value of the point by interpolation
            x = self.xvalues[k]             
            self.spectra[i_interp][1][k] = interpolate_3p(x, x0, x1, x2, y0, y1, y2)

        return (status, message)


    def add_spline(self, index, name=None, kind='cubic', xmin=None, xmax=None, extrapolate=False):
        """ Adds a new spectrum by interpolation of another one, using splines.

            Creates a new spectrum by interpolation (using splines), from anothe one
            which has some (usually many) points missing, which are represented by nans.

            Parameters
            ----------

            index: index of the spectrum of which the interpolated version is requested
            name:  string to be used as name of the new spectrum
                   if not given, a default value will be used
            kind:  kind of interpolation 
                   ('linear', 'nearest', 'zero', 'slinear', 'quadratic, 'cubic')
                   see documentation of scipy.interpolate.interp1d
            xmin,
            xmax:  if given, the interpolated points are calculated only in the window [xmin, xmax]
                   other point will be stored as NaNs

            extrapolate: if set to True, the extrapolation procedure is carried out also 

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = the given window is outside the spectra
                     3 = xmin is greater or equal to xmax
                     4 = too few points are available for interpolation
                     5 = unknown interpolation method
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        if ( (xmin is not None) or (xmax is not None) ):
             if ( (xmin < self.xvalues[0]) or (xmax > self.xvalues[self.npoints-1]) ) :
                 status = 2
                 message = 'xmin or xmax is outside the x-values extremes'
                 return (status, message)
             if ( (xmin >= xmax) ):
                 status = 3
                 message = 'xmax must be greater than xmin'
                 return (status, message)
             
        if (name is None):
            name = self.spectra[index][0] + ' spl'
            if ( (xmin is not None) or (xmax is not None) ):                
                if (xmin is None): str_min = str(self.xvalues[0])
                else:              str_min = str(xmin)
                if (xmax is None): str_max = str(self.xvalues[self.npoints-1])
                else:              str_max = str(xmax) 
                name = name + ' [' + str_min + ', ' + str_max + ']'
        else:
            name = str(name)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum spline'])
        i_interp = len(self.spectra)-1

        # Search for the non-NaN points and store them in two lists
        list_x = []
        list_y = []
        for i in range(self.npoints):
            if not math.isnan(self.spectra[index][1][i]):
                list_x.append(self.xvalues[i])
                list_y.append(self.spectra[index][1][i])

        # Check if there are enough non-NaN values to interpolate
        if   (kind=='cubic'):      N_min_points = 4
        elif (kind=='quadratic'):  N_min_points = 3
        else:                      N_min_points = 2
        if (len(list_x) < N_min_points):
            error = 4
            message = 'not enough points'
            self.spectra.pop(i_interp)
            return (error, message)

        # Define interpolation function (if required, set extrapolation option)
        if extrapolate: fill_value = 'extrapolate'
        else:           fill_value = numpy.nan
        try:
            int_func = interp1d(list_x, list_y, kind=kind, fill_value=fill_value)
        except NotImplementedError:
            error = 5
            message = 'unknown interpolation method \"' + kind + '\"'
            self.spectra.pop(i_interp)
            return (error, message)

        for k in range(self.npoints):
            x = self.xvalues[k]
            if ((x < list_x[0]  and not extrapolate) or ((xmin is not None) and (x < xmin))):
                continue
            if ((x > list_x[-1] and not extrapolate) or ((xmax is not None) and (x > xmax))):
                break
            # Calculate the value of the point by interpolation
            self.spectra[i_interp][1][k] = int_func(x)

        return (status, message)

    
    def add_average(self, index1, index2, name=None):
        """ Adds a new spectrum as the geometric average of other two.

            Parameters
            ----------

            index1,
            index2:  indexes of the two spectra used for the average
            name:    string used as the name of the interpolated spectrum
                     if not given, a default value is used

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
            message: a string containing an error message or 'Ok'           
        """

        status, message = 0, OK

        index1 = int(index1)
        index2 = int(index2)
        if ( (index1 < 0) or (index1 >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index1) + ' does not exist'
            return (status, message)

        if ( (index2 < 0) or (index2 >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index2) + ' does not exist'
            return (status, message)                

        if (name is None):
            name = self.spectra[index1][0] + ', ' + self.spectra[index2][0] + ' ave'
        else:
            name = str(name)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum average'])
        i_average = len(self.spectra)-1
        for i in range(self.npoints):
            self.spectra[i_average][1][i] = numpy.sqrt(
                self.spectra[index1][1][i] * self.spectra[index2][1][i] )

        return (status, message)

        
    def get_index_xvalue(self, x, check=False):
        """ Get the index corresponding to the nearest xvalue.

            Parameters
            ----------

            x:     value to be searched among xvalues
            check: if set to True,  return None if x is outside the xvalues range
                   if set to False, return the minimum or maximum of xvalues instead

            Returns
            -------

            index of the x-value nearest to the given value
            on error None is returned
        """

        status, message = 0, OK

        if check:
            if ((x < self.xvalues.min()) or (x > self.xvalues.max())): 
                return None
        for i in range(self.npoints):
            if (self.xvalues[i] >= x): break

        return i


    def change_yvalue(self, index, xvalue, yvalue=None, index2=None):
        """ Change the yvalue corresponding to a given xvalue, 
            giving the new yvalue or taking it from another spectrum.

            Parameters
            ----------

            index:  index of the spectrum where the yvalue must be changed
            xvalue: position on the x axis where the y value must be changed
            yvalue: value to be stored
            index2: index of the spectrum from which the value
                    to be stored must be taken
                    to be provided if yvalue is not given

            Returns
            -------
            status:  0 = no error
                     1 = wrong spectrum index
                     2 = xvalue is outside allowed range
                     3 = missing yvalue or spectrum index
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index = int(index)
        if ((index < 0) or (index >= len(self.spectra))):
            status = 1
            message = 'spectrum ' + str(index) + ' does not exist'
            return (status, message)

        ix = self.get_index_xvalue(xvalue, check=True)
        if (ix is None):
            status = 2
            message = 'xvalue ' + str(xvalue) + ' is outside allowed range'
            return (status, message)

        if (yvalue is None):
            if (index2 is None):
                status = 3
                message = 'missing yvalue or spectrum index'
                return (status, message)
            else:
                index2 = int(index2)
                if ( (index2 < 0) or (index2 >= len(self.spectra)) ):
                    status = 1
                    message = 'spectrum ' + str(index2) + ' does not exist'
                    return (status, message)
                yvalue = self.spectra[index2][1][ix]

        self.spectra[index][1][ix] = yvalue

        return (status, message)


    def get_first_last(self, index, xmin=None, xmax=None, kind='first'):
        """ Returns the firts (or last) non-NaN value of a spectrum.

            Parameters
            ---------

            index: index of the spectrum of which the first/last non-NaN value is requested
            xmin,
            xmax:  window inside which the search must be carried out
            kind:  'first' -> search the first non-NaN value
                   'last'  -> search the last non-NaN value

            Returns
            -------

            (x,y): x-value and y-value of the first or last non-NaN value of the spectrum
                   on error None is returned
        """

        status, message = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ): return None
        if ( (xmin is not None) and (xmax is not None) and (xmin >= xmax)): return None
        if ( (xmin is None) or (xmin < self.xvalues.min()) ): xmin = self.xvalues.min()
        if ( (xmax is None) or (xmax > self.xvalues.max()) ): xmax = self.xvalues.max()

        x = None
        y = None
        if (kind == 'first'):
            for i in range(self.npoints):
                if (self.xvalues[i] < xmin): continue
                if (self.xvalues[i] > xmax): break
                if not math.isnan(self.spectra[index][1][i]):
                   x = self.xvalues[i]
                   y = self.spectra[index][1][i]
                   break
            return (x,y)
        elif (kind == 'last'):
            for k in range(self.npoints, 0, -1):
                i = k - 1
                if (self.xvalues[i] < xmin): break
                if (self.xvalues[i] > xmax): continue
                if not math.isnan(self.spectra[index][1][i]):
                   x = self.xvalues[i]
                   y = self.spectra[index][1][i]
                   break
            return (x,y)
        else:
            return None
                

    def integrate_spectrum(self, index, xmin=None, xmax=None):
        """ Calculates the integral of a spectrum.

            Parameters
            ----------

            index: index identifying the spectrum to integrate
            xmin,
            xmax:  extremes of the integral

            Returns
            -------

            integral: the value of the integral
                      on error is returned None
        """

        index = int(index)
        if ((index < 0) or (index >= len(self.spectra)) ): return None
        if ((xmin is not None) and (xmax is not None) and (xmin >= xmax)):
            status = 2
            message = 'xmax must be greater than xmin'
            return (status, message)
        if ((xmin is None) or (xmin < self.xvalues.min()) ): xmin = self.xvalues.min()
        if ((xmax is None) or (xmax > self.xvalues.max()) ): xmax = self.xvalues.max()

        sum = 0
        for i in range(self.npoints):
            if (self.xvalues[i] < xmin): continue
            if (self.xvalues[i] > xmax): break
            sum = sum + self.spectra[index][1][i]     

        return sum * self.delta     

    
    def sum_spectrum(self, index, xmin=None, xmax=None):
        """ Sums the y-values of a spectrum.

            Parameters
            ----------

            index: index identifying the spectrum to sum
            xmin,
            xmax:  extremes of the sum

            Returns
            -------

            sum: the value of the sum
                 on error is returned None
        """

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ): return None
        if ( (xmin is not None) and (xmax is not None) and (xmin >= xmax)): return None
        if ( (xmin is None) or (xmin < self.xvalues.min()) ): xmin = self.xvalues.min()
        if ( (xmax is None) or (xmax > self.xvalues.max()) ): xmax = self.xvalues.max()

        sum = 0
        for i in range(self.npoints):
            if (self.xvalues[i] < xmin): continue
            if (self.xvalues[i] > xmax): break
            sum = sum + self.spectra[index][1][i]     

        return sum


    def mean_spectrum(self, index, function=None, xmin=None, xmax=None, purge=True):
        """ Calculates the mean value of a function weighted with a distribution, 
            stored in a histogram.

            Parameters
            ----------

            index:    index identifying the spectrum of which the mean is requested
                      it *must* be of type 'histogram relative' or 'histogram density'
                      otherwise None will be returned
            function: function of which the mean is requested
                      if not given, the identity f(x)=x is used
            xmin,
            xmax:     extremes of the range inside which the mean is requested
            purge:    if True, the spectra created for intermediate calculations 
                      will be removed after being used
                      if False, they will be left in the spectra set

            Returns
            -------

            mean: the value of the mean
                  on error, None is returned
        """

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ): return None
        if (self.spectra[index][2] not in ('histogram relative','histogram density' )): return None

        if ( (xmin is None) or (xmin < self.xvalues.min()) ): xmin = self.xvalues.min()
        if ( (xmax is None) or (xmax > self.xvalues.max()) ): xmax = self.xvalues.max()

        if (function is None): function = lambda x: x

        def new_function(x, y): return function(x)*y
        self.calculate_spectrum_xy(index=index, function=new_function)          
        if   (self.spectra[index][2] == 'histogram density' ):
            mean = self.integrate_spectrum(index=index+1, xmin=xmin, xmax=xmax)
        elif (self.spectra[index][2] == 'histogram relative'):
            mean = self.sum_spectrum(index=index+1, xmin=xmin, xmax=xmax)
        else:
            mean = None

        if purge: self.delete_spectrum(index+1)

        return mean 


#   +-------------------------------------------+
#   | Specific methods for optical applications |
#   +-------------------------------------------+
        

    def add_refractive_index(self, s, index_max=None, index_min=None,
                             index_average=None, name='n', method='transp1'):
        """ Adds a spectrum, as the refractive index calculated from an optical transmission spectrum.

            Calculates the refraction index of a thin film, from am optical transmission spectrum:
            depending on which method is chosen, different information is required

            'transp1' -> substrate refractive index and envelope of transmission minima
            'transp2' -> substrate refractive index and transmission average 
                         (gemetric average between maxima and minima)
            'lowabs'  -> substrate refractive index and envelopes of transmission maxima and minima


            Parameters
            ----------

            s:              refractive index of the substrate
            index_max:      index of the spectrum containg the envelope of maximum transmission
            index_min:      index of the spectrum containg the envelope of minimum transmission
            index_average:  index of the spectrum containg the envelope of average transmission
            name:           name to give to the spectrum of refractive index
            method:         calculation method

            Returns
            -------
            status:  0 = no error
                     1 = wrong spectrum index
                     2 = unknown calculation method
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        if (index_min is not None):
            index_min = int(index_min)
            if ( (index_min < 0) or (index_min >= len(self.spectra)) ):
                status = 1
                message = 'spectrum ' + str(index_min) + ' does not exist'
                return (status, message)
            
        if (index_max is not None):
            index_max = int(index_max)
            if ( (index_max < 0) or (index_max >= len(self.spectra)) ):
                status = 1
                message = 'spectrum ' + str(index_max) + ' does not exist'
                return (status, message)
            
        if (index_average is not None):
            index_average = int(index_average)
            if ( (index_average < 0) or (index_average >= len(self.spectra)) ):
                status = 1
                message = 'spectrum ' + str(index_average) + ' does not exist'
                return (status, message)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum n'])
        i_n = len(self.spectra)-1

        if (method == 'transp1'):
            if (index_min is None):
                status = 3
                message = 'method \"transp1\" needs index of the spectrum of minimum transmission'
                self.spectra.pop(i_n)
                return (status, message)
            for i in range(self.npoints):
                self.spectra[i_n][1][i] = n_transparent_1(s, self.spectra[index_min][1][i])
        elif (method == 'transp2'):
            if (index_average is None):
                status = 3
                message = 'method \"transp2\" needs index of the spectrum of average transmission'
                self.spectra.pop(i_n)
                return (status, message)
            for i in range(self.npoints):
                self.spectra[i_n][1][i] = n_transparent_2(s, self.spectra[index_average][1][i])
        elif (method == 'lowabs'):
            if ( (index_min is None) or (index_max is None) ):
                status = 3
                message = ('method \"lowabs\" needs index of the spectra '
                           + 'of minimum and maximum transmission')
                self.spectra.pop(i_n)
                return (status, message)
            for i in range(self.npoints):
                self.spectra[i_n][1][i] = n_lowabsorption(
                    s, self.spectra[index_max][1][i], self.spectra[index_min][1][i])
        else:
            status = 2
            message = 'unkwnown calculation method'
            self.spectra.pop(i_n)
            return (status, message)

        return (status, message)


    def add_n_brutal(self, thickness, index_peaks=None, name='nb'):
        """ Adds a spectrum, as the refractive index calculated 
            from maxima and/or minima of an optical transmission spectrum.

            Calculates the refraction index of a thin film, 
            from the optical transmission or reflection spectrum
            using the positions of maxima and/or minima

            Parameters
            ----------

            thickness:      thickness of the film in nanometers
            index_peaks:    index of the spectrum containing the maxima 
                            or minima of optical Transmission or Reflection
                            the x-values MUST be expressed in electronvolts
            name:           name to give to the spectrum of refractive index

            Returns
            -------
            status:  0 = no error
                     1 = wrong spectrum index
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        if (index_peaks is not None):
            index_peaks = int(index_peaks)
            if ( (index_peaks < 0) or (index_peaks >= len(self.spectra)) ):
                status = 1
                message = 'spectrum ' + str(index_peaks) + ' does not exist'
                return (status, message)

        if (index_peaks is None):
            status = 3
            message = 'index of the spectrum of minima/maxima of transmission/reflection is needed'
            return (status, message)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum n'])
        index_n = len(self.spectra)-1

        # Search first peak
        for i in range(self.npoints):
            if not math.isnan(self.spectra[index_peaks][1][i]): break
        i0 = i
        E0 = self.xvalues[i0]

        # Search following peaks and calculate n using ajacent peaks
        for i in range(i0 + 1, self.npoints):
            # if there is no peak (nan value) proceed to the next point
            if math.isnan(self.spectra[index_peaks][1][i]): continue
            # if there is a peak (non-nan value) use this x-value and the previous to calculate n
            i1 = i
            E1 = self.xvalues[i1]
            dE = E1 - E0
            ii = int( round( (i1 + i0) / 2 ) ) # index of the point halfway between the peaks
            n = n_brutal(dE, thickness)
            self.spectra[index_n][1][ii] = n
            # store this x-value for next calulation of n
            i0 = i1
            E0 = E1

        return (status, message)


    def add_absorption_coefficient(self, s, thickness, index_n, index_Ta, name='alpha', method='A3'):
        """ Adds a spcetrum representing the absorption coefficient of a thin film, 
            calculated from an optical transmission spectrum.

            Calculates the absorption coefficient of a thin film, from the Transmission spectrum
            information required is the same for each of the two methods implemented

            Parameters
            ----------

            s:              refractive index of the substrate
            thickenss:      thickness of the film
            index_n:        index of the spectrum containing the refractive index of the film
            index_Ta:       index of the spectrum containg the average (fringes free) 
                            transmission spectrum
            name:           name to give to the spectrum of absorption coefficient
            method:         calculation method

            Returns
            -------

            status:  0 = no error
                     1 = wrong spectrum index
                     2 = unknown calculation method
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK

        index_n = int(index_n)
        if ( (index_n < 0) or (index_n >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index_n) + ' does not exist'
            return (status, message)
        
        index_Ta = int(index_Ta)
        if ( (index_Ta < 0) or (index_Ta >= len(self.spectra)) ):
            status = 1
            message = 'spectrum ' + str(index_Ta) + ' does not exist'
            return (status, message)

        self.spectra.append([name, numpy.ones(self.npoints) * numpy.nan, 'spectrum alpha'])
        i_alpha = len(self.spectra)-1

        if (method == 'A3'):
            for i in range(self.npoints):
                self.spectra[i_alpha][1][i] = alpha_A3(s,
                                                       self.spectra[index_n][1][i],
                                                       self.spectra[index_Ta][1][i],
                                                       thickness)
        elif (method == '19'):
            for i in range(self.npoints):
                self.spectra[i_alpha][1][i] = alpha_19(s,
                                                       self.spectra[index_n][1][i],
                                                       self.spectra[index_Ta][1][i],
                                                       thickness)
        else:
            status = 2
            message = 'unkwnown calculation method'
            self.spectra.pop(i_alpha)
            return (status, message)

        return (status, message)


#   +-----------------------------------------+
#   | Specific methods for roughness analysis |
#   +-----------------------------------------+
        

    def add_roughness(self, index, L, method='ISO', dx=1.0, purge=False, output=False):
        """ Extracts waviness and roughness from a given signal, 
            then calculates some roughness parameters.

            Calculates waviness and roughness signals by means of 
            the separation of long- and short-wave components of a given signal, 
            then calculates some roughness parameters: average height, Ra, Rq, Rsk, Rku.

            Parameters
            ---------

            index:    index of the spectrum from which roughness and waviness must be extracted
            L:        cutoff lenght to be used to separate waviness from roughness
            method:   method to use to separate waviness from roughness: 'gaussian', 'ISO', 'weights'
            dx:       used only with the 'ISO' method: distance between points
            purge:    if set to True, remove roughness and waviness spectra after calculation
            output:   if set to True, print

            Initialized data attributes
            ---------------------------

            self.error:  (status, message) used for error checking                               

            Returns
            -------

            (y_mean, Ra, Rq, Rsk, Rku)
            on error, None is return
        """

        self.error = 0, OK

        index = int(index)
        if ( (index < 0) or (index >= len(self.spectra)) ): return None

        # Create waviness spectrum
        w  = int( numpy.rint( L / (2.0*dx) ) )
        name = 'waviness('+str(index)+') L='+str(L)
        result = self.add_smooth(index, name=name, method=method, window=w, dx=dx)
        if (result[0] > 0):
            self.error = result
            return None
        index_waviness = len(self.spectra)-1

        # Create roughness spectrum, substracting waviness from original signal
        name = 'roughness('+str(index)+') L='+str(L)                
        result = self.calculate_spectrum_yy(index,
                                            index_waviness,
                                            name=name,
                                            function=lambda y0,y1: y0-y1)
        if (result[0] > 0):
            self.error = result
            return None
        index_roughness = len(self.spectra)-1

        # Calculate roughness parameters
        R = roughness(self.spectra[index_roughness][1])

        if output:
            print('\n')
            print('Parameter        '.ljust(20), 'Symbol'.ljust(10), 'Value')
            print('-----------------'.ljust(20), '------'.ljust(10), '-----')
            print('Cutoff lenght    '.ljust(20), 'Lc    '.ljust(10), L)
            print('Mean height      '.ljust(20), 'y_mean'.ljust(10), R[0])
            print('Roughness        '.ljust(20), 'Ra    '.ljust(10), R[1])
            print('Roughness        '.ljust(20), 'Rq    '.ljust(10), R[2])
            print('Skewness         '.ljust(20), 'Rsk   '.ljust(10), R[3])
            print('Kurtosis         '.ljust(20), 'Rku   '.ljust(10), R[4])                        

        if purge:
            self.delete_spectrum(index_roughness)                        
            self.delete_spectrum(index_waviness)                             

        return R
        

    def calculate_roughnesses(self, index_list, L_list, method='ISO', dx=1.0,
                              purge=True, output=True, uncert='hmd',
                              fortran=True):
        """ Extracts roughness and waviness spectra from a set of signals.

            Calculates waviness and roughness signals from a set of given signals.

            Parameters
            ---------

            index_List: list of indexes of the spectra from which 
                        roughness and waviness must be extracted
            L_list:     list of cutoff lenghts to be used to separate waviness from roughness
            method:     method to use to separate waviness from roughness: 'gaussian', 'ISO'
            dx:         used only in the 'ISO' method: distance between points
            purge:      if set to True, remove roughness and waviness spectra after calculation
            output:     if set to True, print mean values and uncertainty 
                        (as specified by the 'uncert' parameter) 
                        for all cutoff lenghts
            uncert:     select which value will be printed as uncertainty
                        'hmd':   half of the maximum deviation
                        'stdv':  standard deviation
                        'stdvm': standard deviation of the mean
            fortran:    use precompiled fortran module, if available

            Initialized data attributes
            ---------------------------

            self.Lc:       (numpy array) cutoff lenghts used for the separation 
                           of waviness and roughness
            self.ym_mean:  (numpy array) mean of the values of the mean heights of each spectrum
            self.ym_hmd:   (numpy array) half of maximum deviation of the  values of the mean height
            self.ym_stdv:  (numpy array) standard deviation of the values of the mean height
            self.ym_stdvm: (numpy array) standard deviation of the mean of 
                           the values of the mean height
            self.R_total   list of which each item is a numpy array containing 
                           the joined roughness spectra 
                           for a cutoff lenght

            Returns
            -------
            (status, message) used for error checking                                          
            status: 0 = no error
                    1 = a single index was given, instead of a list
                    2 = less than two spectra indexes were given
                    3 = a single value of Lc was given, instead of a list
                    4 = wrong spectrum index
                    5 = one of the cutoff lenghts is not a number
                    6 = unknown uncertainty type
                    7 = error during roughness calculation
            message: an error message or 'OK'
        """

        status, message = 0, OK

        try:
            n_spectra = len(index_list)
        except TypeError:
            status  = 1
            message = 'a list of indexes must be given'
            return (status, message)
        
        if (n_spectra < 2):
            status  = 2
            message = 'at least two indexes must be given'
            return (status, message)
        
        try:
            n_L = len(L_list)
        except TypeError:
            status  = 3
            message = 'a list of cutoff lenghts was not given'
            return (status, message)
        
        for i in range(n_spectra):
            index_list[i] = int(index_list[i])
            if ( (index_list[i] < 0) or (index_list[i] >= len(self.spectra)) ):
                status = 4
                message = 'spectrum ' + str(index_list[i]) + ' does not exist'
                return (status, message)
            
        for i in range(n_L):
            try:
                L_list[i] = abs(float(L_list[i]))
            except ValueError:
                status = 5
                message = 'Lc value \'' + str(L_list[i]) + '\' is not a valid number'
                return (status, message)
        if (uncert not in ('hmd', 'stdv', 'stdvm')):
            status = 6
            message = 'unknonw uncertainty type'
            return (status, message)

        self.Lc        = numpy.ones(n_L,'d') * numpy.nan

        self.ym_mean   = numpy.ones(n_L,'d') * numpy.nan
        self.ym_hmd    = numpy.ones(n_L,'d') * numpy.nan
        self.ym_stdv   = numpy.ones(n_L,'d') * numpy.nan
        self.ym_stdvm  = numpy.ones(n_L,'d') * numpy.nan

        self.Ra_mean   = numpy.ones(n_L,'d') * numpy.nan
        self.Ra_hmd    = numpy.ones(n_L,'d') * numpy.nan
        self.Ra_stdv   = numpy.ones(n_L,'d') * numpy.nan
        self.Ra_stdvm  = numpy.ones(n_L,'d') * numpy.nan

        self.Rq_mean   = numpy.ones(n_L,'d') * numpy.nan
        self.Rq_hmd    = numpy.ones(n_L,'d') * numpy.nan
        self.Rq_stdv   = numpy.ones(n_L,'d') * numpy.nan
        self.Rq_stdvm  = numpy.ones(n_L,'d') * numpy.nan

        self.Rsk_mean  = numpy.ones(n_L,'d') * numpy.nan
        self.Rsk_hmd   = numpy.ones(n_L,'d') * numpy.nan
        self.Rsk_stdv  = numpy.ones(n_L,'d') * numpy.nan
        self.Rsk_stdvm = numpy.ones(n_L,'d') * numpy.nan                

        self.Rku_mean  = numpy.ones(n_L,'d') * numpy.nan
        self.Rku_hmd   = numpy.ones(n_L,'d') * numpy.nan
        self.Rku_stdv  = numpy.ones(n_L,'d') * numpy.nan
        self.Rku_stdvm = numpy.ones(n_L,'d') * numpy.nan                

        # List of arrays to host the joined roughness signals for each Lc
        self.npoints_total = self.npoints * n_spectra
        self.R_total       = []

        self.Szm = numpy.ones(n_L,'d') * numpy.nan
        self.Sa  = numpy.ones(n_L,'d') * numpy.nan
        self.Sq  = numpy.ones(n_L,'d') * numpy.nan
        self.Ssk = numpy.ones(n_L,'d') * numpy.nan
        self.Sku = numpy.ones(n_L,'d') * numpy.nan

        # Iterate over the given Lc values
        for i in range(n_L):
            self.Lc[i] = L_list[i]

            # Prepare arrays to store the roughness parameters of each signal
            ym_array   = numpy.zeros([n_spectra], 'd')
            Ra_array   = numpy.zeros([n_spectra], 'd')
            Rq_array   = numpy.zeros([n_spectra], 'd')
            Rsk_array  = numpy.zeros([n_spectra], 'd')
            Rku_array  = numpy.zeros([n_spectra], 'd')

            # Prepare a numpy array to store the joined signals for this Lc
            self.R_total.append( numpy.ones(self.npoints_total,'d') * numpy.nan )

            # Index running over the joined signals array
            l = 0

            # Calculate the gaussian weights for this cutoff length
            w  = int( numpy.rint( L_list[i] / (2.0*dx) ) )
            if (method=='ISO'):
                sigma   = numpy.sqrt( numpy.log(2.0) / 2.0 ) * L_list[i] / numpy.pi
            else:
                sigma = None
            weights = calculate_weights(n=w+1, method='gaussian', sigma=sigma)
            if (output):
                print ('\n*** Lc = ' + str(L_list[i]) + ' (w = ' +str(w)
                       + '; sigma = ' + str(sigma) + ')')
                print ('Calculating track number: ',)

            # Iterate over the given signals
            for j in range(n_spectra):
                if (output): print(str(j+1)+'..',)
                # Create waviness spectrum, smoothing the original one with gaussian filter
                name = 'waviness(' + str( index_list[j] ) + ') L=' + str( L_list[i] )
                result = self.add_smooth(index_list[j],
                                         name=name,
                                         method='weights',
                                         weights=weights,
                                         fortran=fortran)
                if (result[0] > 0):
                    print ('\nERROR while processing index='
                           + str(index_list[j]) + ' Lc=' + str(L_list[i]))
                    return result
                index_waviness = len(self.spectra)-1

                # Create roughness spectrum, substracting waviness from original signal
                name = 'roughness(' + str( index_list[j] ) + ') L=' + str( L_list[i] )                
                result = self.calculate_spectrum_yy(index_list[j],
                                                    index_waviness,
                                                    name=name,
                                                    function=lambda y0,y1: y0-y1)
                if (result[0] > 0):
                    print ('\nERROR while processing index='
                           + str(index_list[j]) + ' Lc=' + str(L_list[i]))
                    return result
                index_roughness = len(self.spectra)-1

                # Calculate roughness parameters
                R = roughness(self.spectra[index_roughness][1])                                
                if (R is None):
                    print ('\nERROR while processing index='
                           + str(index_list[j]) + ' Lc=' + str(L_list[i]))
                    error = 7
                    message = 'error during roughness calculation'
                    return (error, message)
                ym_array[j]  = R[0]
                Ra_array[j]  = R[1]
                Rq_array[j]  = R[2]
                Rsk_array[j] = R[3]
                Rku_array[j] = R[4]

                # Createe joined spectrum
                for k in range(self.npoints):
                    self.R_total[i][l]=self.spectra[index_roughness][1][k]
                    l = l + 1

                # If purge requested, remove roughness and waviness spectra for this cutoff length
                if (purge):
                    self.delete_spectrum(index_roughness)
                    self.delete_spectrum(index_roughness-1)

            # Calculate areal parameters
            R = roughness(self.R_total[i])
            if (R is None):
                print('\nERROR while processing joined signal Lc=' + str(L_list[i]))
                error = 7
                message = 'error during roughness calculation'
                return (error, message)
            self.Szm[i] = R[0]                        
            self.Sa[i]  = R[1]
            self.Sq[i]  = R[2]
            self.Ssk[i] = R[3]
            self.Sku[i] = R[4]

            # Calculate mean values of linear parameters         
            self.ym_mean[i]  = numpy.average(ym_array)
            self.ym_hmd[i]   = ( ym_array.max() - ym_array.min() ) / 2.0
            self.ym_stdv[i]  = numpy.std(ym_array, ddof=1)
            self.ym_stdvm[i] = self.ym_stdv[i] / numpy.sqrt(n_spectra)

            self.Ra_mean[i]  = numpy.average(Ra_array)
            self.Ra_hmd[i]   = ( Ra_array.max() - Ra_array.min() ) / 2.0
            self.Ra_stdv[i]  = numpy.std(Ra_array, ddof=1)
            self.Ra_stdvm[i] = self.Ra_stdv[i] / numpy.sqrt(n_spectra)                        

            self.Rq_mean[i]  = numpy.average(Rq_array)
            self.Rq_hmd[i]   = ( Rq_array.max() - Rq_array.min() ) / 2.0
            self.Rq_stdv[i]  = numpy.std(Rq_array, ddof=1)
            self.Rq_stdvm[i] = self.Rq_stdv[i] / numpy.sqrt(n_spectra)

            self.Rsk_mean[i]  = numpy.average(Rsk_array)
            self.Rsk_hmd[i]   = ( Rsk_array.max() - Rsk_array.min() ) / 2.0
            self.Rsk_stdv[i]  = numpy.std(Rsk_array, ddof=1)
            self.Rsk_stdvm[i] = self.Rsk_stdv[i] / numpy.sqrt(n_spectra)

            self.Rku_mean[i]  = numpy.average(Rku_array)
            self.Rku_hmd[i]   = ( Rku_array.max() - Rku_array.min() ) / 2.0
            self.Rku_stdv[i]  = numpy.std(Rku_array, ddof=1)
            self.Rku_stdvm[i] = self.Rku_stdv[i] / numpy.sqrt(n_spectra)                       

            if output:
                if (uncert=='stdv'):
                    ym_str  = unit_manager.print_uncertainty(self.ym_mean[i],  self.ym_stdv[i])
                    Ra_str  = unit_manager.print_uncertainty(self.Ra_mean[i],  self.Ra_stdv[i])
                    Rq_str  = unit_manager.print_uncertainty(self.Rq_mean[i],  self.Rq_stdv[i])
                    Rsk_str = unit_manager.print_uncertainty(self.Rsk_mean[i], self.Rsk_stdv[i])
                    Rku_str = unit_manager.print_uncertainty(self.Rku_mean[i], self.Rku_stdv[i])
                elif (uncert=='stdvm'):
                    ym_str  = unit_manager.print_uncertainty(self.ym_mean[i],  self.ym_stdvm[i])
                    Ra_str  = unit_manager.print_uncertainty(self.Ra_mean[i],  self.Ra_stdvm[i])
                    Rq_str  = unit_manager.print_uncertainty(self.Rq_mean[i],  self.Rq_stdvm[i])
                    Rsk_str = unit_manager.print_uncertainty(self.Rsk_mean[i], self.Rsk_stdvm[i])
                    Rku_str = unit_manager.print_uncertainty(self.Rku_mean[i], self.Rku_stdvm[i])
                else:
                    ym_str  = unit_manager.print_uncertainty(self.ym_mean[i],  self.ym_hmd[i])
                    Ra_str  = unit_manager.print_uncertainty(self.Ra_mean[i],  self.Ra_hmd[i])
                    Rq_str  = unit_manager.print_uncertainty(self.Rq_mean[i],  self.Rq_hmd[i])
                    Rsk_str = unit_manager.print_uncertainty(self.Rsk_mean[i], self.Rsk_hmd[i])
                    Rku_str = unit_manager.print_uncertainty(self.Rku_mean[i], self.Rku_hmd[i])

                Szm_str = unit_manager.print_exp(self.Szm[i], 3)
                Sa_str  = unit_manager.print_exp(self.Sa[i],  3)
                Sq_str  = unit_manager.print_exp(self.Sq[i],  3)
                Ssk_str = unit_manager.print_exp(self.Ssk[i], 3)
                Sku_str = unit_manager.print_exp(self.Sku[i], 3)  

                print('\n')
                print('\n')                               
                print('Parameter        '.ljust(20), 'Symbol'.ljust(10), 'Value'.ljust(10))
                print('-----------------'.ljust(20), '------'.ljust(10), '-----'.ljust(10))
                print('Cutoff lenght    '.ljust(20), 'Lc    '.ljust(10), self.Lc[i])
                print('Mean height      '.ljust(20), 'y_mean'.ljust(10), ym_str)
                print('Roughness        '.ljust(20), 'Ra    '.ljust(10), Ra_str)
                print('Roughness        '.ljust(20), 'Rq    '.ljust(10), Rq_str)
                print('Skewness         '.ljust(20), 'Rsk   '.ljust(10), Rsk_str)
                print('Kurtosis         '.ljust(20), 'Rku   '.ljust(10), Rku_str)

                print('Areal mean height'.ljust(20), 'z_mean'.ljust(10), Szm_str)
                print('Areal roughness  '.ljust(20), 'Sa    '.ljust(10), Sa_str)
                print('Areal roughness  '.ljust(20), 'Sq    '.ljust(10), Sq_str)
                print('Areal skewness   '.ljust(20), 'Ssk   '.ljust(10), Ssk_str)
                print('Areal kurtosis   '.ljust(20), 'Sku   '.ljust(10), Sku_str)

        return status, message

