#!/usr/bin/python

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

""" This module contains functions to plot data on 2D and 3D graphs by means of gnuplot.

This package allows to plot data or mathematical expressions inside python,
using the gnuplot program, in the form of 2D or 3D plots.

.. note:: If the package has been correctly installed, a small demo script can
   be run from the system terminal by typing::

   $ python -m gnuplot_manager.demo

   or inside the python console by typing::

   >>> from gnuplot_manager.demo import main
   >>> main()

   The demo will run without need of input and last about 35 seconds.

Full documentation can be found at:

https://github.com/pietromandracci/gnuplot_manager



List of available functions
===========================

Read the doctrings for a complete description of each function.

Create, modify, and close plot windows
--------------------------------------

*new_plot()*
    create a new plot window
*plot_set()*
    modify some properties of a previously created window
*plot_command()*
    send a command to the gnuplot process
*plot_close()*
    close the plot window and terminate the gnuplot process
*plot_close_all()*
    close all the plot windows and terminate all the gnuplot processes 


Plot data
---------

*plot1d()*
    plot a curve from 1d data
*plot2d()*
    plot a curve from 2d data
*plot3d()*
    plot a curve from 3d data
*plot_box()*
    plot a boxplot from 1d data
*plot_curves()*
    plot several curves at the same time

Plot mathematical functions
---------------------------

*plot_function()*
    plot a mathematical expression
*plot_functions()*
    plot several mathematical expression at once


Print plots to files
--------------------

*plot_print()*
    export a plot to a file

*plot_print_all()*
    export to files all the open plots


Reset, clear and refresh plots
------------------------------

*plot_reset()*
    reset a plot: remove all curves and functions
    and the clear the window 
*plot_reset_all()*
    reset all plot windows
*plot_clear()*
    clear the plot area
*plot_clear_all()*
    clear the plot area of all plots
*plot_replot()*
    refresh the plot window
*plot_replot_all()*
    refresh all the plot windows


Utility functions
-----------------

*plot_label()*
    print a string on the plot
*plot_raise()*
    rise the plot window over the other windows on the screen
*plot_lower()*
    lower the plot window under the other windows on the screen
*plot_raise_all()*
    rise all the plot windows    
*plot_lower_all()*
    lower all the plot windows
*plot_check()*
    print the plot properties
*plot_list()*
    print the properties of all plots


The *_PlotWindow* class
=======================

Each plot window is an instance of the *_PlotWindow* class, 
which has several attributes:

*self.window_number*:   
    an integer number that identifies the plot window, [#window_number]_                               
    mainly used to generate unique names for the data files
*self.gnuplot_process*: 
     gnuplot process (instance of *subprocess.Popen*)    
*self.term_type*:
    the type of gnuplot terminal    
*self.plot_type*:
    a string defining the type of plot : '2D', '3D',
    or *None* if the plot window has been closed
*self.n_axes:*
    number of plot axes (2 for 2D plots, 3 for 3D ones)
*self.xmin*:
    minimum of the x-axis (*None* if not set)
*self.xmax*:
    maximum of the x-axis (*None* if not set)
*self.ymin*:
    minimum of the y-axis (*None* if not set)
*self.ymax*:
    maximum of the y-axis (*None* if not set)
*self.zmin*:
    minimum of the z-axis (*None* if not set)
*self.zmax*:
    maximum of the z-axis (*None* if not set)
*self.persistence*:
    *True* if the plot was opened as persistent
*self.title*:
    the window title (*None* if not given)
*self.filename_out*: 
     name of the file to which gnuplot output is redirected
*self.filename_err*:
     name of the file to which gnuplot errors are redirected     
*self.data_filenames*:
     list containing the names of the datafiles related to the
     curves presently plotted on the window
*self.n_volatiles*:
     number of curves that have been plotted using the *volatile=True*
     argument: they are not listed in *self.data_filenames* since
     there are no associated data files
*self.functions*:
     list containing the function strings [#functions]_
*slef.purge*:
     if True, old data files are removed when new data is plotted
     without the *replot=True* option
*self.error*:
     if there was an error while opening the plot
     an error message is stored here

.. [#window_number] Note that this number is *not* the index that identifies the
   plot window inside the *window_list* variable: in fact the former is fixed,
   while the latter may change when other windows are removed from the list.

.. [#functions] Note that no check is made that function strings given to gnuplot 
   are correct. So even wrong ones (which therefore gnuplot has not plotted)
   are listed here.

.. note:: If you modify the plot by sending commands to gnuplot directly, using
   the *plot_command()* function, some of these attributes, such as the number of curves 
   and the list of data files, may not be updated properly.

The *_PlotWindow* class have some methods also, which are called by the functions
of the *functions.py* module to perform their tasks:

*self._command()*
    method used to send commands to gnuplot
*self._quit_gnuplot()*
    method used to terminate the gnuplot process and close the window
*self._add_functions()*
    method used to add one or more mathematical expression
*self._add_curves()*
    method used to add one or more curves from data

.. note:: Since the package is designed to use the functions in the
   *functions.py* module, these methods are not intended to be called directly.
"""

from .global_variables import *
from .functions import *




