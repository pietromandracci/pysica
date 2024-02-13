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

""" Functions used to plot data on plot windows and to do other tasks. """

from os import remove, path, mkdir, rmdir
from time import sleep

from .global_variables import *
from .errors import *
from .funcutils import correct_filename
from .classes import _PlotWindow


# +--------------------------------------------------------------+
# | Functions to create plot windows and modify their properties |
# +--------------------------------------------------------------+
                                                                                 
def new_plot(*,
             xpos=DEFAULT_XPOS,
             ypos=DEFAULT_YPOS,
             terminal=DEFAULT_TERM,        
             width=DEFAULT_WIDTH,
             height=DEFAULT_HEIGHT,
             gnuplot_default=False,
             plot_type=DEFAULT_TYPE,
             xmin=None, xmax=None,
             ymin=None, ymax=None,
             zmin=None, zmax=None,
             format_x=None, format_y=None, format_z=None,
             logx=False, logy=False, logz=False,
             grid=True,
             style=DEFAULT_STYLE,    
             title=None,
             xlabel=None, ylabel=None, zlabel=None,
             persistence=PERSISTENCE,
             redirect_output=REDIRECT_OUT,
             purge=PURGE_DATA,
             options=None):        
    """ This function returns a plot as an istance of the class _PlotWindow.

        Several gnuplot terminal properties can be set.

        Parameters
        ----------

        (All the parameters are keyword-only)

        xpos:            x position of the plot on the screen
        ypos:            y position of the plot on the screen
        terminal:        terminal type, default is DEFAULT_TERM
        width:           width of the plot window, 
                         default is DEFAULT_WIDTH
        height:          height of the plot window, 
                         default is DEFAULT_HEIGHT  
        gnuplot_default: if True, the terminal type, window size and position,
                         are not set, leaving the gnuplot defaults 
                         (wich are *not* the default values 
                         of the parameters in this function)
        plot_type:       a string defining the type of plot
                         '2D' = 2D data
                         '3D' = 3D data  
        xmin, xmax:      minimum and maximum values of the x-axis scale
        ymin, ymax:      minimum and maximum values of the y-axis scale
        zmin, zmax:      minimum and maximum values of the z-axis scale
                         (this one makes sense for 3D plots only)
        format_x:        string defining the format of the x-axis tick labels [1]_
        format_y:        string defining the format of the y-axis tick labels [1]_
        format_z:        string defining the format of the z-axis tick labels [1]_
        logx:            set logarithmic scale on x-axis
        logy:            set logarithmic scale on y-axis
        logz:            set logarithmic scale on z-axis
        xlabel:          x-axis label (string)
        ylabel:          y-axis label (string)
        zlabel:          z-axis label (string, makes sense for 3D plots only)
        grid:            show a grid on the plot
        style:           style used to plot data
        title:           title of the plot (string)
        persistence:     if set to True, the plot will remain visible
                         even after gnuplot has been closed
        redirect_output: True:  save gnuplot otuput and errors to files
                         False: send them to /dev/stdout and /dev/stderr
                         None:  send them to /dev/null
        purge:           if True, old datafiles are deleted each time 
                         new data is plotted on this window
        options:         a string containing other options for the terminal

        .. [1] read gnuplot documentation for format strings

        Initialized data attributes
        --------------------------- 

        A tuple (status, message) is stored in the *error* data attribute
        of the created _PlotWindow instance. 
        In case of inconsistency of the parameters (e.g. xmin>=xmax) an error
        message is stored there, otherwise (0, 'Ok').

        self.error =  one of the tuples defined in the errors.py module

        Returns
        --------
         
        An instance of the _PlotWindow class
    """
    
    (status, message) = NOERROR

    if (gnuplot_installed is False):
        raise FileNotFoundError('gnuplot program was not found')
        
    # Create an instance of the _PlotWindow class
    plot_window = _PlotWindow(gnuplot_default=gnuplot_default,
                              terminal=terminal,                              
                              width=width,
                              height=height,
                              xpos=xpos,
                              ypos=ypos,
                              options=options,                              
                              plot_type=plot_type,
                              title=title,                              
                              persistence=persistence,
                              redirect_output=redirect_output,
                              purge=purge)
       
    # If there were not errors during the instance creation,
    # do not apply the requested settings
    if not plot_window.error[0]:   
        # Apply the requested plot settings
        (status, message)= plot_set(plot_window,
                                    xmin=xmin, xmax=xmax,
                                    ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax,
                                    format_x=format_x,
                                    format_y=format_y,
                                    format_z=format_z,
                                    logx=logx, logy=logy, logz=logz,
                                    xlabel=xlabel,
                                    ylabel=ylabel,
                                    zlabel=zlabel,
                                    grid=grid,
                                    style=style,
                                    replot=False)
        if status: plot_window.error = (status, message)
    
    return plot_window


def plot_set(plot,
             *,
             resetx=False,
             resety=False,
             resetz=False,
             xmin=None, xmax=None,
             ymin=None, ymax=None,
             zmin=None, zmax=None,
             format_x=None, format_y=None, format_z=None,
             logx=None,     logy=None,     logz=None,
             xlabel=None,   ylabel=None,   zlabel=None,
             grid=None,
             style=None,
             replot=False):
    """ Change some settings of a previously initialized plot. 

        Note that setting a range limit to None (e.g. xmin=None) 
        leaves it untouched: e.g. if you have xrange=[0,1] and call 
        plot_set(<plot>, xmin=None) you are not removing the xmin limit.
        It is possibile to remove limits using the resetx, resety 
        and resetz parameters, which tell the function to remove the actual 
        limits *before* it applies the new ones.
        Example:
        if you have xrange=[0,1] and you want to change it to [-inf,2], 
        you can call plot_set(<plot>, resetx=True, xmax=2)

        Parameters
        ----------

        plot:         the plot on which you want to operate, which must 
                      have been previously created by the new_plot function

        (All the following parameters are keyword-only)
      
        resetx:       reset the x-axis range, before applying other settings
        resety:       reset the y-axis range, before applying other settings
        resetz:       reset the z-axis range, before applying other settings
        xmin, xmax:   min and max value of the x-scale
        ymin, ymax:   min and max value of the y-scale
        zmin, zmax:   min and max value of the z-scale
        format_x:     string describing the format of the x-axis tick labels
                      (see gnuplot documentation for format strings)
        format_y:     string describing the format of the y-axis tick labels
                      (see gnuplot documentation for format strings)
        format_z:     string describing the format of the z-axis tick labels
                      (see gnuplot documentation for format strings)
        logx:         True: set logscale; False: unset; None: do nothing
        logy:         True: set logscale; False: unset; None: do nothing
        logz:         True: set logscale; False: unset; None: do nothing
        xlabel:       x-axis label (string)
        ylabel:       y-axis label (string)
        zlabel:       z-axis label (string, makes sense for 3D plots only)
        grid:         True: set grid
        style:        set a specific plot style
        replot:       if True, replot the window to show changes

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not isinstance(plot, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot.plot_type is None): return ERROR_CLOSED_PLOT
    if (plot.n_volatiles
        and ((logx is not None) or (logy is not None) or (logz is not None))):
        return ERROR_LOGSET_VOLATILE    
    if (replot
        and (not plot.n_volatiles)
        and (not plot.data_filenames) and
        (not plot.functions)):
        return ERROR_NO_REPLOT
    if (style is not None ):
        if ( (plot.plot_type == '2D') and (style not in ALLOWED_STYLES_2D) ):
            (status, message) = ERROR_INVALID_STYLE
            message += ' 2D: \"' + str(style) + '\"'
            return status, message
        if ( (plot.plot_type == '3D') and (style not in ALLOWED_STYLES_3D) ):
            (status, message) = ERROR_INVALID_STYLE
            message += ' 3D: \"' + str(style) + '\"'
            return status, message      
        
    if ( (xmin is not None) and (xmax is not None) and (xmin >= xmax) ):
         (status, message) = ERROR_PLOT_PARAMETERS
         message += ': xmin must be lower than xmax'
         return (status, message)
    if ( (ymin is not None) and (ymax is not None) and (ymin >= ymax) ):
         (status, message) = ERROR_PLOT_PARAMETERS
         message += ': ymin must be lower than ymax'
         return (status, message)
    if ( (zmin is not None) and (zmax is not None) and (zmin >= zmax) ):
         (status, message) = ERROR_PLOT_PARAMETERS
         message += ': zmin must be lower than zmax'
         return (status, message)

    # Reset ranges if required
    if resetx:
        plot.xmin = None
        plot.xmax = None
        plot._command('unset xrange')        
    if resety:        
        plot.ymin = None
        plot.ymax = None
        plot._command('unset yrange')        
    if resetz:
        plot.zmin = None
        plot.zmax = None
        plot._command('unset zrange')
        
    # Update the plot window data attributes
    if (xmin is not None): plot.xmin = xmin
    if (xmax is not None): plot.xmax = xmax
    if (ymin is not None): plot.ymin = ymin
    if (ymax is not None): plot.ymax = ymax    
    if (zmin is not None): plot.zmin = zmin
    if (zmax is not None): plot.zmax = zmax
   
    # Set axis ranges in gnuplot
    if ((plot.xmin is not None) and (plot.xmax is not None)):
        plot._command('set xrange ['
                      + str(plot.xmin) + ':'
                      + str(plot.xmax) + ']')
    if ((plot.xmin is not None) and (plot.xmax is     None)):
        plot._command('set xrange ['
                      + str(plot.xmin) + ':*]')       
    if ((plot.xmin is     None) and (plot.xmax is not None)):
        plot._command('set xrange [*'
                      + ':' + str(plot.xmax) + ']')    
    if ((plot.ymin is not None) and (plot.ymax is not None)):
        plot._command('set yrange [' + str(plot.ymin)
                      + ':' + str(plot.ymax) + ']')
    if ((plot.ymin is not None) and (plot.ymax is     None)):
        plot._command('set yrange [' + str(plot.ymin) + ':*]')
    if ((plot.ymin is     None) and (plot.ymax is not None)):
        plot._command('set yrange [*'
                      + ':' + str(plot.ymax) + ']')
    if ((plot.zmin is not None) and (plot.zmax is not None)):
        plot._command('set zrange ['
                      + str(plot.zmin) + ':' + str(plot.zmax) + ']')
    if ((plot.zmin is not None) and (plot.zmax is     None)):
        plot._command('set zrange ['
                      + str(plot.zmin) + ':*]')
    if ((plot.zmin is     None) and (plot.zmax is not None)):
        plot._command('set zrange [*'
                      + ':' + str(plot.zmax) + ']')

    # Set axis labels    
    if (xlabel is not None):   plot._command('set xlabel   ' + '\"' + str(xlabel)   + '\" ')
    if (ylabel is not None):   plot._command('set ylabel   ' + '\"' + str(ylabel)   + '\" ')
    if (zlabel is not None):   plot._command('set zlabel   ' + '\"' + str(zlabel)   + '\" ')

    # Set tick labels formats    
    if (format_x is not None): plot._command('set format x ' + '\"' + str(format_x) + '\" ')
    if (format_y is not None): plot._command('set format y ' + '\"' + str(format_y) + '\" ')
    if (format_z is not None): plot._command('set format z ' + '\"' + str(format_z) + '\" ')

    # Set logarithmic scales
    if (logx is not None):
        if logx: plot._command('set logscale x')
        else:    plot._command('unset logscale x')
    if (logy is not None):
        if logy: plot._command('set logscale y')
        else:    plot._command('unset logscale y')
    if (logz is not None):         
        if logz: plot._command('set logscale z')
        else:    plot._command('unset logscale z')

    # Set grid
    if (grid is not None):                   
        if grid: plot._command('set grid')
        else:    plot._command('unset grid')

    # Set plot style
    if (style is not None):
        plot._command('set style data ' + style)

    # If requested, do a replot
    if replot: plot._command('replot')

    return (status, message)


# +------------------------------------+
# | Give arbitrary commands to gnuplot |
# +------------------------------------+

def plot_command(plot_window, string):
    """ Send a command to the gnuplot process associated to the plot window.

        NOTE: no check is made that the string is a valid gnuplot command

        Parameters
        ----------

        plot_window: the plot window to which the command must be sent
        string:      the string containing the gnuplot command

        Returns
        -------

        One of the tuples defined in the errors.py module
    """
    (status, message) = NOERROR

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT

    plot_window._command(string)

    return (status, message)


# +------------------------+
# | Functions to plot data |
# +------------------------+

def plot1d(plot_window, data, label=None, volatile=VOLATILE,replot=False):
    """ Plots data on a previously initialized plot 

        Parameters
        ----------

        plot_window:  the plot on which data should be plotted, which must 
                      have been previously created by the new_plot function
        data:         y-values of data points (x-values are automatically 
                      calculated by gnuplot)
        label:        string used to describe that set of data in the plot
        volatile:     if True, the data are not written to file, but sent 
                      to gnuplot as volatile data using the special filename '-'
        replot:       add new plot instead of overwriting an old one

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    if (plot_window.plot_type == '3D'):
        (status, message) = ERROR_WRONG_TYPE
        message += ': 1D data on a 3D plot window'
        return status, message
    
    return plot_window._add_curves([ [None, data, label, None] ], volatile, replot)        



def plot2d(plot_window, x_data, y_data, label=None, volatile=VOLATILE, replot=False):
    """ Plots 2D data on a previously initialized plot 

        Parameters
        ----------

        plot_window:  the plot on which data should be plotted, which must 
                      have been previously created by the new_plot function
        x_data:       numpy array containing the x-values of data points
        y_data:       numpy array containing the y-values of data points
        label:        string used to describe that set of data in the plot
        volatile:     if True, the data are not written to file, but sent 
                      to gnuplot as volatile data using the special filename '-'
        replot:       add new plot instead of overwriting an old one

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    if (plot_window.plot_type == '3D'):
        (status, message) = ERROR_WRONG_TYPE
        message += ': 2D data on a 3D plot window'
        return status, message
    
    return plot_window._add_curves([ [x_data, y_data, label, None] ], volatile, replot)        

        
def plot3d(plot_window, x_data, y_data, z_data, label=None, volatile=VOLATILE, replot=False):
    """ Plots 3D data on a previously initialized plot 

        Parameters
        ----------

        plot_window:  the plot on which data should be plotted, which must 
                      have been previously created by the new_plot function
        x_data:       numpy array containing the x-values of data points
        y_data:       numpy array containing the y-values of data points
        z_data:       numpy array containing the z-values of data points
        label:        string used to describe that set of data in the plot
        volatile:     if True, the data are not written to file, but sent 
                      to gnuplot as volatile data using the special filename '-'
        replot:       add new plot instead of overwerwriting an old one

        Returns
        -------

        One of the tuples defined in the errors.py module
    """
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT 
    if (plot_window.plot_type != '3D'):
        (status, message) = ERROR_WRONG_TYPE
        message += ': 3D data on a 2D plot window'
        return status, message
    
    return plot_window._add_curves([ [x_data, y_data, z_data, label, None] ], volatile, replot)


def plot_box(plot_window, data, width=None, label=None, volatile=VOLATILE, replot=False):
    """ Plots a box plot from the given data on a previously initialized plot 

        Parameters
        ----------

        plot_window:  the plot on which data should be plotted, which must 
                      have been previously created by the new_plot function
        data:         numpy array containing the data to be plotted
        width:        width of the box, if not given is decided by gnuplot
        label:        string used to describe that set of data in the plot
        volatile:     if True, the data are not written to file, but sent 
                      to gnuplot as volatile data using the special filename '-'
        replot:       add new plot instead of overwerwriting an old one

        Returns
        -------

        One of the tuples defined in the errors.py module
    """
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT 
    if (plot_window.plot_type != '2D'):
        (status, message) = ERROR_WRONG_TYPE
        message += ': 3D data on a 2D plot window'
        return status, message
    if (width is not None):
        plot_window._command('set boxwidth ' + str(width))
    
    return plot_window._add_curves([ [None, data, label, 'with boxplot'] ], volatile, replot)

        
def plot_curves(plot_window, data_list, volatile=VOLATILE, replot=False):
    """ Plots several 2D or 3D data on a previously initialized plot 

        Parameters
        ----------

        plot_window:  the plot on which data should be plotted, which must 
                      have been previously created by the new_plot function
        data_list:    a list in one of the following forms:
                      2D data: [ [x1 , y1, label1, options1],     
                                 [x2,  y2, label2, options1]     ... ]
                      3D data: [ [x1,  y1, z1, label1, options1], 
                                 [x2,  y2, z2, label2, options2], ... ] 
                      where:
                      - x1 contains the x-coordinates of the points to plot;
                        for 2D plot windows this can be set to None, 
                        which is useful if the data are 1D, as in the case 
                        of boxplots, or if you want gnuplot to 
                        automatically provide x-values.
                        2D plots with and without x values can be mixed:
                        [ [x1,   y1, label1, options1],     
                          [None, y2, label2, options2],... ]  
                      - y1 contains the y-coordinates of the points to plot 
                      - z1 contains the z-coordinates of the points to plot
                        (for 3D plot windows only)  
                      - label1 is a string, that will be used to identify the plot 
                        in the legend, or None
                      - options1 is a string, containing additional options,
                        or None
                      The form must be consisted with the type of plot 
                      that was defined at the plot window creation, 
                      otherwise an error message is returned
        volatile:     if True, the data are not written to file, but sent 
                      to gnuplot as volatile data using the special filename '-'
        replot:       add new plots instead of overwerwriting an old ones

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    return plot_window._add_curves(data_list, volatile, replot)



# +------------------------------------------+
# | Functions to plot mathematical functions |
# +------------------------------------------+

def plot_function(plot_window, func_string, label=None, replot=False):
    """ Plots a function on a previously initialized plot  

        Parameters
        ----------

        plot_window:  the plot, on which the function should be plotted, which must 
                      have been previously created by the new_plot function
        func_string:  a string containing the mathematical expression
                      using x (or x and y for 3D plots) as independent variable
                      Examples:
                      '3 * sin(x) + x**2'
                      '2 * sin(x**2 + y**2)'
                      
                      The string is sent to gnuplot as it is, without checking 
                      that it is a correct functional expression. If it is not, 
                      gnuplot will provide  an error message, which will be
                      - sent to /dev/stderr, if the plot was opened with the 
                        'redirect_output=False' option (default)
                      - written on a file in '/gnuplot.out/errors/' directory,
                        if the plot was opened with 'redirect_output=False'
                      - discarted if the plot was opened with the 
                        'redirect_output=None' option 
        label:        a string used to describe the function in the plot legend:
                      if it is None, the label is not set, and gnuplot will 
                      automatically manage its value: if you don't want a label 
                      to be shown, set it to "" (empy string)

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    return plot_window._add_functions([ [func_string, label, None] ], replot)


def plot_functions(plot_window, func_list, replot=False):
    """ Plots several functions on a previously initialized plot  

        Parameters
        ----------

        plot_window:  the plot, on which the functions should be plotted, which 
                      must have been previously created by the new_plot function
        func_list:    a list in one the following form:
                      [ [function1, label1, options1], 
                        [function2, label2, options1], ... ]
                      where:
                      - function1 is a string defining the function 
                        to be plotted
                        e.g. '2*x**2+3*cos(x)' or 'sin(x**2+y**2) 
                      - label1 is a string, that will be used to identify 
                        the plot in the legend, or None
                      - options1 is a string, containing additional options
                        or None                       
                      The string is sent to gnuplot as it is, without checking 
                      that it is a correct functional expression. If it is not, 
                      gnuplot will provide  an error message, which will be
                      - sent to /dev/stderr, if the plot was opened with the 
                        'redirect_output=False' option (default)
                      - written on a file in '/gnuplot.out/errors/' directory,
                        if the plot was opened with 'redirect_output=False'
                      - discarted if the plot was opened with the 
                        'redirect_output=None' option 
                      
        Returns
        -------

       One of the tuples defined in the errors.py module
    """

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    return plot_window._add_functions(func_list, replot)


# +----------------------------------+
# | Functions to print plots to file |
# +----------------------------------+

def plot_print(plot_window,
               terminal=DEFAULT_PRINT_TERM,
               dirname=None,
               filename=None,
               options=None):
    """ Export the plot to a file.

        The plot is saved as an image, the type of which depends on the choosen 
        terminal (default is png), in the present working directory.

        Parameters
        ----------

        plot_window:    the plot window which must be printed
        terminal:       type of terminal to use to print the image 
                        (e.g. 'png', 'jpeg', 'gif') mus be one of the
                        terminal listed in the PRINT_TERMINALS global variable
        dirname:        name of the directory where file must be saved;
                        if not given, is the current working directory
        filename:       the file to which the image must be saved
                        if not given, or set to None, a default one is used
        options:        a string with terminal options,
                        read gnuplot documentations for help

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    if ( (not plot_window.n_volatiles)
         and (not plot_window.data_filenames)
         and (not plot_window.functions) ):
        return ERROR_NO_REPLOT    
    
    if terminal not in PRINT_TERMINALS: return ERROR_UNKNOWN_TERM   

    if (filename is None):
        filename = ( PRINT_FILENAME
                     + str(plot_window.window_number)
                     + PRINT_EXT[terminal] )
    else:
        filename = correct_filename(str(filename))
        
    if (dirname is not None):
        dirname = correct_filename(str(dirname))
        if not path.exists(dirname): mkdir(dirname)
        filename = path.join(dirname, filename)
        
    # Save the current terminal settings
    plot_window._command('set terminal push')

    # Set the terminal for printing 
    command_string =  'set terminal '
    command_string += terminal
    if (options is not None):
        command_string += ' ' + str(options)    
    plot_window._command(command_string)

    # Set the filename 
    command_string = 'set output '
    command_string += '\"' + filename  + '\"'

    # Sent the command to change terminal
    plot_window._command(command_string)

    # Print the plot
    plot_window._command('replot')

    # Recover the previous terminal settings
    plot_window._command('set terminal pop')

    return (status, message)


def plot_print_all(terminal=DEFAULT_PRINT_TERM, dirname=None, options=None):
    """ Export all the open plot windows to files.

        Parameters
        ----------

        terminal:       type of terminal to use to print the plots 
                        e.g. 'png', 'jpeg', 'gif'
        dirname:        name of the directory where files must be saved;
                        if not given, is the current working directory
        options:        a string with terminal options,
                        read gnuplot documentations for help

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    
       
    for plot_window in window_list:
        plot_print(plot_window,
                   terminal=terminal,
                   dirname=dirname,
                   options=options)

    return (status, message)    


# +---------------------------------------------+
# | Functions to reset, clear and refresh plots |
# +---------------------------------------------+

def _best_pos(xmin, xmax, d=1):
    """ Return coordinate to plot a dot to force gnuplot plot axis.

        Search for a suitable coordinate on an axis
        to plot a dot on the graph to force gnuplot 
        draw the axes.

        Parameters
        ----------

        xmin: minimum value of the range
        xamx: maximum value of the range
        d:    distance from the extreme of the range
              at which plot is plotted if some conditions
              are true

        Returns
        -------

        The coordinate af the dot on the given axis
    """
    
    # Range is [-inf,+inf]
    if  ((xmin is None) and (xmax is None)): x = d
    # Range is [xmin,+inf]
    elif (xmax is None):
        if (xmin > 0):                       x = xmin + d
        else:                                x = d
    # Range is [-inf,xmax]
    elif (xmin is None):
        if (xmax > 0):                       x = xmax / 2
        else:                                x = xmax - d
    # Range is [xmin,xmax]    
    else:
        if (xmax > 0):                       x = xmax / 2
        else:                                x = xmax - d                                             

    return x   


def plot_reset(plot_window, plot_axes=True):
    """ Resets the plot window.

        Remove all the functions and curves and clear the window,
        if requested, plot the axis also.

        Parameters
        ----------

        plot_window:  the window (instance of the _PlotWindow class)
        plot_axes:    if True, the axes are plotted on the window, 
                      after it has been cleared

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT

    # Reset functions and curves
    if plot_window.purge:
        for filename in plot_window.data_filenames:
            try:
                remove(filename)
            except FileNotFoundError:
                pass    
    plot_window.data_filenames.clear()
    plot_window.functions.clear()
    plot_window.n_volatiles = 0

    # Remove labels
    plot_window._command('unset label')

    # Clear the window
    plot_window._command('clear')
        
    # Force gnuplot plotting the axes
    if plot_axes:
        x = _best_pos(plot_window.xmin, plot_window.xmax)
        if (x <= 0): plot_window._command('unset logscale x')
        
        y = _best_pos(plot_window.ymin, plot_window.ymax)
        if (y <= 0): plot_window._command('unset logscale y')
        
        if (plot_window.plot_type == '3D'):
            z = _best_pos(plot_window.zmin, plot_window.zmax)
            if (z <= 0): plot_window._command('unset logscale z')
            
        if (plot_window.plot_type == '3D'):
            string  = ('splot \"-\" with dots title \"\"\n'
                       + str(x)   + ' ' + str(y)   + ' ' + str(z) + '\n'
                       + str(2*x) + ' ' + str(2*y) + ' ' + str(2*z) + '\n'
                       + 'e\n')
        else:
            string  = ('plot \"-\" with dots title \"\"\n'
                       + str(x) + ' ' + str(y) + '\n'
                       + 'e\n')
        plot_window._command(string)        
    
    return status, message


def plot_reset_all(plot_axes=True):
    """ Resets all the plot windows.

        This function removes all the functions and curves
        and clears all the plot windows.

        Parameters
        ----------

        plot_axes: if True, the axes are plotted on the windows

        Returns
        -------

        One of the tuples defined in the errors.py module

    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    

    for plot_window in window_list:
        plot_reset(plot_window, plot_axes)
        
    return status, message


def plot_clear(plot_window):
    """ Clears the plotting area of the plot window

        Parameters
        ----------

        plot_window: the window to be cleared

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT  
    
    plot_window._command('clear')
   
    return status, message


def plot_clear_all():
    """ Clears the plotting area of all plots 

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    

    for plot_window in window_list:
        plot_window._command('clear')
        
    return status, message

        
def plot_replot(plot_window):
    """ Redraws data on the given plot

        Parameters
        ----------

        plot_window: the window (instance of the _PlotWindow class) to be replotted

        Returns
        -------

        One of the tuples defined in the errors.py module
    """
    
    (status, message) = NOERROR

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    if ( (not plot_window.n_volatiles)
         and (not plot_window.data_filenames)
         and (not plot_window.functions) ):
        return ERROR_NO_REPLOT
    else:
        plot_window._command('replot')
        
    return status, message


def plot_replot_all():
    """ Redraws all the plots 

        Returns
        -------

        One of the tuples defined in the errors.py module
    """
    
    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    
    
    for plot_window in window_list:
        (status1, message1) = plot_replot(plot_window)
        if status1: (status, message) = (status1, message1)

    if status:
        message += ' on some plot windows' 
    
    return status, message


# +-------------------+
# | Utility functions | 
# +-------------------+

def plot_label(plot_window,*, x=1, y=1, z=1, label=None, erase=False, replot=False):
    """ Prints a label on a previously initialized plot.

        The label is not shown immediately, but at the first plotting action.

        Parameters
        ----------

        plot_window:  the plot on which the label should be printed, it must 
                      have been previously created by the new_plot function

        (All the following parameters are keyword-only)       

        x, y, z:      the position on the plot at which the label must 
                      be created, expressed in characters
        label:        a string containing the label to be printed
                      if it is not a string, will be converted to string
                      if it is set to None, no label is shown
                      this can be useful with erase=True to remove all labels
        erase:        if set to True, all previously set labels are removed
                      before setting the new one
        replot:       if set to True, the plot is replotted to show the label
                      immediately
        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR
    
    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    if erase: plot_window._command('unset label')
    if label is not None:
        plot_window._command('set label \"' + str(label) + '\"'
                             + ' at character '
                             + str(x) + ',' + str(y) + ',' + str(z))
                            #+ 'left norotate front')
    plot_window._command('show label')
    if replot: plot_window._command('replot')
    
    return status, message


def plot_raise(plot_window):
    """ Raises the plot window on the desktop

        Parameters
        ----------

        plot_window: the window (instance of the _PlotWindow class) to be rised

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR 

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    plot_window._command('raise')
    
    return status, message    

    
def plot_lower(plot_window):
    """ Lowers the plot window on the desktop

        Parameters
        ----------

        plot_window: the window (instance of the _PlotWindow class)

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT
    
    plot_window._command('lower')
    
    return status, message

    
def plot_raise_all():
    """ Raises all the plot windows on the desktop 

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    
    
    for plot_window in window_list: plot_window._command('raise')
    
    return status, message    


def plot_lower_all():
    """ Lowers all the plot windows on the desktop 

        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS    
    
    for plot_window in window_list: plot_window._command('lower')
    
    return status, message
    
    
# +--------------------------------------------------+
# | Functions to show information about active plots |
# +--------------------------------------------------+ 

def plot_check(plot_window, expanded=False, printout=True, getstring=False):
    """ Prints information about the plot window 

        Parameters
        ----------

        plot_window: the plot (instance of _PlotWindow class)
                     about which info is requested

        expanded:    if set to False, printed information includes
                     - window number
                     - type of terminal
                     - persistence status
                     - purge status
                     - type of window ('2D', 'histogram', '3D')
                     - window title (if given)
                     - number of functions plotted
                     - number of curves plotted
                     - number of volatile curves plotted
                     - x-axis range
                     - y-axis range
                     - z-axis range (for 3D windows only)
                     if set to True, the following information is added
                     - PID of the gnuplot process
                     - names of data files to which the output and errors
                       provided by gnuplot are written;
                       depending on the option 'redirect_output' used at 
                       plot creation, they can be
                        - names of real files
                        - '/dev/stdout' and /dev/stderr' 
                        - 'dev/null'
                     - for each plotted function, the function string [1]_
                     - for each plotted curve, the data file name
        printout:    if True, the info are written on screen
        getstring:   if True, returns a string with the info

        .. [1] Note that no check is made that function strings given to 
               gnuplot are correct. So even wrong ones (which therefore 
               gnuplot has not plotted) are listed here.

        Returns
        -------

        If getstring is True, a string with the plot information
        if getstring is False, or there is an error, 
        one of the tuples defined in the errors.py module is returned
    """

    (status, message) = NOERROR

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT    
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT   
    
    string = ''
    string += 'Window index:         ' + str(window_list.index(plot_window)) + '\n'
    string += 'Window number:        ' + str(plot_window.window_number) + '\n'
    string += 'Terminal type:        ' + '\"' + plot_window.term_type + '\"' + '\n'
    string += 'Persistence:          ' + '\"' + str(plot_window.persistence) + '\"' + '\n'
    string += 'Purge:                ' + '\"' + str(plot_window.purge) + '\"' + '\n'
    string += 'Window type:          ' + '\"' + str(plot_window.plot_type) + '\"' + '\n'   
    string += 'Window title:         ' + '\"' + str(plot_window.title) + '\"' + '\n'
    string += 'Number of functions:  ' + str(len(plot_window.functions)) + '\n'
    string += 'Number of curves:     ' + str(len(plot_window.data_filenames)) + '\n'
    string += 'Number of volatiles:  ' + str(plot_window.n_volatiles) + '\n'
    string += 'X-axis range:         ' + '[' + str(plot_window.xmin) + ',' + str(plot_window.xmax) + ']' + '\n'
    string += 'Y-axis range:         ' + '[' + str(plot_window.ymin) + ',' + str(plot_window.ymax) + ']' + '\n'
    if (plot_window.plot_type == '3D'):
        string += 'Z-axis range:         ' + '[' + str(plot_window.zmin) + ',' + str(plot_window.zmax) + ']' + '\n'
    if expanded:
        string += 'Gnuplot process PID:  ' + str(plot_window.gnuplot_process.pid) + '\n'
        string += 'Gnuplot output file:  ' + '\"' + plot_window.filename_out + '\"' + '\n'
        string += 'Gnuplot errors file:  ' + '\"' + plot_window.filename_err + '\"' + '\n'
        if plot_window.functions:                 
            string += 'Functions\n'
            for i in range(len(plot_window.functions)):
                string += '#' + str(i).rjust(3) + ': \"' + plot_window.functions[i] +'\"' + '\n'
        if plot_window.data_filenames: 
            string += 'Curves\n'            
            for i in range(len(plot_window.data_filenames)):
                string += '#' + str(i).rjust(3) + ': \"' + plot_window.data_filenames[i] +'\"' + '\n'
    if printout: print(string)

    if getstring: return string
    else:         return (status, message)


def plot_list(expanded=False, printout=True, getstring=False):
    """ Prints a list of all active plots 

        Parameters
        ----------

        expanded:    if set to False, printed information includes
                     - window nunmber
                     - window title (if given)
                     - type of window ('2D', 'histogram', '3D')
                     - number of curves plotted
                     - number of functions plotted
                     - for each curve, the name of associated data file
                     if set to True, the following information is added
                     - _PlotWindow instance
                     - gnuplot process (Popen instance)
                     - PID of the gnuplot process
                     - names of data files to which the output and errors
                       provided by gnuplot are written;
                       depending on the option 'redirect_output' used at 
                       plot creation, they can be
                        - names of real files
                        - '/dev/stdout' and /dev/stderr' 
                        - 'dev/null'
                     - for each plotted curve, the data file instance
        printout:    if True, the info are written on screen
        getstring:   if True, returns a string with the info

        Returns
        -------

        If getstring is True, a string with the plot information
        if getstring is False, or there is an error, 
        one of the tuples defined in the errors.py module is returned
    """

    (status, message) = NOERROR

    if not window_list: return ERROR_NO_PLOTS
   
    string = ''
    for plot_window in window_list:
        string += plot_check(plot_window, expanded, False, True)
        string += '\n'

    if printout: print(string)
     
    if getstring:
        return string
    else:
        return status, message


# +---------------------------------+
# | Functions to close active plots |
# +---------------------------------+
            
def plot_close(plot_window, keep_output=False, delay=None):
    """ Closes a plot windows and exits from the related gnuplot process 

        If the plot window was created calling the new_plot function with 
        the *purge=True* argument (which is the default) the data files 
        of the presently plotted curves are deleted.

        Parameters
        ----------

        plot_window:  the plot to be closed (instance of _PlotWindow class)
        keep_output:  do not remove the files to which gnuplot output and 
                      errors have been redirected, even if the window
                      was created with the *purge=True* argument                 
        delay:        an int or float number, defining the time to
                      wait (in seconds) before deleting data files;
                      this can be useful if a plot is created with
                      the *persist=True* option and then immediately closed,
                      since the gnuplot process could not be able to plot the
                      data before the datafiles are deleted
        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR

    if not isinstance(plot_window, _PlotWindow): return ERROR_NOT_A_PLOT 
    if (plot_window.plot_type is None): return ERROR_CLOSED_PLOT    

    # If window has purge attribute, delete the files associated to this plot
    if plot_window.purge:
        # Remove data files but, if requested, wait a while before doing it 
        if (delay is not None): sleep(delay)
        for filename in plot_window.data_filenames:
            try:
                remove(filename)
            except FileNotFoundError:
                (status, message) = ERROR_FILE
                message += ('one or more datafiles for window # '
                            + str(plot_window.window_number)
                            + ' were not found')
        if not keep_output:        
            # Remove gnuplot output file
            if (plot_window.filename_out not in ('/dev/stdout', '/dev/null')):
                try:
                    remove(plot_window.filename_out)
                except FileNotFoundError:
                    (status, message) = ERROR_FILE
                    message += ('the gnuplot output file for window # '
                                + str(plot_window.window_number)
                                + ' was not found')
            # Remove gnuplot error file                
            if (plot_window.filename_err not in ('/dev/stderr', '/dev/null')):
                try:
                    remove(plot_window.filename_err)
                except FileNotFoundError:
                    (status, message) = ERROR_FILE
                    message += ('the gnuplot error file for window # '
                                + str(plot_window.window_number)
                                + ' was not found')
                    
    # Shut down the associated gnuplot process
    plot_window._quit_gnuplot()

    # Set to None the plot_type attribute
    plot_window.plot_type = None
    
    # Remove the _PlotWindow instance from the plot list
    window_list.remove(plot_window)
  
    return status, message


def plot_close_all(keep_output=False, delay=None, purge_dir=PURGE_DIR):
    """ Closes all plots.

        For plot windows which were created calling the new_plot function with 
        the *purge=True* argument (which is the default), the data files 
        of the presently plotted curves are deleted.

        Parameters
        ----------

        keep_output:  do not remove the files to which gnuplot output and 
                      errors have been redirected, even if the window
                      was created with the *purge=True* argument
        delay:        an int or float number, defining the time to
                      wait (in seconds) before deleting data files;
                      this can be useful if a plot is created with
                      the persist=True option and then immediately closed,
                      since the gnuplot process could not be able to plot the
                      data before the datafiles are deleted
        purge_dir:    remove the directories in which files are stored; 
                      this works only if they are empty and is uneffective
                      if keep_output has bee set to True
        Returns
        -------

        One of the tuples defined in the errors.py module
    """

    (status, message) = NOERROR
    
    while window_list:
        plot_window = window_list[-1]
        (status1, message1) = plot_close(plot_window, keep_output, delay)
        # If there was an error closing the window, store it
        if status1: (status, message) = (status1, message1)

    if (purge_dir and not keep_output):
        # Try to remove data file subdirectory
        try:
            rmdir(DIRNAME_DATA)
        except FileNotFoundError:
            pass
        except OSError:
            (status, message) = ERROR_FILE
        # Try to remove subdirectory with gnuplot output            
        try:
            rmdir(DIRNAME_OUT)
        except FileNotFoundError:
            pass
        except OSError:
            (status, message) = ERROR_FILE
        # Try to remove subdirectory with gnuplot errors            
        try:
            rmdir(DIRNAME_ERR)
        except FileNotFoundError:
            pass
        except OSError:
            (status, message) = ERROR_FILE
        # Try to remove main directory
        try:
            rmdir(DIRNAME)
        except FileNotFoundError:
            pass
        except OSError:
            (status, message) = ERROR_FILE
                
    if status:
        message += ': some files or directories could not be deleted'                
  
    return status, message
