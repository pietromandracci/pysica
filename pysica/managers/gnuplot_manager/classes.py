# COPYRIGHT (c) 2020-2022 Pietro Mandracci

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

""" The _PlotWindow class, used to define a plot window. """

from subprocess import Popen, PIPE, DEVNULL, TimeoutExpired
from os import remove, path, mkdir

from .global_variables import *
from .errors import *
from .funcutils import *


# +---------+
# | Classes |
# +---------+
       
class _PlotWindow:
    """ This class defines a single plot window. """
    
    def __init__(self,
                 gnuplot_default=False,
                 terminal=DEFAULT_TERM,
                 width=DEFAULT_WIDTH,
                 height=DEFAULT_HEIGHT,
                 xpos=DEFAULT_XPOS,
                 ypos=DEFAULT_YPOS,
                 options=None,
                 plot_type=DEFAULT_TYPE,
                 title=None,
                 persistence=PERSISTENCE,
                 redirect_output=REDIRECT_OUT,
                 purge=PURGE_DATA):
        """ An instance of this class contains a plot window.

            Parameters
            ----------
            gnuplot_default: if it is set to True, the terminal 
                             type, window size and position,
                             are not set, leaving the gnuplot defaults 
                             (which are *not* the default values 
                             of the parameters in this method)
            terminal:        terminal type, must be one of
                             the strings in the TERMINALS tuple
            width:           width of the plot window
            height:          height of the plot window
            xpos:            x position of the plot on the screen
            ypos:            y position of the plot on the screen
            options:         a string containing other options for the terminal
            plot_type:       a string defining the type of plot
                             must be one of the strings in PLOT_TYPES
            title:           a string used as the title of the plot
            persistence:     if True, the window on the screen will not close 
                             after the gnuplot process is terminated
            redirect_output: True:  save gnuplot otuput and errors to files
                             False: send then to /dev/stdout and /dev/stderr
                             None:  send them to /dev/null
            purge:           if True, old datafiles are deleted each time 
                             new data is plotted on this window

            Initialized data attributes
            ---------------------------

            self.window_number:   integer number that identifies the plot window
            self.gnuplot_process: gnuplot process (instance of subprocess.Popen)
            self.term_type:       the type of gnuplot terminal started 
                                  (string, e.g. 'x11' or 'wxt')
            self.persistence:     True if the window is persistent
            self.filename_out:    file to which gnuplot output is redirected 
            self.filename_err:    file to which gnuplot error messages
                                  are redirected
            self.plot_type:       a string defining the type of plot
            self.n_axes:          number of plot axes (2 for 2D plots, 3 for 3D ones)
            self.xmin:            minimum of x-axis range (initialized to None)
            self.xmax:            maximum of x-axis range (initialized to None)
            self.ymin:            minimum of y-axis range (initialized to None)
            self.ymax:            maximum of y-axis range (initialized to None)
            self.zmin:            minimum of z-axis range (initialized to None)
            self.zmax:            maximum of z-axis range (initialized to None)
            self.title:           the window title (None if not given)
            self.data_filenames:  list containing the names of the data files
            self.functions:       list containing the functions plotted
            self.n_volatiles      number of curves plotted as volatile
                                  (passing the data inline, without writing to file)
            self.purge:           if True, old datafiles are deleted each time 
                                  new data is plotted
            self.error:           one of the tuples in the errors.py module
        """

        self.error = NOERROR
        
        # Check that terminal and plot type are among the known ones;
        # otherwise, use defaults and store an error message in self.error
        if terminal in TERMINALS:
            self.term_type = terminal
        else:
            self.term_type = DEFAULT_TERM
            (status, message) = ERROR_UNKNOWN_TERM
            message += (' \"' + str(terminal) + '\"'
                        + ', using default  \"'
                        + DEFAULT_TERM + '\"')
            self.error = (status, message)
        if plot_type in PLOT_TYPES: 
            self.plot_type = plot_type
        else:
            self.plot_type = DEFAULT_TYPE
            (status, message) = ERROR_UNKNOWN_TYPE
            message += (' \"' + str(plot_type) + '\"'
                        + ', using default \"'
                        + DEFAULT_TYPE + '\"' )  
            self.error = (status, message)

        # Store the persistence status
        self.persistence = persistence

        # Store the purge status
        self.purge = purge
        
        # Store the window title (may be None)
        self.title = title        

        # Store the number of axis     
        if (self.plot_type == '3D'): self.n_axes = 3
        else:                        self.n_axes = 2
        
        # Assign a unique number to this plot window
        # to be used to create unique names for datafiles
        self.window_number = None
        for i in range(len(window_list)):
            if (i != window_list[i].window_number):
                self.window_number = i
                break
        if (self.window_number is None): self.window_number = len(window_list)
        
        # Add this plot window to the list of the active ones
        window_list.append(self)
        
        # Initialize the lists where data filenames 
        # and functon strings will be stored
        self.data_filenames = []
        self.functions      = []
        self.n_volatiles    = 0

        # Initialize the axis ranges
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None            
                
        # Start the gnuplot process to plot on this window
        if persistence: exec_list = ['gnuplot','-p']
        else:           exec_list = ['gnuplot']        
        if redirect_output is True:
            # If redirection of gnuplot output to files is requested
            # Create the directories to save gnuplot output and errors
            if not path.exists(DIRNAME):     mkdir(DIRNAME)
            if not path.exists(DIRNAME_OUT): mkdir(DIRNAME_OUT)
            # Create the filenames
            if not path.exists(DIRNAME_ERR): mkdir(DIRNAME_ERR)
            self.filename_out = path.join( DIRNAME_OUT, FILENAME_OUT
                                           + '_w_' + str(self.window_number) )
            self.filename_err = path.join( DIRNAME_ERR, FILENAME_ERR
                                           + '_w_' + str(self.window_number) )
            if self.title is not None:
                self.filename_out += '(' + correct_filename(self.title) + ')'
                self.filename_err += '(' + correct_filename(self.title) + ')'
            # Start the gnuplot process with output redirected to files
            self.gnuplot_process = Popen(exec_list,
                                         stdin=PIPE,
                                         stdout=open(self.filename_out, 'w'),
                                         stderr=open(self.filename_err, 'w'),
                                         universal_newlines=True)
        elif redirect_output is False:
            # Start the gnuplot process without redirection
            self.filename_out = '/dev/stdout'
            self.filename_err = '/dev/stderr'               
            self.gnuplot_process = Popen(exec_list,
                                         stdin=PIPE,
                                         universal_newlines=True)            
        else:
            # Start the gnuplot process silenced
            self.filename_out = '/dev/null'
            self.filename_err = '/dev/null'        
            self.gnuplot_process = Popen(exec_list,
                                         stdin=PIPE,
                                         stdout=DEVNULL,
                                         stderr=DEVNULL,
                                         universal_newlines=True)
        if not gnuplot_default:
            command_string = ( 'set terminal ' + str(self.term_type)
                               + ' size '      + str(width) + ',' + str(height)
                               + ' position '  + str(xpos)  + ',' + str(ypos)   )
            if (options is not None):
                command_string += ' ' + str(options)
            self._command(command_string)
        if title is not None:
            self._command('set title \"' + correct_filename(self.title) + '\"')
            
            
    def _command(self, command_string):
        """ Send a command to the gnuplot process.           

            Parameters
            ----------

            command_string: string containing a gnuplot command
                            NOTE: no check is made that the string 
                            is a valid gnuplot command
        """

        self.gnuplot_process.stdin.write(command_string + EOL)
        self.gnuplot_process.stdin.flush()

        
    def _quit_gnuplot(self):
        """ Terminate the gnuplot process associated to this window. 

            Returns
            -------

            The last output provided by the gnuplot process
        """

        try:
            gnuplot_last_output = self.gnuplot_process.communicate(
                input='quit', timeout=TIMEOUT_QUIT)
        except TimeoutExpired:
            self.gnuplot_process.kill()
            gnuplot_last_output = self.gnuplot_process.communicate()

        return gnuplot_last_output                    
        

    def _add_functions(self, function_list, replot=False):
        """ Add one or more functions to a plot window.

            Parameters
            ----------

            function_list: a list in the following form:
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
            replot:       if True, the functions are plotted without erasing 
                          previously plotted ones
            Returns
            -------
            
            One of the tuples defined in the errors.py module
        """

        (status, message) = NOERROR

        if (self.n_volatiles and replot): return ERROR_REPLOT_VOLATILE
        if not function_list: return ERROR_NO_FUNCTION 
        if (replot and (not self.data_filenames) and (not self.functions)):
            return ERROR_NO_REPLOT
        # Check consistency of the whole list
        for i in range(len(function_list)):
            if ( len(function_list[i]) != 3):
                (status, message) = ERROR_PLOT_PARAMETERS
                if len(function_list) > 1:
                    message += ' in list item # ' + str(i)
                message += ': 3 items expected, given ' + str(len(function_list[i]))
                return (status, message)
                   
        # If it is not a replot, reset the plot window   
        if not replot:
            # If purge status is True, remove old datafiles
            if self.purge:
                for filename in self.data_filenames:
                    try:
                        remove(filename)
                    except FileNotFoundError:
                        pass
            # Clear the list of datafiles and functions
            # and reset the number of volatiles
            self.data_filenames.clear()
            self.functions.clear()
            self.n_volatiles = 0
            
        # Beginning of the command string to send to gnuplot
        if replot:                     command_string = 'replot '
        elif (self.plot_type == '3D'): command_string = 'splot '
        else:                          command_string = 'plot '

        command_list = []
              
        for function_item in function_list:
            # Initialize the command string for this fucntion
            function_command_string = ''
            # Add the function to the function list
            self.functions.append(str(function_item[0]))
            
            # Read the label, if given
            if (function_item[1] is None): label = None
            else:                          label = str(function_item[1])
            
            # Add  to the command string the plot command for this function
            function_command_string += str(function_item[0])
            if (label is not None):
                function_command_string += ' title \"' + str(label) + '\"'
            if (function_item[2] is not None):
                function_command_string += ' ' + str(function_item[2])
                
            command_list.append(function_command_string)

        command_string += ', '.join(command_list)
                        
        # Send command string to gnuplot   
        self._command(command_string)

        return status, message        
        

    def _add_curves(self, data_list, volatile=VOLATILE, replot=False):
        """ Add one or more curves from data to the plot window. 

            Parameters
            ----------

            data_list:    a list in one of the following forms:
                          2D data: [ [x1 ,  y1, label1, options1],     
                                     [x2,   y2, label2, options2],... ]
                          3D data: [ [x1,   y1, z1, label1, options1], 
                                     [x2,   y2, z2, label2, options2], ... ] 
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
                          - options1 is a string, containing additional options
                            or None
                          The form must be consisted with the type of plot 
                          that was defined at the plot window creation, 
                          otherwise an error message is returned
            volatile:     if True, the data are not written to file, but sent 
                          to gnuplot as volatile data using the special filename '-'
            replot:       if True, the curves are plotted without erasing 
                          previously plotted ones

            Returns
            -------

            One of the tuples defined in the errors.py module
        """

        (status, message) = NOERROR

        if (self.n_volatiles and replot): return ERROR_REPLOT_VOLATILE      
        if not data_list: return ERROR_NO_DATA
        if (replot and (not self.data_filenames) and (not self.functions)):
            return ERROR_NO_REPLOT
        
        # Check consistency of the each data list item
        for i in range(len(data_list)):
            # Check the number of data given
            if (len(data_list[i]) != self.n_axes + 2):
                (status, message) = ERROR_PLOT_PARAMETERS
                if len(data_list) > 1: message += ' in list item # ' + str(i)
                message += ( ': ' +  str(self.n_axes + 2)
                             + ' items expected, given '
                             + str(len(data_list[i])) )
                return (status, message)            
            
            # Check if xdata are missing in this item            
            if (data_list[i][0] is None):
                xmissing = True
                len_x = None
                # Missing x-data are not allowed in 3D plots
                if (self.plot_type == '3D'):
                    (status, message) = ERROR_WRONG_TYPE
                    if (len(data_list) > 1): message += ' in list item # ' + str(i)
                    message += ': missing x-data in a 3D plot window'
                    return status, message
            else:
                xmissing = False
                
            # Transform single numbers to lists
            if not xmissing:
                try:
                    len_x = len(data_list[i][0])
                except TypeError:
                    data_list[i][0] = [ data_list[i][0] ]
                    len_x = 1
            try:
                len_y = len(data_list[i][1])
            except TypeError:
                data_list[i][1] = [ data_list[i][1] ]
                len_y = 1
            if (self.n_axes == 3):
                # Transform single z-numbers to lists
                try:
                    len_z = len(data_list[i][2])
                except TypeError:
                    data_list[i][2] = [ data_list[i][2] ]
                    len_z = 1
                    
            # Check consistency of the data lengths
            if ( (not xmissing) and (len_x != len_y) ):
                (status, message) = ERROR_PLOT_PARAMETERS
                if len(data_list) > 1:  message += ' in list item # ' + str(i)
                message += ': x and y data have different sizes'
                return status, message
            elif ( (self.n_axes == 3) and (len_z != len_y) ):
                (status, message) = ERROR_PLOT_PARAMETERS
                if len(data_list) > 1:  message += ' in list item # ' + str(i)
                message += ': y and z data have different sizes'
                return status, message               
            elif (len_y == 0):
                (status, message) = ERROR_PLOT_PARAMETERS
                if len(data_list) > 1: message += ' in list item # ' + str(i)
                message += ': zero size data given'
                return status, message

        # Create the data files directory, if it doesn't exist,
        # unless a volatile plot has been requested
        if not volatile:
            if not path.exists(DIRNAME):      mkdir(DIRNAME)
            if not path.exists(DIRNAME_DATA): mkdir(DIRNAME_DATA)
                   
        # If it is not a replot, reset the plot window   
        if not replot:
            # If purge status is True, remove old datafiles
            if self.purge:
                for filename in self.data_filenames:
                    try:
                        remove(filename)
                    except FileNotFoundError:
                        pass
            # Clear the list of datafiles and functions
            # and reset the number of volatiles
            self.data_filenames.clear()
            self.functions.clear()
            self.n_volatiles = 0            
        
        # Beginning of the command string to send to gnuplot
        if replot:                     command_string = 'replot '
        elif (self.plot_type == '3D'): command_string = 'splot '
        else:                          command_string = 'plot '

        # Initialize the list in which commands for each curve
        # will be stored
        command_list = []
        
        # Save the old number of curves before adding the new ones
        n_curves = len(self.data_filenames)       
        for data_item in data_list:
            # Initialize the string for this curve
            curve_command_string = ''
            
            # Read the curve label, if given
            if (data_item[self.n_axes] is None):
                label = None
            else:
                label = str(data_item[self.n_axes])

            # Read additional options, if given
            if (data_item[self.n_axes+1] is None):
                options = None
            else:
                options = str(data_item[self.n_axes+1])
                        
            # Define the unique filename for the data file
            # or use the gnuplot special filename '-' for volatile data
            if volatile:
                filename = FILENAME_VOLATILE
            else:
                curve_number = n_curves + i
                name_string = FILENAME_DATA 
                name_string += '_w' + str(self.window_number)
                if (self.plot_type == '3D'): name_string += '_3D'
                elif (data_item[0] is None): name_string += '_1D'
                else:                        name_string += '_2D'        
                if self.title is not None:
                    name_string += '(' + correct_filename(self.title) + ')'
                name_string += '_c' + str(curve_number)
                if label is not None:
                    name_string += '(' + correct_filename(label) + ')' 
                name_string += FILENAME_DATA_EXT
                filename = path.join( DIRNAME_DATA, name_string )

            # Add filename to list, or increase volatile data counter
            if volatile:
                self.n_volatiles += 1
            else:
                self.data_filenames.append(filename)

            # Save data to data files, unless they are volatile
            if not volatile:
                if (self.plot_type == '3D'):
                    data_file_3d( data_item[0],
                                  data_item[1],
                                  data_item[2],
                                  filename )

                elif (data_list[i][0] is None):
                    data_file_1d( data_item[1],
                                  filename )

                else:
                    data_file_2d( data_item[0],
                                  data_item[1],
                                  filename )

            # Add plot command to the command list
            curve_command_string += '\"' + filename + '\"'
            if (label is None):
                # If we don't set a title, gnuplot automatically
                # uses the data file name, so we set "" as title
                curve_command_string += ' title \"\"'
            else:
                curve_command_string += ' title \"' + label + '\"'
            if (options is not None):
                curve_command_string += ' ' + options
                
            command_list.append(curve_command_string)

        command_string += ', '.join(command_list)

        # For volatile data, add the datastrings          
        if volatile:
            command_string += EOL
            for data in data_list:
                if (self.plot_type == '3D'):
                    command_string += data2string_3d(data[0],
                                                   data[1],
                                                   data[2])

                elif (data[0] is None):
                    command_string += data2string_1d(data[1])

                else:
                    command_string += data2string_2d(data[0],
                                                   data[1])
                        
        # Send command string to gnuplot   
        self._command(command_string)

        return status, message


