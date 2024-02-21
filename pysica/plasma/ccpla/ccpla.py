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

""" Simulation of a capacitively coupled plasma discharge. """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

# Modules from the standard Python library
import sys
import math
from random import random
from optparse import OptionParser
import tkinter as tk
#import tkMessageBox

# Modules provided by the Python community
import numpy

# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.plasma.ccpla.ccpla_defaults import *

# Import required modules, classes, and functions
from pysica.managers.io.io_screen import wait_input, clear_screen
from pysica.managers import unit_manager
from pysica.plasma.ccpla import ccpla_gui
from pysica.plasma.ccpla.discharge import reactors, target_particles, moving_particles, particle_mover
from pysica.plasma.ccpla.ccpla_init import *
from pysica.plasma.ccpla.ccpla_print import *

# The gnuplot_installed flag is set to False if gnuplot is not installed
from pysica.managers.gnuplot_manager import gnuplot_installed


# +---------------------------------------------------+
# | Manage how to write messages and exit the program |
# +---------------------------------------------------+

class CcplaMessageWindow(tk.Frame):
    def __init__(self, message, error, master=None):
        tk.Frame.__init__(self, master)

        if error:
            self.back         = 'light yellow'
            self.fore         = 'red'
            self.button_label = 'Exit'
        else:
            self.back         = 'light grey'
            self.fore         = 'blue'
            self.button_label = 'Dismiss'
        self.message = tk.Message(self,
                                  text=message,
                                  background=self.back,
                                  foreground=self.fore,
                                  width=1500,
                                  justify=tk.CENTER,
                                  relief=tk.RIDGE)
        self.message.grid()
        self.button  = tk.Button(self, text=self.button_label, command=self.master.destroy)
        self.button.grid()
        self.grid()
 
def create_message_window(message, error):
    root = tk.Tk()
    if error: root.title("Ccpla Error")
    else:     root.title("Ccpla Message")
    root.resizable(False, False)

    message_window = CcplaMessageWindow(message, error, master=root)

    root.mainloop()

def exit_ccpla(message=None, gui=False, error=False):
    if (message is None): message = "Exiting program " + SCRIPTNAME
    if not gui:
        sys.exit(EOL + message + EOL)
    else:
        create_message_window(message, error)
        sys.exit()

def print_message(message, gui=False, error=False):
    if gui:
        create_message_window(message, error)
    else:
        if error: message = 'ERROR: ' + message
        print('\n' + message + '\n')
             
# +------------------------+
# | Command line arguments |
# +------------------------+

usage = "usage: %prog [options]"        # String to be given as program usage description
parser = OptionParser(usage)

# Use parallel Fortran module
help_string = "use multicore parallel Fortran module"
parser.add_option("-m", "--multicore", action="store_true",  dest="cpu_multicore", default=False, 
                  help=help_string)

# Do not start simulation, only print parameters values and exit
help_string = "show parameters values and exit"
parser.add_option("-p", "--print-only", action="store_true",  dest="print_only", default=False, help=help_string)

# Only save default values 
help_string = "write default parameters to file \"" + FILENAME_DEFAULTS + "\" and exit"
parser.add_option("-s", "--save-defaults", action="store_true",  dest="save_defaults", default=False, 
                  help=help_string)

# Batch mode
help_string = "suppress all input from user"
parser.add_option("-b", "--batch-mode", action="store_true", dest="batch_mode", default=batch_mode,
                  help=help_string)

# GUI mode
help_string = "start GUI"
parser.add_option("-g", "--gui-mode", action="store_true", dest="gui_mode", default=gui_mode,
                  help=help_string)

# Width of text window in GUI mode
help_string = ("Width of text window in GUI ["
               + str(TEXT_WINDOW_WIDTH_MIN) + ".." + str(TEXT_WINDOW_WIDTH_MAX)
               + "] (default=" + str(text_window_width) + ")")
parser.add_option("-W", "--text-window-width", action="store", type="int", dest="text_window_width",
                  default=text_window_width, help=help_string)

# Height of text window in GUI mode
help_string = ("Height of text window in GUI ["
               + str(TEXT_WINDOW_HEIGHT_MIN) + ".." + str(TEXT_WINDOW_HEIGHT_MAX)
               + "] (default=" + str(text_window_height) + ")")
parser.add_option("-H", "--text-window-height", action="store", type="int", dest="text_window_height",
                  default=text_window_height, help=help_string)

# Font dimension in text window in GUI mode
help_string = ("Font size in GUI text window ["
               + str(TEXT_WINDOW_FONT_SIZE_MIN) + ".." + str(TEXT_WINDOW_FONT_SIZE_MAX)
               + "] (default=" + str(text_window_font_size) + ")")
parser.add_option("-F", "--text-window-font", action="store", type="int", dest="text_window_font_size",
                  default=text_window_font_size, help=help_string)

# Redirect output
help_string = "redirect screen output to file \'" + FILENAME_OUTPUT_LOG + "\'"
parser.add_option("-o", "--redirect-output", action="store_true", dest="redirect_output", default=redirect_output,
                  help=help_string)        

# Redirect error messages
help_string = "redirect error messages to file \'" + FILENAME_ERROR_LOG + "\'"
parser.add_option("-e", "--redirect-errors", action="store_true", dest="redirect_error_messages",
                  default=redirect_error_messages,
                  help=help_string)

# Verbosity level
help_string = "verbosity level [0..3] (default="+str(verbosity)+")"
parser.add_option("-v", "--verbosity", action="store", type="int", dest="verbosity",  default=verbosity,
                  help=help_string)

# Debug level for python code
help_string = "Python debug level [0..2] (default="+str(debug_level_python)+")"
parser.add_option("-d", "--debug-level-python", action="store", type="int", dest="debug_lev",
                  default=debug_level_python, help=help_string)

# Debug level for Fotran code
help_string = "Fortran debug level [0..3] (default="+str(debug_level_fortran)+")"
parser.add_option("-D", "--debug-level-fortran", action="store", type="int", dest="debug_lev_for",
                  default=debug_level_fortran, help=help_string)        
        
# Plot cross sections before starting simulation
help_string = "plot cross sections graphs before start"
parser.add_option("-x", "--graph-xsec", action="store_true",  dest="plot_xsec", default=False, help=help_string)

(cl_options, args) = parser.parse_args()

if cl_options.redirect_error_messages:
    sys.stderr = open(FILENAME_ERROR_LOG, 'w')
    
if cl_options.redirect_output:
    sys.stdout = open(FILENAME_OUTPUT_LOG, 'w')
    

# +-----------------------------+   
# | Error checking and warnings |
# +-----------------------------+

# This option must change an import directive in the discharge.particle_mover module
# so it cannot be passed as a variable to a function inside the module
# in this way it can be imported from this  module

if (cl_options.verbosity not in list(range(4))):
    exit_ccpla('ERROR: verbosity must be in range 0..3' + EOL,
               gui=cl_options.gui_mode, error=True)
    
if (cl_options.debug_lev not in list(range(3))):
    exit_ccpla('ERROR: Python debug level must be in range 0..2' + EOL,
               gui=cl_options.gui_mode, error=True)
    
if (cl_options.debug_lev_for not in list(range(4))):
    exit_ccpla('ERROR: Fortran debug level must be in range 0..3' + EOL,
               gui=cl_options.gui_mode, error=True)
    
if cl_options.gui_mode:
    if not gnuplot_installed:        
        exit_ccpla('ERROR: gnuplot is not installed on your system\n'
                   + 'GUI mode is not available', gui=True)            
    cl_options.verbosity = 0
    cl_options.batch_mode = True
    
if ( (cl_options.text_window_width < TEXT_WINDOW_WIDTH_MIN)
     or (cl_options.text_window_width > TEXT_WINDOW_WIDTH_MAX) ):
    exit_ccpla('ERROR: text window width must be in range '
               + str(TEXT_WINDOW_WIDTH_MIN)
               + '..'
               + str(TEXT_WINDOW_WIDTH_MAX)
               + EOL,
               gui=cl_options.gui_mode,
               error=True)
    
if ( (cl_options.text_window_height < TEXT_WINDOW_HEIGHT_MIN)
     or (cl_options.text_window_height > TEXT_WINDOW_HEIGHT_MAX) ):
    exit_ccpla('ERROR: text window height must be in range '
               + str(TEXT_WINDOW_HEIGHT_MIN)
               + '..'
               + str(TEXT_WINDOW_HEIGHT_MAX)
               + EOL,
               gui=cl_options.gui_mode,
               error=True)
    
if ( (cl_options.text_window_font_size < TEXT_WINDOW_FONT_SIZE_MIN)
     or (cl_options.text_window_font_size > TEXT_WINDOW_FONT_SIZE_MAX) ):
    exit_ccpla('ERROR: text window font size must be in range '
               + str(TEXT_WINDOW_FONT_SIZE_MIN)
               + '..'
               + str(TEXT_WINDOW_FONT_SIZE_MAX)
               + EOL,
               gui=cl_options.gui_mode,
               error=True)

if ((cl_options.print_only) and (cl_options.gui_mode)):
    exit_ccpla('ERROR: -p option is not compatible with GUI mode', gui=True, error=True)
    
if (cl_options.plot_xsec):
    if not gnuplot_installed:
        exit_ccpla('ERROR: gnuplot is not installed on your system\n'
                   + 'plotting of cross section is not possible',
                   gui=cl_options.gui_mode,
                   error=True)
    elif cl_options.gui_mode:        
        cl_options.plot_xsec = False
        print_message(message='WARNING: -x option is disabled in GUI mode\n'
                      + 'cross section graphs can be accessed from the \"Parameters\" menu entry', 
                      gui=True,
                      error=False)        
    elif cl_options.batch_mode:
        cl_options.plot_xsec = False        
        print_message(message='WARNING: plotting of cross sections is disabled in batch mode', 
                      gui=False,
                      error=False)

cl_options.text_window_font = (TEXT_WINDOW_FONT_TYPE, str(cl_options.text_window_font_size))

# +----------------------+
# | Print copyright info |
# +----------------------+

if (not cl_options.gui_mode): print(GPL_MESSAGE)

if (not cl_options.batch_mode): wait_input()
    

# +-------------------------+
# | Read configuration file |
# +-------------------------+

parameters        = ConfigurationOptions()
(status, message) = initialize_parameters(parameters,
                                          verbose=(cl_options.verbosity > 0),
                                          saveonly=cl_options.save_defaults,
                                          filename_config=FILENAME_CONFIG,
                                          filename_defaults=FILENAME_DEFAULTS)
if (status != 0):
        exit_ccpla(EOL + message + EOL, gui=cl_options.gui_mode, error=True)
if (cl_options.save_defaults):
        exit_ccpla('\nDefault values of parameters were saved to file: \"' + FILENAME_DEFAULTS + '\"\n',
                   gui=cl_options.gui_mode, error=False)

# +------------------------------------+
# | Initialize reactor properties data |
# +------------------------------------+

ccp = reactors.CcpProperties(parameters.distance,
                             parameters.length,
                             parameters.V_bias,
                             parameters.frequency,
                             parameters.phase,
                             parameters.N_cells,
                             parameters.lateral_loss)


# +---------------------+
# | Load gas properties |
# +---------------------+

if (cl_options.verbosity > 0): print('\nCreating neutrals ensamble from file \"' + FILENAME_NEUTRALS + '\" ...')
neutrals = target_particles.TargetParticles( parameters.N_sigma,
                                             parameters.N_sigma_ions,
                                             parameters.T_neutrals,
                                             parameters.p_neutrals,
                                             parameters.min_scattered,
                                             parameters.isactive_recomb,
                                             filename=FILENAME_NEUTRALS )
(status, message) = neutrals.read_error
if (status != 0):
    ERROR = '\nERROR in file \"'+FILENAME_NEUTRALS+'\": '
    exit_ccpla(ERROR + message + EOL, gui=cl_options.gui_mode, error=True)

if (cl_options.verbosity > 0): print('\nReading gases properties from file \"' + FILENAME_NEUTRALS + '\" ...')    
(status, message) = neutrals.read_properties( FILENAME_NEUTRALS, '\t', 
                                              parameters.e_min_sigma,
                                              parameters.e_max_sigma,
                                              parameters.e_min_sigma_ions,
                                              parameters.e_max_sigma_ions, 
                                              debug=(cl_options.debug_lev>1) )
if (status != 0):
    ERROR = '\nERROR in file \"'+FILENAME_NEUTRALS+'\": '
    exit_ccpla(ERROR + message + EOL, gui=cl_options.gui_mode, error=True)

#neutrals.min_scattered = parameters.min_scattered


# +---------------------+
# | Load cross sections |
# +---------------------+

(status, message) = initialize_cross_sections(neutrals, cl_options, parameters.isactive_recomb)
if (status != 0): exit_ccpla(message + EOL, gui=cl_options.gui_mode, error=True)


# +-----------------------------------+
# | Electrons and ions initialization |
# +-----------------------------------+

# Define the ensamble of charged particles moving in the plasma: electrons and ions
if (cl_options.verbosity > 0): print('\nDefining charged particles ensambles ...')
charges = moving_particles.MovingParticles(neutrals.types+1, parameters.Nmax_particles, parameters.start_weight, parameters.rescale_factor)

(status, message) = initialize_ensambles(charges, neutrals, parameters, cl_options)
if (status!=0): exit_ccpla(message + EOL, gui=cl_options.gui_mode, error=True)


# +-----------------------------------+
# | If it was requested, call the GUI |
# +-----------------------------------+

if cl_options.gui_mode:
    message = ccpla_gui.main(charges, neutrals, ccp, parameters, cl_options)
    exit_ccpla()


# +---------------------------------------+
# | Initialization of save directory name |
# +---------------------------------------+

if (parameters.save_delay > 0):
    initialize_filenames(neutrals, parameters)
    create_save_dir(parameters)

       
# +-------------------------------------------------+
# | Print some information on simulation parameters |
# +-------------------------------------------------+

if (cl_options.verbosity > 0):
    print()
    print(print_simulation_information(neutrals, ccp, parameters, cl_options))
    print()
    if (not cl_options.batch_mode): wait_input()        
    print(print_gas_information(neutrals))
    print() 
    if (not cl_options.batch_mode): wait_input()

# If requested by command-line option, exit the script after printing information
if cl_options.print_only: exit_ccpla()


# +--------------------------+
# | Data file initialization |
# +--------------------------+

if (parameters.save_delay > 0):
    charges.initialize_savefiles(parameters.filename_stat_ele,
                                 parameters.filename_stat_ion,
                                 parameters.filename_distrib_ele,
                                 parameters.filename_distrib_ion,
                                 parameters.filename_epos_z,
                                 parameters.filename_ipos_z,
                                 append=False, sep='\t', ext=EXT)
    neutrals.initialize_savefile(parameters.filename_stat_neu, append=False, sep='\t', ext=EXT)
    ccp.initialize_savefiles(parameters.filename_I, parameters.filename_V, append=False, sep='\t', ext=EXT)
    print(print_filenames(charges, neutrals, ccp))
    print() 
    if (not cl_options.batch_mode): wait_input()

        
# +----------------+
# | Main iteration |
# +----------------+

if (cl_options.verbosity > 0): print('\nSTARTING SIMULATION')

# Initialize time and main loop index
#charges.time = 0
#charges.dt   = parameters.dt

# Data save counters and flag
i_save_data = 1
i_save_dist = 1
save_dist   = True

while (charges.n_active(0) > 0):

    if (cl_options.debug_lev > 0): 
        print('t= ' + unit_manager.print_unit(charges.time,'s') + ' V= ' + unit_manager.print_unit(ccp.V,'V',3))
        print('E      min, max, mean, sigma [eV] ', charges.e_min(0), charges.e_max(0), 
                                                    charges.e_average(0), charges.e_sigma(0))
        print('theta  min, max, mean, sigma [deg]', math.degrees(charges.theta_min(0)), 
                                                    math.degrees(charges.theta_max(0)), 
                                                    math.degrees(charges.theta_average(0)),
                                                    math.degrees(charges.theta_sigma(0))) 
        print()

    if (cl_options.verbosity > 2):
        time_before          = charges.time
        n_active_el_before   = charges.n_active(0)
        electric_bias_before = ccp.V

    particle_mover.move_particles(charges, neutrals, ccp, parameters, cl_options)

    if (cl_options.verbosity > 2):
        clear_screen()
        print(print_runtime_info(charges,
                                 neutrals,
                                 ccp,
                                 parameters,
                                 cl_options,
                                 time_before,
                                 electric_bias_before,
                                 n_active_el_before))
        
    # Save data to files, if required
    if (parameters.save_delay > 0):
        # Data are seved every n cycles, where n = parameters.save_delay (if not zero)
        if ( (i_save_data >= parameters.save_delay) and (charges.n_active(0) > 0) ):
            # If required, energy and z distributions are not saved always,
            # but every n data saves only, to preserve disk space
            # where n = parameters.save_delay_dist
            if (parameters.save_delay_dist > 1):
                if (i_save_dist >= parameters.save_delay_dist):
                    save_dist = True
                    i_save_dist = 0
                else:
                    save_dist = False
                i_save_dist += 1
            charges.save_data_to_files(save_edf=save_dist, save_z=save_dist)
            neutrals.save_data_to_files(charges.time)
            ccp.save_data_to_files(charges.time)
            i_save_data = 0
        i_save_data += 1

    # Exit the loop if a maximum duration was provided and it has been reached
    if ( (parameters.sim_duration > 0) and (charges.time >= parameters.sim_duration) ): break
        
if (cl_options.verbosity > 0): 
    if (charges.time < parameters.sim_duration):
        print('\n\nSIMULATION INTERRUPTED (NO MORE ELECTRONS)\n\n')
    else:
        print('\n\nSIMULATION COMPLETED\n\n')
