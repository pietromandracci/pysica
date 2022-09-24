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

""" Default values of several global variables """

from math import pi

# +-------------------------------+
# | Configuration file parameters |
# +-------------------------------+

class ConfigurationOptions:
    """ Class used to append configuration variables """

    def __init__(self):
        # Default values of parameters that can be modified by user via configuration file
        self.Nmax_particles   = 3000        # Max number of electrons allowed during simulation
        self.rescale_factor   = 10.0        # Factor for rescaling of particle numbers and weights
        self.sim_duration     = 1.0E-6      # Requested total simulation duration [s]
        self.save_delay       = 1           # Data save periodicity (0 = never, 1 = every cycle, 2 = every 2 cycles, etc..) 
        self.dt_output        = 1.0E-9      # Time between data outptut [s]
        self.dt               = 0.0         # Timestep for simulation [s] (must be lower than dt_output):
                                            #        if zero, dt is variable and calculated by the program
        self.dt_var           = True        # (boolean) Decide if a variable timestep should be used
        self.isactive_recomb  = False       # (boolean) Decide if electron/ion recombination processes must be activated
        self.lateral_loss     = False       # (boolean) Decide if electrons and ions will be lost when they reach the border
        self.maxcollfreq      = 0.01        # Maximum allowed collision frequency: used for the calculus of dt
                                            #         in variable dt mode
        self.min_scattered    = 100         # Minimum number of scattering (including null) processes
                                            #        required to apply many particles method
        self.T_neutrals       = 300         # Temperature of neutrals [K]
        self.p_neutrals       = 10          # Total pressure of gases [Pa]
        self.start_ion_deg    = 1.0E-17     # Starting ionization degree (before plasma discharge is lit)
        self.e_min_sigma      = 0.1         # Minimum energy for which cross sections for electron-neutral collisions
                                            #        will be interpolated [eV]
        self.e_max_sigma      = 1.0E3       # Maximum energy for which cross sections for electron-neutral collisions
                                            #        will be interpolated [eV]
        self.e_min_sigma_ions = 0.1         # Minimum energy for which cross sections for ion-neutral collisions
                                            #         will be interpolated [eV]
        self.e_max_sigma_ions = 1.0E2       # Maximum energy for which cross sections for ion-neutral collisions
                                            #        will be interpolated [eV]
        self.min_ion_density  = 0.0         # Minimum ion density required to activate electron/ion recombination processes
        self.N_sigma          = 1000        # Number of electron cross-section values to be calculated by interpolation
        self.N_sigma_ions     = 100         # Number of ion cross-section values to be calculated by interpolation
        self.distance         = 0.01        # Distance between electrodes [m]
        self.length           = 0.1         # Lateral length of square electrodes [m]
        self.V_bias           = 100         # Electric bias applied to electrodes [V]
        self.frequency        = 0.0         # Electric field frequency [Hz] (zero means DC)
        self.phase            = 90          # Electric field phase at t=0 [deg] (if frequency==0, this is ignored)
        self.N_cells          = 100         # Number of cells used in PIC scheme
        self.dot_points       = False       # (boolean) Decide if points must be used in plots instead of crosses
        self.n_max_points     = 1000        # Maximum number of points allowed in hisotric plots
        self.decimation_factor= 10          # When max number of points is reached, a point is removed every this number
            

# +---------------------------------+
# | Command line options parameters |
# +---------------------------------+

# Default values of parameters that can be modified via command-line options
verbosity               = 1      # How much should the program bore with text-based output
plot_xsec               = False  # Plot cross sections at start
debug_level_python      = 0      # For debugging purposes
debug_level_fortran     = 0      # 
print_only              = False  # Only print information; do not start iterations
save_defaults           = False  # Save default values to a configuration file and exit
gui_mode                = False  # Start the GUI
batch_mode              = False  # Suppress any input from user
redirect_output         = False  # Save output to a log file instead of on showing it stdout
redirect_error_messages = False  # Save python error messages to a log file instead of showing them on stderr
text_window_width       = 160    # Width of the output text window in GUI
text_window_height      = 39     # Hight of the output text window in GUI
text_window_font_size   = 12     # Size of font used in the output text window in GUI


            
# +----------------------+
# | Immutable parameters |
# +----------------------+

# Directories, filemames and separators
EXEC_PATH                 = '/usr/bin/'
EDITOR_NAME               = 'mousepad'
#EDITOR_NAME               = 'gedit'
BROWSER_NAME              = 'firefox'
DOC_URL                   = 'https://github.com/pietromandracci/pysica/tree/master/doc/plasmapro/ccpla_manual.rst'
SCRIPTNAME                = 'ccpla'                       # Name of the script
EXT                       = '.csv'                        # Extension to use for saved data files
EOL                       = '\n'                          # End-of-line char for data files
SEP                       = '\t'
SAVE_DIRECTORY_NAME       = SCRIPTNAME + '.out'           # Directory into which simulation data is saved
FILENAME_DEFAULTS         = SCRIPTNAME + '.defaults'      # Name of the file to which default parameters are saved with -s option
FILENAME_CONFIG           = SCRIPTNAME + '.conf'          # Configuration file 
FILENAME_NEUTRALS         = SCRIPTNAME + '.neutrals'      # Filename for neutrals properties
FILENAME_OUTPUT_LOG       = SCRIPTNAME + '_output.log'    # File to which output will be saved
FILENAME_ERROR_LOG        = SCRIPTNAME + '_errors.log'    # File to which python error messages will be saved
NAME_STAT_ELE             = '_means_ele'                  # Name of output file: time evolution of electron parameters
NAME_STAT_ION             = '_means_ion'                  # First part of name of output file: time evolution of ion parameters
NAME_STAT_NEU             = '_means_neu'                  # Name of output file: time evolution of several parameters
NAME_DISTRIB_ELE          = '_eedf'                       # Name of output file: electrons energy distribution
NAME_I                    = '_current'                    # Name of output file: electric current
NAME_V                    = '_potential'                  # Name of output file: electric potential spatial distribution
NAME_DISTRIB_ION          = '_iedf'                       # First part of name of output files: ion energy distribution
NAME_Z_ELECTRONS          = '_zpos_ele'                   # Name of output file: electron positions, y component
NAME_Z_IONS               = '_zpos_ion'                   # First part of name of output files: ion positions, z component 
NAME_ION_ELASTIC          = '_sigma4ions_Ar_elastic'
NAME_ION_CHARGE_EX        = '_sigma4ions_Ar_charge-exchange'

# Simulation parameters
NMAXPARTICLES             = 100000                        # Maximum allowed number of computational particles
MAX_RESCALE_FACTOR        = 1.0E6                         # Maximum allowed rescale factor
DEFAULT_RESCALE_FACTOR    = 10.0                          # Default scaling factor for rescaling of particles weight
DEFAULT_MIN_SCATTERED     = 100                           # Default minimum number of scattering (including null) processes
                                                          #        required to apply many particles method
START_WEIGHT              = 1.0                           # Initial value of computational weight for electrons and ions
NMAXCELLS                 = 1000                          # Maximum allowed number of cells in PIC scheme
NMINCELLS                 = 10                            # Minimum allowed number of cells in PIC scheme
MINTRATIO                 = 100                           # Minimum allowed ratio between electric field period 
                                                          #       and simulation timstep (dt)
MAX_IONIZATION_DEGREE     = 1.0E-4                        # Maximum allowed initial value of ionization degree

MAX_DISSOCIATION_TYPES    = 3                             # Maximum allowed number of dissociation channels for each neutral type
MAX_EXCITATION_TYPES      = 3                             # Maximum allowed number of excitation channels for each neutral type
MOLECULE_TYPES            = ('a', 'm', 'p')               # Types of molecules
INITIAL_VELOCITY          = 0.0                           # / m*s**-1 starting speed of electrons
INITIAL_ANGLE             = pi/6                          # / rad     starting angle between electrons velocity and z-axis direction

# GUI parameters
TEXT_WINDOW_WIDTH_MIN     = 120                           # Width of the output text window: minumum value
TEXT_WINDOW_WIDTH_MAX     = 200                           # Width of the output text window: maximum value
TEXT_WINDOW_HEIGHT_MIN    = 20                            # Height of the output text window: minumum value
TEXT_WINDOW_HEIGHT_MAX    = 80                            # Height of the output text window: maximum value
TEXT_WINDOW_FONT_TYPE     = 'mono'                        # Type of font used in the output text window
TEXT_WINDOW_FONT_SIZE_MIN = 6                             # Size of font used in the output text window
TEXT_WINDOW_FONT_SIZE_MAX = 18                            # Size of font used in the output text window
TEXT_WINDOW_FORE          = 'black'                       # Foregroun color of text in the output text window
TEXT_WINDOW_BACK          = 'light grey'                  # Background color in the output text window

FONT_MESSAGE              = ('monospace','12')
FONT_MESSAGE_SMALL        = ('monospace','10')
FONT_MESSAGE_BIG          = ('monospace','14')

IDLE_TIME                 = 1 #ms                         # Idle time betweem each main iteration cycle (in ms)
                                                          #     this is necessary to allow responsivity of the GUI
IDLE_TIME_START           = 100 #ms                       # Time between activation of main and about windows at program start
IDLE_TIME_SHUTDOWN        = 100 #ms                       # Time allowed to quit the kernele beefore killing it
IDLE_TIME_GRAPHS          = 500 #ms                       # Time allowed to reset plots before lowering them
MIN_DT_OUTPUT_EXP         = -15                           # log10 of the minimum allowed value of time between text output 
MAX_DT_OUTPUT_EXP         = -6                            # log10 of the maximum allowed value of time between text output 
MIN_PLOT_DELAY            = 1                             # Min value of number of text output cycles between each update of plots
MAX_PLOT_DELAY            = 100                           # Max value of number of text output cycles between each update of plots
RES_PLOT_DELAY            = 10                            # Minimum amount of which the value is modified in the ruler
DEF_PLOT_DELAY            = 1                             # Default value of number of text outputs between each update of plots

# Plot parameters
N_MAX_OUTPUT              = 1000                          # Maximum number of data values registered for historic plots
MAX_DECIMATION_FACTOR     = 20                            # Maximum allowed value of the decimation factor
N_BINS_LOG                = 100                           # Number of bins over which x log scale is activated in EEDF and IEDF plots
DEFAULT_TERMINAL          = 'x11'                         # Gnuplot terminal type
PERSIST                   = 0                             # Persistency of gnuplot graphs, after end of script (1=Yes, 0=No)
SCREEEN_WIDTH             = 1024                          # Default screen width
SCREEEN_HEIGHT            = 768                           # Default screen height
VERTICAL_SHIFT            = 25                            # Vertical shift toward top of screen
HBORDER                   = 12                            # Horizontal border of plot windows
VBORDER                   = 50                            # Vertical border of plot windows (including window header)
LABEL_XPOS                = 1                             # x position of the label with time value
LABEL_YPOS                = 1                             # y position of the label with time value
DEL_DATA_FILES_DELAY      = 1                             # Time to wait (in seconds) before deleting data files
                                                          # when closing cross section plots
# Debug parameters
KERNEL_DEBUG_LEV          = 1                             # Minimum debug level to get kernel debug messages
