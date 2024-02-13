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

""" This module defines the global variables. """

from os import path
from subprocess import run, PIPE

# Check if the gnuplot program exists
# and set the gnuplot_installed flag accordingly
# If the 'which' program is not found, set to None
try: 
    gnuplot_installed = not run(['which', 'gnuplot'],
                                           stdout=PIPE,
                                           stderr=PIPE).returncode
except FileNotFoundError:    
    gnuplot_installed = None

# Allowed gnuplot terminal types and default one
DEFAULT_TERM       = 'x11'
TERMINALS          = {'x11', 'wxt'}
DEFAULT_PRINT_TERM = 'png'
PRINT_TERMINALS    = {'png', 'jpeg', 'eps', 'gif', 'svg',
                      'latex', 'postscript','pdfcairo',
                      'dumb'}
PRINT_EXT          = { 'png':        '.png',
                       'jpeg':       '.jpg',
                       'eps':        '.eps',
                       'gif':        '.gif',
                       'svg':        '.svg',
                       'latex':      '.latex',
                       'postscript': '.ps',
                       'pdfcairo':   '.pdf',
                       'dumb':       '.txt'
                      }

# Default dimensions and position on the screen of new plot windows
DEFAULT_WIDTH      = 800
DEFAULT_HEIGHT     = 600
DEFAULT_XPOS       = 480
DEFAULT_YPOS       = 100

# Plot window types and default one
DEFAULT_TYPE       = '2D'
PLOT_TYPES         = {'2D', '3D'}

# Plotting styles and default one
DEFAULT_STYLE      = 'points'
ALLOWED_STYLES_2D  = {'points', 'dots', 'lines', 'linespoints','impulses', 
                      'steps', 'histeps', 'fsteps', 'boxes', 'histograms',
                      'boxplot'}
ALLOWED_STYLES_3D  = {'points', 'dots', 'lines', 'linespoints'}

# Separator char and end-of-line char
SEP                = '\t'
EOL                = '\n'

# Name of the output directory
DIRNAME            = 'gnuplot.out'

# Default name of print file
PRINT_FILENAME     = 'output_window#'

# Names of subdirs for data and gnuplot output
DIRNAME_DATA       = path.join(DIRNAME, 'data')
DIRNAME_OUT        = path.join(DIRNAME, 'output')
DIRNAME_ERR        = path.join(DIRNAME, 'errors')

# Strings used to construct filenames
FILENAME_VOLATILE  = '-'
FILENAME_DATA      = 'gnuplot_data'
FILENAME_DATA_EXT  = '.csv'
FILENAME_OUT       = 'gnuplot_out'
FILENAME_ERR       = 'gnuplot_err'

# Characters that are not allowed in filenames
INVALID_CHARS      = {'\000', '/', '\\', '*', '!', '|','\"', '\'',
                      '>', '<', ',', '\n', '\t'}
# The invalid characters are substituted with this one
SUBSTITUTE_CHAR    = '_'

# Default values of the purge options
# If True, old data files are deleted by default each time new data is plotted
# and when a window is closed
PURGE_DATA         = True
# If True, the directories are removed (if empty) when
# all the plots are closed
PURGE_DIR          = True

# Time  (in seconds) allowed for the gnuplot process
# to quit nicely, before being killed
TIMEOUT_QUIT       = 1

# REDIRECT_OUT=True:  save gnuplot otuput and errors to files
# REDIRECT_OUT=False: leave it to /dev/stdout and /dev/stderr
# REDIRECT_OUT=None:  send them to /dev/null
REDIRECT_OUT       = False

# Default persistence: if True, plots remain visible after gnuplot is closed
PERSISTENCE        = False

# Defines the default way to pass data to gnuplot:
# True:  pass data inline as a string
# False: write data to data files
VOLATILE           = False

# List of all the active plots, as instances of the _PlotWindow class
window_list        = []
