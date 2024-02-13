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

""" Parameters used in several modules of this package.

    This module contains some parameters, which are used in several modules
    and subpackages.
"""

# Parameters used for stadard behavior when plotting via matplotlib or gnuplot
#PLOT_INTERFACE  = 'gnuplot'
PLOT_INTERFACE  = 'pylab'
PLOT_SYMBOL     = '+'
PLOT_COLOR      = 'red'
PLOT_SYMBOLS    = [ '.',     'o',    '+',    'x',     '^',    'v',       '>',   '<',      'D'      ]
N_PLOT_SYMBOLS  = len(PLOT_SYMBOLS)
PLOT_COLORS     = [ 'black', 'grey', 'blue', 'green', 'cyan', 'magenta', 'red', 'orange', 'yellow' ]
N_PLOT_COLORS   = len(PLOT_COLORS)
FILL_COLOR      = 'grey'
PLOT_LINE       = 'None'

# Standard answer for error checking
OK	        = 'Ok'

# Standard characters used to separate data rows and columns
EOL		= '\n'	   # End-of-line character
SEP		= '\t'	   # Separator character for data input/output files
ZERO		= 1.0E-9

GPL_MESSAGE     = ( '\n'
                    + 'COPYRIGHT (c) 2020-2024 Pietro Mandracci\n'
                    + '\n'
                    + 'This program is free software: you can redistribute it and/or modify\n'
                    + 'it under the terms of the GNU General Public License as published by\n'
                    + 'the Free Software Foundation, either version 3 of the License, or\n'
                    + '(at your option) any later version.\n'
                    + '\n'
                    + 'This program is distributed in the hope that it will be useful,\n'
                    + 'but WITHOUT ANY WARRANTY; without even the implied warranty of\n'
                    + 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n'
                    + 'GNU General Public License for more details.\n'
                    + '\n'
                    + 'You should have received a copy of the GNU General Public License\n'
                    + 'along with this program.  If not, see <http://www.gnu.org/licenses/>.\n'
                    + '\n' )
