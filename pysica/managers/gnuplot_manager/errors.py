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

""" This module defines the error messages used in this package. """

OK = 'Ok'

NOERROR = (0, OK)

# Errors on plot windows
ERROR_NOT_A_PLOT      = (10, 'object is not a plot window')
ERROR_CLOSED_PLOT     = (11, 'trying to operate on a closed plot window')
ERROR_NO_PLOTS        = (12, 'no open plots')
ERROR_UNKNOWN_TERM    = (13, 'unknown terminal type')
ERROR_UNKNOWN_TYPE    = (14, 'unknown plot type')

# Errors about plotting functions or data
ERROR_NO_REPLOT       = (20, 'no functions or curves to replot')
ERROR_NO_FUNCTION     = (21, 'no functions to plot')
ERROR_NO_DATA         = (22, 'no data to plot')
ERROR_INVALID_STYLE   = (23, 'invalid plot style')
ERROR_WRONG_TYPE      = (24, 'invalid data type')
ERROR_REPLOT_VOLATILE = (25, 'replot impossible with volatile curves')
ERROR_LOGSET_VOLATILE = (26, 'logscale toggle impossible with volatile curves')

# Errors about parameters inconsistency (e.g. xmin>xmax)
ERROR_PLOT_PARAMETERS = (30, 'parameters error')

# Errors about file operations
ERROR_FILE            = (40, 'file error')





