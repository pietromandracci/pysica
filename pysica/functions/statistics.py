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

"""
    PYthon tools for SImulation and CAlculus: statistical functions.

    This module contains some statistical functions.

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import numpy

# +----------------------+
# | Statistics functions |
# +----------------------+
          

def calculate_bins_number(n_data, method='sqrt', width=None, stdv=None):
    """ Calculates the number of bins to create an histogram

        Parameters
        ----------

        n_data:      number of data
        method:      method used to calculate the  number of intervals
                        'sqrt':  n_bins = sqrt(n_points)
                        'log2':  n_bins = 1 + log2(n_points)
                        'root3': n_bins = n_points**(1/3)
                        'norm':  n_bins = n_data / (3.5 * stdv / n_points**(1/3) )
        with:        width of the data interval (i.e. maximum value - minimum value)
                     used only if the 'norm' method is requested
        stdv:        standard deviation of the data
                     used only if the 'norm' method is requested

        Returns
        -------
        number of intervals, on error None is returned
    """

    if   (method == 'sqrt'):  return int( numpy.sqrt(n_data) )
    elif (method == 'log2'):  return int( 1 + numpy.log2(n_data) )
    elif (method == 'root3'): return int(2.0 * pow(n_data, 1.0/3.0))
    elif (method == 'norm'):
        if ( (width is None) or (stdv is None) ): return None
        else:
            delta  = 3.5 * stdv / pow(n_data, 1.0/3.0)
            return int(width / delta)
    else: return None
