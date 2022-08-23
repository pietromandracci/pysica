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

""" PYthon tools for SImulation and CAlculus: utility functions

    This subpackage contains some modules useful to manage the reding or writing of raw data or configuration variablees

    This sub-package contains some modules which can be used for managing data input/output from/to files, 
    print of physical quantities whith automatic managing of the units prefixes

    * subpackages
            - gnuplot manager        -> plot graphs by means of gnuplot

    * modules
            - data_manager           -> read/write data from/to ASCII files
            - unit_manager           -> prints physical quantities, managing the unit prefixes 
                                        (e.g. 1.23456E-9 A -> 1.23 nA)

    Documentation is also available in the docstrings.
"""
