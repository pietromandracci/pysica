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

""" Tools for the simulation of plasma discharges: fortran module

    This sub-package contains some modules, compiled from Fortran using f2py, which can be used for simulation 
    of the motion of particles (mainly electrons) in cold plasmas using leap-frog and PIC schemes

    * modules

            - fmodule		-> functions used to simulate the motion of charged particles in a cold plasma discharge (single core)
            - fmodule_parallel	-> functions used to simulate the motion of charged particles in a cold plasma discharge (milti core)

    Documentation is also available in the docstrings.
"""
