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
    Tools for the simulation of plasma discharges.

    This package contains some modules and subpackages for the simulation of low pressure CCP discharges.

    Sub-packages
    ------------

    *discharge*
        modules used to define particles ensambles (charged and neutrals), reactors properties
        and to simulate the charged particles motion


    Modules
    -------

    *ccpla_defaults*
        some constants used in the package
    *ccpla_init*
        intialization functions
    *ccpla_print*
        functions used to print simulation info in text windows
    *ccpla_gui*
        functions and classes used to run the simulaion within a GUI
    *cppla_kernel*
        numerical kernel which is called from the ccpla_gui module and runs in parallel to it
    *ccpla_plot*
        functions used to plot simulation data

    Scripts
    -------

    *ccpla*
        script designed to simulate the motion of charged particles (electrons and ions) in 
        a capacitivley coupled plasma discharge
    *ccpla_analysis*
        analysis of output from ccpla script 
        (ccpla must have been run with the option -s to save simulation results on ascii files)

    Documentation is also available in the docstrings.
"""

