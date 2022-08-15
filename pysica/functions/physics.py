# COPYRIGHT 2020 by Pietro Mandracci

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
    PYthon tools for SImulation and CAlculus: physical functions 
    
    This module contains some functions used to calculate or to convert physical quantities.

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import numpy

from ..constants import *


# +--------------------+
# | Physical functions |
# +--------------------+

def mean_speed_maxwell(temperature, mass):
    """ Calculates the mean speed of the molecules in a gas at thermal equilibrium 

        Parameters
        ----------

        temperature:    gas temperature / K
        mass:           molecular mass / atomic units

        Returns
        -------

        mean speed / m s**-1

    """

    return numpy.sqrt(8 * K_BOLTZMANN * temperature / (numpy.pi * mass * ATOMIC_UNIT_MASS))


def number_density(pressure, temperature):
    """ Calculates the number density of a gas at given pressure and temperature.

        Parameters
        ----------

        pressure:       gas pressure / Pa
        temperature:    gas temperature / K

        Returns
        -------

        number density: numer density of gas molecules / m**-3
    """

    return pressure / (K_BOLTZMANN * temperature)


def pressure_conversion(pressure_in, unit_in, unit_out):
    """ Converts a pressure measure between different units.

        Parameters
        ----------
                
        pressure_in:    pressure value
        unit_in:        units in which the input pressure is expressed 
                        ('Pa', 'mbar', 'mtorr', 'torr', 'atm', 'psi')
        unit_out:       units in which the output pressure must be expressed 
                        ('Pa', 'mbar', 'mtorr', 'torr')

        Returns
        -------

        pressure_out:   pressure value expressed in converted units
                        if the string given as input or output does not match any of the 
                        pressure units, the value -1 will be returned
    """
                
    if (unit_out=='mtorr'):
        if   (unit_in=='Pa'):    return pressure_in * ATM_PRESSURE_TORR * 1000 / ATM_PRESSURE_PA
        elif (unit_in=='mbar'):  return pressure_in * ATM_PRESSURE_TORR * 1000 / ATM_PRESSURE_PA * 100
        elif (unit_in=='mtorr'): return pressure_in
        elif (unit_in=='torr'):  return pressure_in * 1000
        elif (unit_in=='atm'):   return pressure_in * ATM_PRESSURE_TORR * 1000
        elif (unit_in=='psi'):   return pressure_in * ATM_PRESSURE_TORR * 1000 / ATM_PRESSURE_PSI
        else:                    return -1 


    elif (unit_out=='torr'):
        if   (unit_in=='Pa'):    return pressure_in * ATM_PRESSURE_TORR / ATM_PRESSURE_PA
        elif (unit_in=='mbar'):  return pressure_in * ATM_PRESSURE_TORR / ATM_PRESSURE_PA * 100
        elif (unit_in=='mtorr'): return pressure_in / 1000.0
        elif (unit_in=='torr'):  return pressure_in 
        elif (unit_in=='atm'):   return pressure_in * ATM_PRESSURE_TORR 
        elif (unit_in=='psi'):   return pressure_in * ATM_PRESSURE_TORR / ATM_PRESSURE_PSI
        else:                    return -1        

    elif (unit_out=='Pa'):
        if   (unit_in=='Pa'):    return pressure_in
        elif (unit_in=='mbar'):  return pressure_in * 100.0
        elif (unit_in=='mtorr'): return pressure_in * ATM_PRESSURE_PA   / (ATM_PRESSURE_TORR * 1000)
        elif (unit_in=='torr'):  return pressure_in * ATM_PRESSURE_PA / ATM_PRESSURE_TORR
        elif (unit_in=='atm'):   return pressure_in * ATM_PRESSURE_PA
        elif (unit_in=='psi'):   return pressure_in * ATM_PRESSURE_PA / ATM_PRESSURE_PSI
        else:                    return -1

    elif (unit_out=='mbar'):
        if   (unit_in=='Pa'):    return pressure_in / 100.0
        elif (unit_in=='mbar'):  return pressure_in
        elif (unit_in=='mtorr'): return (pressure_in * ATM_PRESSURE_PA / 100
                                         / (ATM_PRESSURE_TORR * 1000))
        elif (unit_in=='torr'):  return pressure_in * ATM_PRESSURE_PA / 100 / ATM_PRESSURE_TORR
        elif (unit_in=='atm'):   return pressure_in * ATM_PRESSURE_PA / 100
        elif (unit_in=='psi'):   return pressure_in * ATM_PRESSURE_PA / 100 / ATM_PRESSURE_PSI
        else:                    return -1

    elif (unit_out=='atm'):
        if   (unit_in=='Pa'):    return pressure_in * 1.0 / ATM_PRESSURE_PA
        elif (unit_in=='mbar'):  return pressure_in * 1.0 / ATM_PRESSURE_PA * 100
        elif (unit_in=='mtorr'): return pressure_in * 1.0 / ATM_PRESSURE_TORR / 1000
        elif (unit_in=='torr'):  return pressure_in * 1.0 / ATM_PRESSURE_TORR
        elif (unit_in=='atm'):   return pressure_in
        elif (unit_in=='psi'):   return pressure_in * 1.0 / ATM_PRESSURE_PSI 
        else:                    return -1 

    elif (unit_out=='psi'):
        if   (unit_in=='Pa'):    return pressure_in * ATM_PRESSURE_PSI / ATM_PRESSURE_PA
        elif (unit_in=='mbar'):  return pressure_in * ATM_PRESSURE_PSI / ATM_PRESSURE_PA * 100
        elif (unit_in=='mtorr'): return pressure_in * ATM_PRESSURE_PSI / (ATM_PRESSURE_TORR * 1000)
        elif (unit_in=='torr'):  return pressure_in * ATM_PRESSURE_PSI / ATM_PRESSURE_TORR
        elif (unit_in=='atm'):   return pressure_in * ATM_PRESSURE_PSI
        elif (unit_in=='psi'):   return pressure_in
        else:                    return -1 

    else:
        return -1

