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

""" Classes for the simulation of neutral and charged particles in a plasma discharge.

    This module contains some classes defining the characteristics of reactors and the physical quantities involved
    in their usage, such as electric fields, electrodes dimensions and so on.

    Documentation is also available in the docstrings.
"""

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import math
import numpy

# +-------------+
# | CCP Reactor |
# +-------------+

class CcpProperties:
        """"A class defining the physical properties of a capacitively coupled cold plasma reactor."""

        def __init__(self, d, l, V, f, phi, n_cells, lateral_loss):
                """Initialises the properties of a CCP plasma reactor.

                        Parameters
                        ----------

                        d:            distance between the electrodes
                        l:            lateral length the electrodes
                        V:            potential difference applied between the electrodes: peak value
                        f:            frequency of the bias applied between the electrodes
                        phi:          starting phase of the potential difference
                        n_cells:      number of cells used for PIC scheme
                        lateral_loss: if set to True, electrons and ions that reach the x or y borders will be lost,
                                      otherwise will be recovered


                        Initialized data attributes
                        ---------------------------

                        self.distance           distance between the electrodes
                        self.length             lateral length the electrodes
                        self.area               area of the electrodes
                        self.volume             plasma volume
                        self.lateral_loss:      if set to True, electrons and ions that reach the x or y borders will be lost,
                                                otherwise will be recovered
                        self.E_peak             intensity of the average electric between the electrodes: peak value
                        self.V_peak             intensity of the external electric bias applied between the electrodes: peak value
                        self.frequency          frequency of the external electric bias applied between the electrodes
                        self.pulsation          pulsation of the external electric bias applied applied between the electrodes
                        self.period             period of the external electric bias applied applied between the electrodes
                        self.phase              phase of the external electric field applied to the particles
                        self.average_current    average current measured at the electrodes 
                                                it is defined as a rank-zero numpy array in order to get working f2py 
                                                automatic conversion between fortran and python
                        self.potential          array used to store electric potential values at cells boundaries in PIC scheme
                        self.charge_density     array used to store electric charge density values at cells boundaries in PIC scheme
                        self.grid_points        array used to store the posizion along z-axix of the grid nodes in PIC scheme

                """
                
                self.distance           = d
                self.length             = l
                self.area               = l*l
                self.volume             = self.area * d
                self.lateral_loss       = lateral_loss
                self.E_peak             = V / d
                self.V_peak             = V 
                self.V                  = self.V_peak * math.sin(phi)
                self.frequency          = f
                self.pulsation          = 2 * math.pi * f
                if (f > 0):
                        self.period     = 1.0 / f
                else:
                        self.period     = numpy.inf
                self.phase              = phi
                self.average_current    = numpy.zeros(1, 'd')          # Needs to be an array since it will be passed to f2py
                self.potential          = numpy.zeros(n_cells+1, 'd')
                self.charge_density     = numpy.zeros(n_cells+1, 'd')
                self.grid_points        = numpy.zeros(n_cells+1, 'd')
                self.delta_grid         = d / n_cells
                for i in range(n_cells+1):
                        self.grid_points[i] = i * self.delta_grid

