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

""" PYthon tools for SImulation and CAlculus: random generation utilities.

    This module contains some functions to generate random number with specific distributions.

    Documentation is also available in the docstrings.
"""

#import math
import numpy
from random import random

from ..constants import *
from .pdf import *

def random_maxwell_velocity(v_mean, modulus=False):
    """ Generate the x,y,z components of a velocity vector following Maxwell distribution

        Parameters
        ----------

        v_mean:   mean value of the velocity module
        modulus:  if True, the modulus of the velocity vector is returned also

        Returns
        -------

        if modulus==False: (vx, vy, vz):   components of the random velocity vector
        if modulus==False: (vx, vy, vz,v): components and modulus of the random velocity vector
    """

    v = random_maxwell_speed(v_mean)

    costheta = 2 * random() - 1
    sintheta = numpy.sqrt(1 - costheta * costheta)
    phi      = 2 * numpy.pi * random()
    vx       = v * sintheta * numpy.cos(phi)
    vy       = v * sintheta * numpy.sin(phi)
    vz       = v * costheta

    if modulus: return (vx, vy, vz, v)
    else:       return (vx, vy, vz)      

    
def random_maxwell_speed(v_mean):
    """ Generate a random speed according to the Maxwell distribution.

        Parameters
        ----------
        
        v_mean: mean value of the velocity module       

        Returns
        -------
        
        random speed value
    """
    
    xmax = 4.0
    ymax = 4.0 / (numpy.sqrt(numpy.pi) * numpy.e)
    f = 0
    y = 1
    while(f < y):
        x = random() * xmax
        y = random() * ymax
        f = pdf_maxwell_reduced(x)

    return sqrt(pi) / 2 * v_mean * x
