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

""" PYthon tools for SImulation and CAlculus: probability density functions.

    This module contains some probability distribution functions.   

    Documentation is also available in the docstrings.
"""

#import math
from numpy import sqrt, exp, pi

from ..constants import *

#+------------------------------------------+
#| Probability distribution functions (pdf) |
#+------------------------------------------+

def pdf_uniform(x, lower, upper):
    """ Uniform probability distribution function

        f(x) = 1 / (upper-lower)

        Parameters
        ----------
 
        x:      variable
        low:    lower possible value of x
        upp:    upper possible value of x

        Returns
        -------

        pdf_uniform: value of the pdf
    """
        
    return 1.0 / ( upper - lower  )


def pdf_normal(x, mu, sigma):
    """ Normal probability distribution function

        f(x) = 1 / (sqrt(2*pi)*sigma) * exp(- (x-mu)**2 / (2*sigma**2) )

        Parameters
        ----------
 
        x:      variable
        mu:     average
        sigma:  standard deviation

        Returns
        -------

        pdf_normal: value of the pdf
    """
        
    return (1.0 / ( sqrt(2.0*pi) * sigma )
            * exp(- (x - mu)**2 / (2.0*sigma*sigma)))


def pdf_maxwell_energy(e, kt):
    """ Maxwell probability distribution function for the particle energy

        f(e) = 2/sqrt(pi) * (kt)**(3/2) * sqrt(e) * exp(-e/kt)

        Parameters
        ----------
 
        e:      energy
        kt:     temperature in energy units 
                (must be the same units used for e)

        Returns
        -------

        pdf_maxwell_energy: value of pdf
    """

    return 2.0 * sqrt( e / (pi*kt*kt*kt) ) * exp(-e/kt)


def pdf_druyvesteyn_energy_old(e, A, B):
    """ Druyvesteyn probability distribution function for the particle energy

        f(e) = A * sqrt(e) * exp(-B*e)

        Parameters
        ----------
 
        e:      energy
        A:      parameter
        B:      parameter

        Returns
        -------

        pdf_druyvesteyn_energy_old: value of pdf
    """

    return A * sqrt(e) * exp(-B*e*e)


def pdf_druyvesteyn_energy(e, B, M, E, l):
    """ Druyvesteyn probability distribution function for 
        the electron energy
        taken from: 
        FRIEDMAN A.A, KENNEDY L.A, "PLasma Physics and Engineering", 
        Taylor and Francis, 2004, pag. 184

        f(e) = B * exp(- 3*m_e/M * ( e / (q_e*E*l) )**2 )
 
        Parameters
        ----------
 
        B:      normalization coefficient
        e:      electron energy
        M:      background gas molecules mass
        E:      electric field intensity
        l:      mean free path        

        Returns
        -------

        pdf_druyvesteyn_energy: value of pdf
    """

    return (B * sqrt(e)
            * exp(- 3*ELECTRON_MASS / M * ( e / (ELECTRON_CHARGE*E*l) )**2 ))
        

def pdf_margenau_energy(e, B, M, E, w, l):
    """ Margenau probability distribution function for the electron energy
        taken from: 
        FRIEDMAN A.A, KENNEDY L.A, "PLasma Physics and Engineering", 
        Taylor and Francis, 2004, pag. 184

        f(e) = B * exp(- 3*m_e/M * ( (e**2+e*m*w**2*l**2) / (q_e*E*l)**2 )
 
        Parameters
        ----------
 
        B:      normalization coefficient
        e:      electron energy
        M:      background gas molecules mass
        E:      electric field intensity
        w:      electric field pulsation
        l:      mean free path        

        Returns
        -------

        pdf_margenau_energy: value of the pdf
    """

    return (B * sqrt(e)
            * exp(- 3*ELECTRON_MASS / M
                  * ((e**2+ELECTRON_CHARGE*ELECTRON_MASS*w*w*l*l)
                     / (ELECTRON_CHARGE*E*l)**2)))


def pdf_maxwell_speed(v, m, kt):
    """ Maxwell probability distribution function for the particle speed

        f(v) = sqrt(2/pi * (m/kt)**3) * v**2 * exp(-m*v**2 / (2*kt))

        Parameters
        ----------
 
        v:   speed
        m:   mass of the particle
        kt:  temperature in energy units 
             (must be compatible to the ones used for v and m)

        Returns
        -------

        pdf_maxwell_speed: value of the pdf
    """

    return sqrt( 2/pi * (m/kt)**3 ) * v*v * exp(- m*v*v / (2*kt) )


def pdf_maxwell_reduced(x):
    """ Maxwell probability distribution function for the adimensional particle speed

        f(x) = 4/sqrt(pi) * x**2 * exp(-x**2)

        where x = 2/sqrt(pi) v / v_mean

        Parameters
        ----------
 
        x:     adimensional speed x = 2/sqrt(pi) v / v_mean
        Returns
        -------

        pdf_maxwell_reduced: value of the pdf
    """
    
    return 4 / sqrt(pi) * x*x * exp(- x*x )


def pdf_maxwell_velocity_component(vx, m, kt):
    """ Maxwell probability distribution function for a component of the particle velocity vector

        f(vx) = (m / (2*pi*kt))**1/2 * exp(-m*v**2 / (2*kt))

        Parameters
        ----------
 
        vx:  component of the velocity vector
        m:   mass of the particle
        kt:  temperature in energy units 
             (must be compatible to the ones used for v and m)

        Returns
        -------

        pdf_maxwell_speed: value of the pdf
    """

    return sqrt(m / (2*pi*kt)) * exp(- m*vx*vx / (2*kt) )
