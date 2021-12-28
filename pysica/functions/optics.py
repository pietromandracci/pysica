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

""" PYthon tools for SImulation and CAlculus: calculation of optical absorption in thin films

    This module contains some functions that may be useful for the calculation
    of some optical quantities (such as reflactive index or absorption coefficient) in thin films.

    Documentation is also available in the docstrings.
"""

# +-------------------------+
# | Import required modules |
# +-------------------------+

# Mudules from the standard Python library
import math

# Modules from the Python community
import numpy

#
# All functions, except "alpha_brutal" and "n_brutal", from: R. Swanepoel, J. Phys. E: Sci. Instrum, vol. 16 (1983) pp. 1214-1222
#
# s    -> refractive index of the substrate
# n    -> refractive index of the film
# t    -> film thickness
# TMax -> value of transmittance maxima envelope
# Tmin -> value of transmittance minima envelope
# Ta   -> value of geometric average of transmittance maxima and minima envelope 


# +------------------------------+
# | Refractive index calculation |
# +------------------------------+

def n_transparent_1(s, Tmin):
        M  = 2*s / Tmin  -  (s*s  +  1)  /  2
        return numpy.sqrt( M + numpy.sqrt(M*M - s*s) )

# Valid in transparent region using Ta
def n_transparent_2(s, Ta):
        H  = 4*s*s / ( Ta*Ta * (s*s + 1) ) - (s*s + 1) / 2
        return numpy.sqrt( H + numpy.sqrt(H*H - s*s) )

def n_lowabsorption(s, TMax, Tmin):
        N  = 2*s * (TMax-Tmin) / (TMax*Tmin) + (s*s + 1) / 2
        return numpy.sqrt( N + numpy.sqrt(N*N - s*s) )

# Calculate n from the distance between adjacent maxima/minima of transmission or reflection spectrum
# t  -> film thickness in nm
# dE -> energy separation between Maxima (or minima, since they should be equal) in eV
def n_brutal(dE, t):
        return 1239.84 / (2*dE*t) 

# +------------------------------------+
# | Absorption coefficient calculation |
# +------------------------------------+

# Calculate the absorption coefficient
# t  -> film thickness in nm
# n  -> refrative index of the thin film
# s  -> refrative index of the substrate
# Ta -> optical transmittance (fringes free, i.e. geometric average of maxima and minima envelopes) 
# T0 -> 1-R, where R is the optical reflectance


def alpha_A3(s, n, Ta, t):
        R1 = ( (1-n) / (1+n) )**2
        R2 = ( (n-s) / (n+s) )**2
        R3 = ( (s-1) / (s+1) )**2
        P  = (R1-1) * (R2-1) * (R3-1)   
        
        Q  = 2 * Ta * (R1*R2 + R1*R3 - 2*R1*R2*R3)
        x  = ( P + numpy.sqrt( P*P + 2*Q*Ta*(1-R2*R3) ) ) / Q
        return -numpy.log(x) / t


def alpha_19(s, n, Ta, t):
        G1 = 128*n**4*s*s / Ta**2
        G2 = n*n * (n*n - 1)**2 * (s*s - 1)**2
        G3 = (n*n - 1)**2 * (n*n - s*s)**2
        G  = G1 + G2 + G3       
        
        x0 = numpy.sqrt( G*G - (n*n - 1)**6 * (n*n-s**4)**2 )
        x  = numpy.sqrt(G - x0) / ((n - 1)**3 * (n - s*s))
        return -numpy.log(x) / t


def alpha_brutal(T0, Ta, t):
        return numpy.log(T0/Ta) / t
