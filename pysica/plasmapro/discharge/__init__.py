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
	Tools for the simulation of plasma discharges: definition of particles ensambles and simulation of particle motion

        The modules in this subpackage define classes to store data about electrons, ions and neutrals.
        They also provide functions which call fortran functions compiled into the fortran package 
        to simulate the motion and interaction of these particles inside plasma discharges,
 
       	* modules
		- reactors 	         -> defines characteristics of plasma reactors, such as electrodes dimensions, bias, etc...
                - moving_particles       -> defines the ensables of electrons and ions
		- target_particles       -> defines neutral particles and their cross sections for impact with electrons and ions

                - particles_data_manager -> read properties of neutral particles from a ascii file
                - particle_mover         -> provides functions to simulate the particles motion and interaction in a cold plasma
                                            it uses the modules contained in the fortran package
"""
