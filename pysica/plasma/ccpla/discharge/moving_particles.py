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

""" Classes for the simulation of charged particles in a plasma discharge.

    This module contains a class to manage charged particles (electrons and ions) moving inside a plasma discharge.
    It provides methods to extract some characteristics of their ensables, such as mean velocity or kinetic energy,
    etc...

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import math
import numpy
from pysica.parameters import *
from pysica.constants import *
from pysica.plasma.ccpla.ccpla_defaults import *
from pysica.managers.io.io_files import write_to_disk


# +------------------+
# | Moving particles |
# +------------------+

class MovingParticles:
        """ A class defining a collection of particles moving in a cold plasma"""
        
        def __init__(self, types, Nmax_particles, start_weight=START_WEIGHT, rescale_factor=DEFAULT_RESCALE_FACTOR):
                """Initialises the properties of a collection of particles moving in a cold plasma .

                   The number of particle types must be given (e.g. 3 if electrons and 2 types of ions are considered)
                   The masses and charges of all particle types are set to the electrons values.
                   All the particles are set as inactive (non existing).
                   The starting position is set to zero.

                        Parameters
                        ----------

                        types:                  number of particle types
                        Nmax_particles:         maximum number of particles which will be possible to manage (for each type)
                                                (dimension of the array in which particles data will be stored)
                        start_weight:           initial value of the particles computational weight 
                        fdt:                    requested value of f*dt in variable dt mode
                        rescale_factor:         the factor by which the number of particles will be divided and the weight 
                                                       will be multiplied if they exceed the maximum allowed number
                                                       by default it is the same for all particle types, 
                                                       but you can set a different value
                                                       for each type by manually changing the values of the 
                                                       self.rescale_factor array



                        Initialized data attributes
                        ---------------------------

                        self.types:             number of particles types
                        self.names:             array containg particle names (may be used for some output)
                        self.n:                 dimension of the array in which particles data are stored
                        self.charge:            particle charge / C (initialized to the electron charge)
                        self.mass:              particle mass / kg  (initialized to the electron mass)
                        self.cm_ratio:          particle charge / mass ratio / C*kg**-1
                        self.number_density:    number density of this particle type / m**-3
                        self.weight:            weight of the particle for simulation purposes,
                                                        i.e. the number of real particles represented by the simulation particle
                        self.rescale_factor:    factor by which the number of particles (of a given type) 
                                                        will be divided and their weight will be multiplied
                                                        if their number exceeds the maximum allowed
                        self.active:            an array identifying which particles are existing at a given time
                                                        True  = active (existing) particle
                                                        False = inactive (non existing) particle
                        self.restart_lf         an array identifying which particles must restart leap-frog scheme
                                                        at next iteration
                                                        True  = restart at next iteration
                                                        False = do not restart

                        self.z                          position along direction parallel to electric field / m
                        self.r                          position along direction normal to electric field / m
                        self.v:                         velocity module / m*s**-1
                        self.vz:                        z components of velocities / m*s**-1
                        self.vx:                        x components of velocities / m*s**-1
                        self.vy:                        y components of velocities / m*s**-1
                        self.vr:                        r components of velocities (x-y plane) / m*s**-1
                        self.theta:                     angle between velocity direction and z axis / rad
                        self.dt:                        timestep for particles motion
#                        self.fdt:                       maximum value of electron collision freq * dt
#                                                            used to calculate dt in variable dt mode
                        self.time:                      actual value of simulation time
                        self.tau_mis:                   approximated value of time between electron collisions
                        self.p_coll:                    approximated value of elctron collision probability

                """
               
                self.types              = int(types)
                self.n                  = int(Nmax_particles)
                self.names              = numpy.zeros(self.types, 'U8')
                self.charge             = numpy.ones(self.types, 'd') * ELECTRON_CHARGE
                self.mass               = numpy.ones(self.types, 'd') * ELECTRON_MASS
                self.cm_ratio           = self.charge / self.mass
                self.number_density     = numpy.zeros(self.types, 'd')
                self.debye_length       = numpy.zeros(self.types, 'd')
                self.plasma_frequency   = numpy.zeros(self.types, 'd')
                self.weight             = numpy.ones(self.types, 'd') * start_weight
                self.rescale_factor     = numpy.ones(self.types, 'd') * rescale_factor
                self.active             = numpy.zeros( (self.types, Nmax_particles), 'bool', order='F')
                self.restart_lf         = numpy.ones ( (self.types, Nmax_particles), 'bool', order='F')
                self.z                  = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.x                  = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.y                  = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.v                  = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.vz                 = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')#self.v*numpy.cos(self.theta)
                self.vx                 = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.vy                 = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.vr                 = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.theta              = numpy.zeros( (self.types, Nmax_particles), 'd',    order='F')
                self.dt                 = 0.0
                self.time               = 0.0
                self.tau_mis            = 0.0
                self.p_coll             = 0.0

                for i in range(self.types):
                        self.names[i] = 'Q_'+str(i)


        def reset(self, start_weight=START_WEIGHT, rescale_factor=DEFAULT_RESCALE_FACTOR):
                """Reset ensamble properties (only mutable ones) to initial values"""
                
                self.number_density.fill(0)
                self.debye_length.fill(0)
                self.plasma_frequency.fill(0)
                self.weight.fill(start_weight)
                self.rescale_factor.fill(rescale_factor)
                self.active.fill(False)
                self.restart_lf.fill(True)
                self.z.fill(0)
                self.x.fill(0)
                self.y.fill(0)
                self.v.fill(0)
                self.vz.fill(0)
                self.vx.fill(0)
                self.vy.fill(0)
                self.vr.fill(0)
                self.theta.fill(0)
                self.dt                 = 0.0
                self.time               = 0.0
                self.tau_mis            = 0.0
                self.p_coll             = 0.0

                        
        def n_active(self, i_type):
                """Returns the number of active particles of the type i_type.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------
                        
                        number of active particles of type i_type
                """

                return self.active[i_type].sum()


        def v_max(self, i_type):
                """Calculates the maximum of active particles speeds.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        self.v:         array of velocity modules / m*s**-1
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------
                        
                        maximum speed of active particles of type i_type / m*s**-1
                """

                if self.active[i_type].any():
                        return self.v[i_type][self.active[i_type]][numpy.argmax( self.v[i_type][self.active[i_type]] )]
                else:
                        return 0.0
                


        def v_min(self, i_type):
                """Calculates the minimum of active particles speeds.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        self.v:         array of velocity modules / m*s**-1
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------

                        return: minimum speed of active particles of type i_type  / m*s**-1
                """

                if self.active[i_type].any():
                        return self.v[i_type][self.active[i_type]][numpy.argmin( self.v[i_type][self.active[i_type]] )]
                else:
                        return 0.0


        def v_average(self, i_type):
                """Calculates the average of active particles speeds.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        
                        self.v:         array of velocity modules / m*s**-1
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------

                        average speed of active particles of type i_type  / m*s**-1

                """
                if self.active[i_type].any():
                        return numpy.average( self.v[i_type][self.active[i_type]] )
                else:
                        return 0.0


        def v2_average(self, i_type):
                """Calculates the average of the square of active particles velocities.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        
                        self.v:         array of velocity modules / m*s**-1
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------

                        average of squared speeds of active particles of type i_type  / (m*s**-1)**2
                """
                
                if self.active[i_type].any():
                        return numpy.average( ( self.v[i_type][self.active[i_type]] )**2 )
                else:
                        return 0.0


        def v_sigma(self, i_type):
                """Calculates the standard deviation of active particles speeds.

                        Parameters
                        ----------      

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        
                        self.v:         array of velocity modules / m*s**-1
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle
                        Returns
                        -------
                        
                        standard deviation of velocity modules of active particles of type i_type  / m*s**-1
                """

                if (self.n_active(i_type) > 1):
                        return numpy.sqrt( sum( ( self.v[i_type][self.active[i_type]] - self.v_average(i_type) )**2 ) / \
                                           (self.n_active(i_type)-1) )
                else:
                        return 0.0


        def all_energies(self, i_type):
                """Returns an array with particles kinetic energies, including inactive ones.

                        Parameters
                        ----------

                        i_type: index that refers to to the type of particle

                        self:   record of particles data

                        self.v: array of velocity modules / m*s**-1

                        Returns
                        -------

                        array of kinetic energies of (active and inactive) particles of type i_type  / eV
                """

                if (self.charge[i_type] != 0):
                        return self.mass[i_type] / (2*numpy.abs(self.charge[i_type])) * self.v[i_type]**2                         
                else:
                        return numpy.zeros(self.n)


                
        def energies(self, i_type):
                """Returns an array with active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:      index that refers to to the type of particle

                        self:        record of particles data

                        self.v:      array of velocity modules / m*s**-1

                        self.active: an array identifying which particles are existing
                                     True  = active particle
                                     False = inactive (non existing) particle

                        Returns
                        -------

                        array of kinetic energies of active particles of type i_type  / eV
                """

                if (self.charge[i_type] != 0):
                        return self.mass[i_type] / (2*numpy.abs(self.charge[i_type])) * self.v[i_type][self.active[i_type]]**2 
                else:
                        return numpy.zeros(self.n_active(i_type))



        def e_average(self, i_type):
                """Calculates the average of active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:  index that refers to to the type of particle

                        self:    record of particles data

                        Returns
                        -------

                        average kinetic energy of active particles of type i_type  / eV
                """

                if ( (self.charge[i_type] != 0) and self.active[i_type].any() ):
                        return self.mass[i_type] / (2*numpy.abs(self.charge[i_type])) * self.v2_average(i_type)                         
                else:
                        return 0.0                        

        def kT(self, i_type):
                """Calculates the peak value (for Boltzmann distribution) of active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:  index that refers to to the type of particle

                        self:    record of particles data

                        Returns
                        -------

                        k*T / eV
                """

                if (self.charge[i_type] != 0):
                        return 2.0 / 3.0 * self.e_average(i_type) # <E> = 3/2 kT                        
                else:
                        return 0.0
                

        def e_min(self, i_type):
                """Calculates the minimum of active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:  index that refers to to the type of particle

                        self:    record of particles data

                        Returns
                        -------

                        minimum kinetic energy of active particles of type i_type  / eV
                """
                
                if ( (self.charge[i_type] != 0) and self.active[i_type].any() ):
                        return self.mass[i_type] / (2*numpy.abs(self.charge[i_type])) * self.v_min(i_type)**2                        
                else:
                        return 0.0



        def e_max(self, i_type):
                """Calculates the maximum of active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:  index that refers to to the type of particle

                        self:    record of particles data

                        Returns
                        -------

                        maximum kinetic energy of active particles of type i_type  / eV
                """
                
                if ( (self.charge[i_type] != 0) and self.active[i_type].any() ):
                        return self.mass[i_type] / (2*numpy.abs(self.charge[i_type])) * self.v_max(i_type)**2                        
                else:
                        return 0.0


        def e_sigma(self, i_type):
                """Calculates the standard deviation of active particles kinetic energies.

                        Parameters
                        ----------

                        i_type:  index that refers to to the type of particle

                        self:    record of particles data

                        Returns
                        -------

                        standard deviation of kinetic energies of active particles of type i_type  / eV
                """

                if (self.n_active(i_type) > 1): # In order to have a finite standard deviation, need at least 2 particles
                        return numpy.sqrt( sum( (self.energies(i_type) - self.e_average(i_type) )**2 ) / (self.n_active(i_type)-1) )
                else:
                        return 0.0



        def theta_average(self, i_type):
                """Calculates the average of active particles angles.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle
        
                        self:           record of particles data
                        
                        self.theta:     angle between velocity direction and z-axis / rad
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------

                        average angle of active particles of type i_type  / deg
                """
                
                if self.active[i_type].any():
                        return numpy.average(self.theta[i_type][self.active[i_type]])
                else:
                        return 0.0


        def theta_sigma(self, i_type):
                """Calculates the standard deviation of active particles angles.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle
        
                        self:           record of particles data
                        
                        self.theta:     angle between velocity direction and z-axis / rad
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------

                        standard deviation of angles of active particles of type i_type  / rad
                """

                if (self.n_active(i_type) > 1):
                        return numpy.sqrt( sum( (self.theta[i_type][self.active[i_type]] - self.theta_average(i_type) )**2 ) / \
                                (self.n_active(i_type)-1.0) )                       
                else:
                        return 0.0
                

        def theta_max(self, i_type):
                """Calculates the maximum of active particles angles.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        self.theta:     array of angles
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------
                        
                        maximum angle / rad
                """
                if self.active[i_type].any():
                        return self.theta[i_type][self.active[i_type]][numpy.argmax(self.theta[i_type][self.active[i_type]])]
                else:
                        return 0.0


        def theta_min(self, i_type):
                """Calculates the minimum of active particles angles.

                        Parameters
                        ----------

                        i_type:         index that refers to to the type of particle

                        self:           record of particles data
                        self.theta:     array of angles
                        self.active:    an array identifying which particles are existing
                                                True  = active particle
                                                False = inactive (non existing) particle

                        Returns
                        -------
                        
                        minimum angle of active particles of type i_type  / rad
                """

                if self.active[i_type].any():
                        return self.theta[i_type][self.active[i_type]][numpy.argmin(self.theta[i_type][self.active[i_type]])]
                else:
                        return 0.0


        def initialize_savefiles(self,
                                 filename_mean_ele, filename_mean_ion,
                                 filename_distrib_ele, filename_distrib_ion,
                                 filename_zpos_ele, filename_zpos_ion,
                                 append=False, sep=SEP, ext='.csv'):
                """Initializes the data files to which simulation data about moving particles will be saved

                   Parameters
                   ----------
                   filename_mean_ele:     name of the file to which the mean electron quantities will be saved:
                                          energy, number density, etc.
                   filename_distrib_ele:  name of the file to which the electron energy distrobution function
                                          will be saved
                   filename_zpos_ele:     name of the file to which the electron z coordinates will be saved
                   filename_mean_ion:     root used to form the names of the files to which the mean ion quantities will be saved:
                                          energy, number density, etc.
                   filename_distrib_ion:  root used to form the names of the files to which the ion energy 
                                          distribution functions will be saved.
                   filename_zpos_ion:     root used to form the names of the files to which the z coordinates 
                                          of the ions will be saved.
                   append:                if True, open the file for appendig new data to the end;
                                          the header will not be re-written
                   sep:                   character or sequence to separate data columns
                   ext:                   extension to give to filenames (deafult is '.csv')


                   Initialized data attributes
                   ---------------------------

                   self.sep:                   character or sequence to separate data columns
                   self.f_mean_ele_name:       name of the file to which values of mean electron quantities will be saved
                   self.f_mean_ion_names:      list of the names of the files to which values of mean ion quantities will be saved
                   self.f_distrib_ele_name:    name of the file to which eedf will be saved
                   self.f_distrib_ion_names:   list of the names of the files to which iedf of each ion type will be saved
                """
                self.sep         = sep
                self.f_mean_ele_name = str(filename_mean_ele) + ext
                
                # Open data file where to save time evolution of mean electron quantities
                if  append:
                        data_file = open(self.f_mean_ele_name,'a')
                else:
                        data_file = open(self.f_mean_ele_name,'w')                        
                        # Write header
                        data_file.write('\"t[ns]\"'         + self.sep )
                        data_file.write('\"N electrons\"'   + self.sep ) # number of active electrons
                        data_file.write('\"w\"'             + self.sep ) # weight
                        data_file.write('\"En[eV]\"'        + self.sep )
                        data_file.write('\"sEn[eV]\"'       + self.sep )
                        data_file.write('\"MinEn[eV]\"'     + self.sep )
                        data_file.write('\"MaxEn[eV]\"'     + self.sep )
                        data_file.write('\"theta[deg]\"'    + self.sep )
                        data_file.write('\"stheta[deg]\"'   + self.sep )
                        data_file.write('\"Mintheta[deg]\"' + self.sep )
                        data_file.write('\"Maxtheta[deg]\"' + self.sep )
                        data_file.write('\"dt[fs]\"'        + self.sep )
                        data_file.write('\"tau[fs]\"'       + self.sep )
                        data_file.write('\"coll_prob[%]\"'  + self.sep )
                        data_file.write('\"n_e[m**-3]\"'    + self.sep )
                        data_file.write(EOL)                        
                data_file.close()

                # Open data files where to save time evolution of mean ion quantities for each ion type
                self.f_mean_ion_names = []
                for i in range(1, self.types):
                        filename = str(filename_mean_ion) + '_' + self.names[i] + ext
                        self.f_mean_ion_names.append(filename)
                        if  append:
                                data_file = open(filename,'a')
                        else:
                                data_file = open(filename,'w')                        
                                # Write header
                                data_file.write( '\"t[ns]\"'         + self.sep )
                                data_file.write( '\"N \"'            + self.sep ) # number of active ions
                                data_file.write( '\"w\"'             + self.sep ) # weight
                                data_file.write( '\"En[eV]\"'        + self.sep )
                                data_file.write( '\"sEn[eV]\"'       + self.sep )
                                data_file.write( '\"MinEn[eV]\"'     + self.sep )
                                data_file.write( '\"MaxEn[eV]\"'     + self.sep )
                                data_file.write( '\"theta[deg]\"'    + self.sep )
                                data_file.write( '\"stheta[deg]\"'   + self.sep )
                                data_file.write( '\"Mintheta[deg]\"' + self.sep )
                                data_file.write( '\"Maxtheta[deg]\"' + self.sep )
                                data_file.write( '\"n [m**-3]\"'     + self.sep )
                                data_file.write(EOL)                        
                        data_file.close()                

                # Open data file where to save evolution of eedf
                self.f_distrib_ele_name = str(filename_distrib_ele) + ext
                # Open data file where to save evolution of eedf
                if (append):
                        data_file = open(self.f_distrib_ele_name,'a')
                else:
                        data_file = open(self.f_distrib_ele_name,'w')
                data_file.close()                        
                        
                # Open data files where to save evolution of iedf for each ion type
                self.f_distrib_ion_names = []
                for i in range(1, self.types):
                        filename = str(filename_distrib_ion) + '_' + self.names[i] + ext
                        self.f_distrib_ion_names.append(filename)                        
                        if append:                                
                                data_file = open(filename,'a')
                        else:
                                data_file = open(filename,'w')
                        data_file.close()

                # Open data file where to save evolution of electron z coordinates
                self.f_zpos_ele_name = str(filename_zpos_ele) + ext
                if (append):
                        data_file = open(self.f_zpos_ele_name,'a')
                else:
                        data_file = open(self.f_zpos_ele_name,'w')
                data_file.close()

                # Open data files where to save evolution of ion z coordinates
                self.f_zpos_ion_names = []
                for i in range(1, self.types):
                        filename = str(filename_zpos_ion) + '_' + self.names[i] + ext
                        self.f_zpos_ion_names.append(filename)                        
                        if append:                                
                                data_file = open(filename,'a')
                        else:
                                data_file = open(filename,'w')
                        data_file.close()                
                
                        
        def save_data_to_files(self, save_edf=True, save_z=True):
                """Saves actual data values to files """

                # Save data on electron mean quantities
                data_file = open(self.f_mean_ele_name, 'a')                             
                data_file.write( str( 1E9*self.time                      ) + self.sep )
                data_file.write( str( self.n_active(0)                   ) + self.sep ) # number of active electrons
                data_file.write( str( self.weight[0]                     ) + self.sep ) # weight
                data_file.write( str( self.e_average(0)                  ) + self.sep )
                data_file.write( str( self.e_sigma(0)                    ) + self.sep )
                data_file.write( str( self.e_min(0)                      ) + self.sep )
                data_file.write( str( self.e_max(0)                      ) + self.sep )
                data_file.write( str( math.degrees(self.theta_average(0))) + self.sep )       
                data_file.write( str( math.degrees(self.theta_sigma(0))  ) + self.sep )
                data_file.write( str( math.degrees(self.theta_min(0))    ) + self.sep )
                data_file.write( str( math.degrees(self.theta_max(0))    ) + self.sep )
                data_file.write( str( 1.0E15*self.dt                     ) + self.sep )
                data_file.write( str( 1.0E15*self.tau_mis                ) + self.sep )
                data_file.write( str( self.p_coll                        ) + self.sep )
                data_file.write( str( self.number_density[0]             ) + self.sep )                 
                data_file.write(EOL)
                data_file.close()
                
                # Save data on ion mean quantities
                for i in range(1, self.types):
                        i_file = i-1
                        data_file = open(self.f_mean_ion_names[i_file], 'a')                             
                        data_file.write( str( 1E9*self.time                      ) + self.sep )
                        data_file.write( str( self.n_active(i)                   ) + self.sep ) # number of active ions
                        data_file.write( str( self.weight[i]                     ) + self.sep ) # weight
                        data_file.write( str( self.e_average(i)                  ) + self.sep )
                        data_file.write( str( self.e_sigma(i)                    ) + self.sep )
                        data_file.write( str( self.e_min(i)                      ) + self.sep )
                        data_file.write( str( self.e_max(i)                      ) + self.sep )
                        data_file.write( str( math.degrees(self.theta_average(i))) + self.sep )       
                        data_file.write( str( math.degrees(self.theta_sigma(i))  ) + self.sep )
                        data_file.write( str( math.degrees(self.theta_min(i))    ) + self.sep )
                        data_file.write( str( math.degrees(self.theta_max(i))    ) + self.sep )
                        data_file.write( str( self.number_density[i]             ) + self.sep )                 
                        data_file.write(EOL)
                        data_file.close()

                if save_edf:                
                        # Save actual eedf
                        data_file = open(self.f_distrib_ele_name,'a')
                        self.e_distr = self.energies(0)
                        for j in range(self.n_active(0)):
                                data_file.write( str(self.e_distr[j]) + self.sep )
                        data_file.write( EOL )
                        data_file.close()
                        # Save actual iedf for all ion types
                        for i in range(1,self.types):
                                i_file = i-1
                                data_file = open(self.f_distrib_ion_names[i_file], 'a')
                                self.e_distr = self.energies(i)
                                for j in range(self.n_active(i)):
                                        data_file.write( str(self.e_distr[j]) + self.sep )
                                data_file.write( EOL )
                                data_file.close()

                if save_z:                
                        # Save actual elecron z coordinates
                        data_file = open(self.f_zpos_ele_name,'a')
                        for j in range(self.n_active(0)):
                                data_file.write( str(self.z[0][self.active[0]][j]) + self.sep )
                        data_file.write( EOL )
                        data_file.close()                        
                        # Save actual ion z coordinates for all ion types
                        for i in range(1, self.types):
                                i_file = i-1
                                data_file = open(self.f_zpos_ion_names[i_file], 'a')
                                #self.e_distr = self.energies(i)
                                for j in range(self.n_active(i)):
                                        data_file.write( str(self.z[i][self.active[i]][j]) + self.sep )
                                data_file.write( EOL )
                                data_file.close()                
