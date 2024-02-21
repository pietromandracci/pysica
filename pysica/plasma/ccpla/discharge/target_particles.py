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
        Classes for the simulation of neutral and charged particles in a plasma discharge.

        This module contains a class for managing the neutral particles (mainly atoms and molecules) 
        which are the targets of collisions by moving particles, usually electrons and ions.
        These target particles are considered fixed in space, since their speed is much lower respect to moving particles.
        Cross sections needed to treat the collision processes are imported from ASCII files.

        How collisions are treated
        --------------------------

        To decide if an electron (or ion) scatters whith a neutral (by any type of scattering process: elastic, 
        ionization, dissociation of any type) during a time interval dt, we consider that the 
        probability for (at least) a scattering event to happen during a time interval dt is
        p = 1 - exp( f(v) * dt ) = 1 - exp( r(v) * v * dt ) = 1 - exp( s(v) * n * v * dt )
        where: 
                v           -> electron speed
                n           -> total neutral density (sum of densities of all neutral types)
                f(v)=r(v)*v -> scattering frequency (number of scattering events experienced by a single electron in a unit time)
                r(v)=s(v)*n -> total scattering rate (number of scattering events for a unit length of 
                               path travelled by the electron)
                s(v)        -> total cross section (sum of all cross sections for all scattering types) 
                               at speed v

        Once decided that a scattering has happened, to decide of which type it was, we have to use
        relative probabilities and probability limits.
        The relative probability p_i that the scattering was of type i is
                p_i = r_i(v)/r(v)
        where:
                i      -> identifies the type of scattering process
                          0                <= i < self.types-1      -> elastic scattering (with neutral type i)
                          self.types       <= i < 2*self.types-1    -> ionization 
                                                                       (of neutral type i-self.types)
                          2*self.types     <= i < 2*self.types   
                                                  + n_excitations   -> excitation
                                                                       (neutral type and dissociation
                                                                       process depending on the value of i)
                          2*self.types-1 +     
                          + n_excitations  <= i                     -> dissociation 
                                                                       (neutral type and dissociation
                                                                       process depending on the value of i)

                r_i(v) -> scattering rate for process i at speed v
                r(v)   -> total scattering rate for all processes at speed v
                of course, the sum of all p_i gives 1.
                To decide the type of scattering, we extract a random numbed rand between 0 and 1, than compare the
                        number with the limits PL_i defined as follows
                        PL_0 = p_0
                        PL_i = PL_(i-1) + p_i           for i>0
                The type of scattering is given by the value of i for which it is true that
                        PL_(i-1) < rand < PL_i          whith i>0
                        or
                        0        < rand < PL_0          if i=0

        Recombination
        -------------

        The same considerations hold for recombination of electrons with ions, by the three main recombination channels:
        dissociative recombination, three body recombination and radiative recombination.
        In this case, we must substuitute the neutral density n with the ion density n_ion. However, while the former
        can be considered approximately constant if the ionization degree is very low, the latter is changing in time
        as the ionization process proceeds.

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
from pysica.functions.physics import number_density, mean_speed_maxwell
from pysica.managers.data_manager import DataGrid
from pysica.managers.gnuplot_manager import gnuplot_installed, new_plot, plot2d, plot_curves, plot_close, plot_close_all
from pysica.managers.io.io_files import write_to_disk
from pysica.plasma.ccpla.discharge.particles_data_manager import read_file_ntypes, read_file_properties


XSEC_PLOT_PERSIST = True

# +------------------+
# | Target particles |
# +------------------+

class TargetParticles:
        """A class defining a collection of particles, to be considered as the targets of collisions by electrons and ions."""

        def __init__(self, N_sigma_electrons, N_sigma_ions, temperature, total_pressure, min_scattered=DEFAULT_MIN_SCATTERED,
                     isactive_recomb=True, neutral_types=None, filename=None):
                """ Initialises the properties of target particles.

                        Initialises the properties of target particles.
                        If the number of target particle types is not specified, it will be read from a ASCII file containing
                        particles properties.

                        Parameters
                        ----------

                        N_sigma_electrons: number of energies for which cross section for electrons scattering values are tabulated
                        N_sigma_ions:      number of energies for which cross section for ions scatterig values are tabulated
                        temperature:       temperature of gases in K
                        total_pressure:    total pressure of gases in Pa
                        isactive_recomb:   set if electron/ion recombination processes must be activated or not
                        neutral_types:     number of neutral types, if this value is equal to None
                                              than the number of neutrals will be read from 
                                              the ASCII file with name filename
                        filename:          name of the ASCII file containing neutrals properties
                                              if neutral_types=None this parameter cannot be left blank

                        Initialized data attributes
                        ---------------------------

                        self.read_error:                used for error checking
                                                        (status, message)
                                                                status: 0 = no error
                                                                        1 = file not found
                                                                        2 = no gases were given in the file
                                                                message: an error message or 'OK'
                        self.temperature:               temperature of the gas mixture / K
                        self.total_pressure:            total pressure of the gas mixture / Pa
                        self.isactive_recomb:           set if electron/ion recombination processes must be activated or not
                        self.types:                     how many types of target atoms/molecules are involved
                        self.n_sigma:                   number of velocity values for which electron impact cross sections 
                                                               are tabulated
                        self.n_sigma_ions:              number of velocity values for which ion impact cross sections are tabulated
                """

                self.read_error = (0, 'OK')

                # If the number of gases was not given, try to read it from file
                if (neutral_types is None):
                        if (filename is None):
                                self.read_error = (3, "Filename was not given")
                                return
                        else:
                                self.read_error = read_file_ntypes(filename, self)
                                if (self.read_error[0] != 0): return
                else:
                        self.types              = neutral_types

                # Characteristics of the plasma discharge
                self.temperature                = temperature
                self.total_pressure             = total_pressure
                self.isactive_recomb            = isactive_recomb

                # Number of velocity values used to discretize cross sections for electron and ion impact with neutrals
                self.n_sigma                    = N_sigma_electrons
                self.n_sigma_ions               = N_sigma_ions
                self.min_scattered              = min_scattered


        # +-------------------------------------------------------+
        # | Get energies for which cross section values are known |
        # +-------------------------------------------------------+

        def sigma_energies(self):
                """Returns an array with the energies (in eV) for which cross section values for electrons are tabulated.

                        Returns
                        -------

                        array of the energies (expressed in eV) for which the cross section values are tabulated

                """

                return - 0.5 * ELECTRON_MASS / ELECTRON_CHARGE * self.sigma_velocities**2


        def sigma_energies_ions(self, ion_type):
                """Returns an array with the energies (in eV) for which cross section values for ions are tabulated.

                        Parameters
                        ----------

                        ion_type: the type of ion for the energies are requested


                        Returns
                        -------

                        array of the energies (expressed in eV) for which the cross section values are tabulated

                """

                return - 0.5 * self.mass[ion_type] * ATOMIC_UNIT_MASS / ELECTRON_CHARGE * self.sigma_velocities_ions[ion_type]**2


        # +--------------------------------+
        # | Get the properties of neutrals |
        # +--------------------------------+

        def read_properties(self, filename, sep, e_min_sigma, e_max_sigma, e_min_sigma_ions, e_max_sigma_ions, debug=False):
                """Reads the neutrals characteristics from an ascii file and calculates some properties.

                        Parameters
                        ----------

                        filename:         name of file containing neutral gases characteristics
                        sep:              character used in the ASCII file to separate values (usually the tab char)
                        e_min_sigma:      minimum energy for which cross sections for electron scattering are tabulated
                        e_max_sigma:      maximum energy for which cross sections for electron scattering are tabulated
                        e_min_sigma_ions: minimum energy for which cross sections for ion scattering are tabulated
                        e_max_sigma_ions: maximum energy for which cross sections for ion scattering are tabulated

                        Initialized data attributes
                        ---------------------------
                        self.names:                      names of the gas molecules (max 8 char)
                        self.types_atoms:                number of atom types
                        self.types_molecules:            number of molecule types (not atoms)
                        self.molecule_type:              string describing the type of molecule (1 char)
                                                                'a'= atom, 'm'= nonpolar molecule, 'p'= polar molecule
                                                                actually polar and non-polar molecules are treated in the 
                                                                same way, but this may change in the future
                        self.mass:                       masses of gas molecules / u
                        self.mean_v:                     mean speeds of gas molecules (all at the same temperature)
                        self.gas_flow:                   array of gas fluxes / arbitrary units
                        self.partial_pressure:           array of gases partial pressures / Pa
                        self.number_density:             number density of neutral species
                        self.secondary_emission:         secondary emission coeff. i.e. average number of e- emitted 
                                                                when an ion is collected at one electrode
                        self.dissociation_rate           dissociation rate for each molecule type and its dissociation types
                        self.dissociation_rate_constant  dissociation rate constant for each molecule type and its dissociation types

                        self.energy_loss:                fraction of energy transferred from an electron to a neutral particle
                                                                during elastic scattering
                        self.energy_loss_ions:           fraction of energy transferred from an ion to a neutral particle
                                                                during elastic scattering
                        self.v_ratio:                    electron speed reduction during elastic scattering
                        self.v_ratio_ions:               ion speed reduction during elastic scattering
                        self.ionization_energy:          first ionization energy / eV
                        self.excitation_tpyes            number of excitation processes for each neutral type
                        self.excitation_energy:          excitation energies for processes to be considered / eV
                        self.dissociation_types:         number of dissociation processes for each neutral type
                        self.dissociation_energy:        dissociation energies for processes to be considered / eV

                        self.sigma_velocities:           array of velocity values for which electron impact cross sections 
                                                                are tabulated
                        self.sigma_velocities_ions:      matrix of velocity values for which cross sections are tabulated 
                                                                for each ion type

                        self.sigma_elastic:              cross section table for electron impact elastic scattering 
                        self.sigma_ionization:           cross section table for electron impact ionization
                        self.sigma_excitation            cross section table for electron impact excitation
                        self.sigma_dissociation:         cross section table for electron impact dissociation   

                        self.sigma_recombination_diss:   cross section table for electron impact dissociative recombination 
                                                                of ionized molecules
#                       self.sigma_recombination_3body:  cross section table for electron impact three body recombination 
#                                                                of ionized atoms/molecules
#                       self.sigma_recombination_rad:    cross section table for electron impact radiative recombination 
#                                                                of ionized atoms/molecules

                        self.sigma_elastic_ions:         cross section table for ion impact elastic scattering 
                        self.sigma_charge_exchange_ions: cross section table for ion impact elastic scattering 

                        self.collisions_null:            number of electron null collisions
                        self.collisions_elastic:         number of electron elastic collisions
                        self.collisions_ionization:      number of electron ionization collisions
                        self.collisions_excitation:      number of electron excitation collisions
                        self.collisions_dissociation:    number of electron dissociation collisions
                        self.collisions_recombination:   number of electron/ion recombination collisions
                        self.collisions_total_electron:  total number of electron collisions, of any type
                        self.dissociation_rate:          dissociation rate of neutral particles, 
                                                               for each type and dissociation process
                        self.dissociation_rate_constant: dissociation rate constant of neutral particles, 
                                                               for each type and dissociation process

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = could not open file
                                         2 = missing values
                                         3 = numerical value expected
                                         4 = unknown molecule type
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                # Characteristics of gases
                self.names                      = numpy.zeros(self.types, 'U8')
                self.molecule_type              = numpy.zeros(self.types, 'U1')
                self.mass                       = numpy.zeros(self.types, 'd')
                self.mean_v                     = numpy.zeros(self.types, 'd')
                self.gas_flow                   = numpy.zeros(self.types, 'd')
                self.partial_pressure           = numpy.zeros(self.types, 'd')
                self.number_density             = numpy.zeros(self.types, 'd')
                self.secondary_emission         = numpy.zeros(self.types, 'd')

                # Characteristics for interaction with electrons
                self.energy_loss                = numpy.zeros(self.types, 'd')
#                self.v_ratio                    = numpy.zeros(self.types, 'd')
                self.ionization_energy          = numpy.zeros(self.types, 'd')
                self.excitation_types           = numpy.zeros(self.types, 'b')
                self.excitation_energy          = numpy.zeros( (self.types, MAX_EXCITATION_TYPES), 'd', order='F' )
                self.dissociation_types         = numpy.zeros(self.types, 'b')
                self.dissociation_energy        = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES), 'd', order='F' )

                # Characteristics for interaction with ions
                self.energy_loss_ions           = numpy.zeros( (self.types, self.types), 'd', order='F' )
#                self.v_ratio_ions               = numpy.zeros( (self.types, self.types), 'd', order='F' )

                # Arrays of velocity values for which cross sections are tabulated
                self.sigma_velocities           = numpy.zeros(self.n_sigma, 'd')        
                self.sigma_velocities_ions      = numpy.zeros( (self.types, self.n_sigma_ions), 'd', order='F') 

                # Arrays of cross sections for electron impact processes
                self.sigma_elastic              = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')
                self.sigma_ionization           = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')
                self.sigma_excitation           = numpy.zeros( (self.types, MAX_EXCITATION_TYPES, self.n_sigma), 
                                                             'd', order = 'F')
                self.sigma_dissociation         = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES, self.n_sigma), 
                                                             'd', order = 'F')

                self.sigma_recombination_diss   = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES, self.n_sigma), 
                                                             'd', order = 'F')
                self.sigma_recombination_3body  = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')
                self.sigma_recombination_rad    = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')

                # Arrays of cross sections for ion impact processes
                self.sigma_elastic_ions         = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')
                self.sigma_charge_exchange_ions = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')

                # Quantities that will be calculated during the simulation
                # collision counters for electron/neutral collisions and electron/ion recombinations
                # and dissociation rate / rate constant counters
                self.collisions_null            = numpy.zeros(1, 'd') # Must be an array since it will be passed to f2py
                self.collisions_elastic         = numpy.zeros(self.types, 'd')
                self.collisions_ionization      = numpy.zeros(self.types, 'd')
                self.collisions_excitation      = numpy.zeros( (self.types, MAX_EXCITATION_TYPES),   'd', order='F' )
                self.collisions_dissociation    = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES), 'd', order='F' )
                self.collisions_recombination   = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES), 'd', order='F' )
                self.collisions_total_electron  = 0
                self.dissociation_rate          = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES), 'd', order='F' )
                self.dissociation_rate_constant = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES), 'd', order='F' )               

               
                # Read neutral properties from ascii file
                (status, message) = read_file_properties(filename, sep, self, debug)
                if (status != 0): return (status, message)

                # Calculate the number of types of atoms and molecules
                self.types_atoms = 0
                for i in range(self.types):
                        if (self.molecule_type[i] == 'a'): self.types_atoms += 1
                self.types_molecules = self.types - self.types_atoms

                # Calculate neutral properties
                self.partial_pressure   = self.gas_flow * self.total_pressure / self.gas_flow.sum()
                self.number_density     = number_density(self.partial_pressure, self.temperature)
                
                # Average electron energy lost during an elastic scattering
                self.energy_loss        = 2 * ELECTRON_MASS / (self.mass * ATOMIC_UNIT_MASS)
                # Average ratio of electron speed measured after and before an elastic scattering
#                self.v_ratio            = numpy.sqrt(numpy.abs(1-self.energy_loss))                     
                self.mean_v             = mean_speed_maxwell(self.temperature, self.mass)

                for i in range(self.types):
                        for j in range(self.types):
                                # Average ion energy lost during an elastic scattering
                                self.energy_loss_ions[i][j] = 2 * self.mass[i] * self.mass[j] / (self.mass[i] + self.mass[j])**2
                                # Average ratio of ion speed measured after and before an elastic scattering
#                                self.v_ratio_ions[i][j]     = numpy.sqrt(1-self.energy_loss_ions[i][j])

                # Calculate speed values at which cross sections for electrons will be interpolated
                v_min_sigma = math.sqrt( 2 * e_min_sigma * (-ELECTRON_CHARGE / ELECTRON_MASS) )
                v_max_sigma = math.sqrt( 2 * e_max_sigma * (-ELECTRON_CHARGE / ELECTRON_MASS) )
                self.sigma_velocities = numpy.linspace(v_min_sigma, v_max_sigma, self.n_sigma )

                # Calculate speed values at which cross sections for each ion type will be interpolated
                # kinetic energy values are stored in eV, so they are multiplied by electron charge
                # to transform them in joule (ELECTRON_CHARGE is a negative number)
                for i in range(self.types):
                        v_min_sigma = math.sqrt( 2 * e_min_sigma_ions * (-ELECTRON_CHARGE / (self.mass[i] * ATOMIC_UNIT_MASS)) )
                        v_max_sigma = math.sqrt( 2 * e_max_sigma_ions * (-ELECTRON_CHARGE / (self.mass[i] * ATOMIC_UNIT_MASS)) )
                        self.sigma_velocities_ions[i] = numpy.linspace(v_min_sigma, v_max_sigma, self.n_sigma_ions)

                return (status, message)


        # +--------------------------------------------------------------------+
        # | General method used to read cross section values from a ascii file |
        # +--------------------------------------------------------------------+

        def _read_xsec(self, filename, sep, qm_ratio, target_index, ion_index=None, plot=False, title=''):
                """Reads a table of cross section values as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form
                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment
                                        where energies are in eV and cross section values are in m**2
                                        
                        sep:            character used to separate values (usually the tab char)

                        qm_ratio:       the charge to mass ratio of incident particle

                        ion_index:      integer value indentifying the type of ion (i.e. the type of neutral that was ionized) 
                                        that scatters with the target, 
                                        if it is None, than the particle is an electron

                        target_index:   integer value indentifying the type of target

                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation

                        title:          a string containing the title to be printed on the plot


                        Initialized data attributes
                        ---------------------------

                        self.sigma_return       numpy array containing the cross section values that were read from the file


                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                         3 = ion index out of bounds
                                message: a string containing an error message or 'Ok'
                """


                status, message  = 0, 'OK'              # This variable is used for error checking      

                if (ion_index is None):
                        self.sigma_return= numpy.zeros(self.n_sigma, 'd')
                else:
                        if ( (ion_index < 0) or (ion_index > self.types-1)):
                                status = 3
                                message = 'ion type ' + str(ion_index) + ' is out of bounds'
                                return(status, message)                                
                        self.sigma_return= numpy.zeros(self.n_sigma_ions, 'd')

                # Read the tabulated values of elastic cross-section into a numpy array
                sigma = DataGrid()
                (status, message) = sigma.read_file(filename, n_columns=2, sep=sep, transpose=True)
                if (status != 0):
                        return (status, message)
                else: 
                        sigma.n_data = len(sigma.data_array[0])

                # Get the absolute value of mass to charge ratio
                if (qm_ratio < 0): qm_ratio = - qm_ratio

                # Check that all cross-section values are positive
                for j in range(sigma.n_data):
                        if sigma.data_array[0][j] < 0:  
                                status  = 1
                                message = 'negative energy value found in row # '+str(j+1)
                                return (status, message)
                        if sigma.data_array[1][j] < 0:
                                status  = 2
                                message = 'negative cross-section value found in row # '+str(j+1)
                                return (status, message)

                # Change from energies [eV] to velocities [m/s]
                sigma.data_array[0] = numpy.sqrt( 2 * qm_ratio * sigma.data_array[0] )

                # Calculate a interpolated cross-section table
                if (ion_index is None): 
                        self.sigma_return = numpy.interp(self.sigma_velocities,
                                                         sigma.data_array[0],
                                                         sigma.data_array[1])
                else:
                        self.sigma_return = numpy.interp(self.sigma_velocities_ions[ion_index],
                                                         sigma.data_array[0],
                                                         sigma.data_array[1])

                # Plot cross-section
                if (plot and gnuplot_installed): 
                        xsec_graph = new_plot(title=title,
                                              xlabel='Energy / eV',
                                              ylabel='Cross section / m**2',
                                              format_x='%g',
                                              format_y='%g',
                                              logx=True,
                                              grid=True,
                                              persistence=XSEC_PLOT_PERSIST,
                                              redirect_output=True)
                       
                        data_list = []
                        
                        # These are the values given in the data file
                        data_list.append( [0.5 / qm_ratio * sigma.data_array[0]**2,
                                                sigma.data_array[1],
                                                'Tabulated values',
                                                None ] )
                        
                        # These are the values calculated by interpolation of the ones in the data file
                        if (ion_index is None):
                                data_list.append( [ self.sigma_energies(),
                                                         self.sigma_return,
                                                         'Interpolated values',
                                                         None ] )
                        else:
                                data_list.append( [ self.sigma_energies_ions(ion_index),
                                                         self.sigma_return,
                                                         'Interpolated values',
                                                         None ] )
                       
                        plot_curves(xsec_graph, data_list, volatile=True)
                       
                        # Close the plots, to quit the associated gnuplot processes
                        # plots remain visible thanks to the "persist" option
                        plot_close(xsec_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)

                        del(xsec_graph)

                return (status, message)


        # +------------------------------------------------------------------+
        # | Methods to read cross sections for electrons/neutrals scattering |
        # +------------------------------------------------------------------+

        def read_xsec_electrons_elastic(self, filename, sep, neutral_index, plot=False):
                """Reads the table of elastic cross section values as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                        
                        sep:            character used to separate values (usually the tab char)
                        neutral_index:  integer value indentifying the type of neutral for which cross section
                                        shuld be read

                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation


                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for elastic scattering of e- on '+self.names[neutral_index]
                qm_ratio = - ELECTRON_CHARGE / ELECTRON_MASS

                (status, message) = self._read_xsec(filename, sep, qm_ratio, neutral_index, None, plot, title)

                if (status == 0): self.sigma_elastic[neutral_index] = self.sigma_return

                return (status, message)


        def read_xsec_electrons_ionization(self, filename, sep, neutral_index, check=False, plot=False):
                """Reads the table of ionization cross section values as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form

                                                        #Comment
                                                        <energy> <separator> <value>    #Comment
                                                        <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                                        

                        sep:            character used to separate values (usually the tab char)

                        neutral_index:  integer value indentifying the type of neutral for which cross section
                                        shuld be read

                        check:          if True, check that all cross section values for energies lower
                                        than ionization threshold are zero

                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation


                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for e- impact ionization of '+self.names[neutral_index]
                qm_ratio = - ELECTRON_CHARGE / ELECTRON_MASS

                (status, message) = self._read_xsec(filename, sep, qm_ratio, neutral_index, None, plot, title)

                if (status == 0): self.sigma_ionization[neutral_index] = self.sigma_return

                # Check that cross section is zero for energies lower than ionization threshold
                if check:
                        # Calculate speed corresponding to ionization treshold
                        v_min = math.sqrt( 2.0 * self.ionization_energy[neutral_index] * qm_ratio )
                        i=0
                        while(self.sigma_velocities[i]<=v_min): 
                                self.sigma_ionization[neutral_index, i]=0
                                i=i+1

                return (status, message)


        def read_xsec_electrons_excitation(self, filename, sep, neutral_index, exc_type, check=False, plot=False):
                """Reads the table of excitation cross section values as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of the file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                                

                        sep:            character used to separate values (usually the tab char)

                        neutral_index:  integer value indentifying the type of neutral for which cross section
                                                shuld be read

                        exc_type:       integer value identifying the type of excitation process

                        check:          if True, check that all cross section values for energies lower
                                        than dissociation threshold are zero

                        plot:           if True, plot the cross section values as a function of energy,
                                                before and after interpolation

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for e- impact excitation of '+ self.names[neutral_index]+' process #'+str(exc_type)
                qm_ratio = - ELECTRON_CHARGE / ELECTRON_MASS

                (status, message) = self._read_xsec(filename, sep, qm_ratio, neutral_index, None, plot, title)

                if (status == 0): self.sigma_excitation[neutral_index, exc_type] = self.sigma_return

                # Check that cross section is zero for energies lower than excitation threshold
                if check:
                        # Calculate speed corresponding to excitation treshold
                        v_min = math.sqrt( 2.0 * self.excitation_energy[neutral_index, exc_type] * qm_ratio )
                        i=0
                        while(self.sigma_velocities[i]<=v_min): 
                                self.sigma_excitation[neutral_index, exc_type, i]=0
                                i=i+1                   

                return (status, message)
        
        

        def read_xsec_electrons_dissociation(self, filename, sep, neutral_index, diss_type, check=False, plot=False):
                """Reads the table of dissociation cross section values as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of the file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                                

                        sep:            character used to separate values (usually the tab char)

                        neutral_index:  integer value indentifying the type of neutral for which cross section
                                                        shuld be read

                        diss_type:      integer value identifying the type of dissociation process

                        check:          if True, check that all cross section values for energies lower
                                        than dissociation threshold are zero

                        plot:           if True, plot the cross section values as a function of energy,
                                                before and after interpolation

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for e- impact dissociation of '+ self.names[neutral_index]+' process #'+str(diss_type)
                qm_ratio = - ELECTRON_CHARGE / ELECTRON_MASS

                (status, message) = self._read_xsec(filename, sep, qm_ratio, neutral_index, None, plot, title)

                if (status == 0): self.sigma_dissociation[neutral_index, diss_type] =self.sigma_return

                # Check that cross section is zero for energies lower than dissociation threshold
                if check:
                        # Calculate speed corresponding to dissociation treshold
                        v_min = math.sqrt( 2.0 * self.dissociation_energy[neutral_index, diss_type] * qm_ratio )                      
                        i=0
                        while(self.sigma_velocities[i]<=v_min): 
                                self.sigma_dissociation[neutral_index, diss_type, i]=0
                                i=i+1                   

                return (status, message)


        def calculate_total_xsec_electrons(self, debug=False):
                """Calculates main statistical parameters from cross sections.

                        Calculate, for the different processes (elastic scattering, ionization, dissociation types): 
                                - total cross sections
                                - collision frequencies
                                - relative scattering probabilities
                                - probability boundaries
                        (see the module docstring for more information on these quantities and their usage)

                        Parameters
                        ----------

                        Initialized data attributes
                        ---------------------------

                        self.sigma_total_elastic:       sum of elastic scattering cross sections for all neutral types
                        self.sigma_total_ionization:    sum of ionization cross sections for all neutral types
                        self.sigma_total_excitation:    sum of excitation cross sections for all neutral types
                                                                and all excitation processes for each neutral type
                        self.sigma_total_dissociation:  sum of dissociation cross sections for all neutral types
                                                                and all dissociation processes for each neutral type
                        self.sigma_total_global:        sum of all cross section for all processes

                        self.frequency_elastic:         collision frequencies for elastic scattering
                        self.frequency_ionization:      collision frequencies for ionization scattering
                        self.frequency_excitation:      collision frequencies for excitation scattering
                        self.frequency_dissociation:    collision frequencies for dissociation scattering
                        self.frequency_total_global:    sum of all collision frequencies, for all prossible processes
                        self.frequency_null_collision:  null collision frequency, defined as the difference between
                                                                the maximum value (over all velocities) of the 
                                                                total collision frequency
                                                                and the collision frequency

                        self.p_relative_elastic:        relative probability of elastic scattering
                        self.p_relative_ionization:     relative probability of ionization
                        self.p_relative_excitation:     relative probabilities excitation processes
                        self.p_relative_dissociation:   relative probabilities dissociation processes
                        self.p_relative_null_collision: relative probability of null collision                                  

                        self.nmax_limits:               maximum allowed number of probability limits
                                                        i.e. max number of processes
                        self.probability_limits:        limits for calculation of scattering type
                        self.probability_processes:     names of the processed associated to limits
                        self.probability_index_exc      index of the first probability limit for excitation
                        self.probability_index_diss     index of the first probability limit for dissociation

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                self.sigma_total_elastic        = numpy.zeros( self.n_sigma, 'd')
                self.sigma_total_ionization     = numpy.zeros( self.n_sigma, 'd')
                self.sigma_total_excitation     = numpy.zeros( self.n_sigma, 'd')
                self.sigma_total_dissociation   = numpy.zeros( self.n_sigma, 'd')
                self.sigma_total_global         = numpy.zeros( self.n_sigma, 'd')

                self.frequency_elastic          = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')
                self.frequency_ionization       = numpy.zeros( (self.types, self.n_sigma), 'd', order = 'F')
                self.frequency_excitation       = numpy.zeros( (self.types, MAX_EXCITATION_TYPES, self.n_sigma), 
                                                               'd', order = 'F')                
                self.frequency_dissociation     = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES, self.n_sigma), 
                                                               'd', order = 'F')

                self.frequency_total_global     = numpy.zeros( self.n_sigma, 'd')
                self.frequency_null_collision   = numpy.zeros( self.n_sigma, 'd')

                self.p_relative_elastic         = numpy.zeros( (self.types, self.n_sigma), 'd', order='F')
                self.p_relative_ionization      = numpy.zeros( (self.types, self.n_sigma), 'd', order='F')
                self.p_relative_excitation      = numpy.zeros( (self.types, MAX_EXCITATION_TYPES, self.n_sigma), 
                                                                'd', order = 'F')                
                self.p_relative_dissociation    = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES, self.n_sigma), 
                                                                'd', order = 'F')
                self.p_relative_null_collision  = numpy.zeros( self.n_sigma, 'd')

                self.nmax_limits                = self.types * (2 + MAX_EXCITATION_TYPES + MAX_DISSOCIATION_TYPES)
                self.probability_limits         = numpy.zeros( (self.nmax_limits, self.n_sigma), 'd', order = 'F')
                self.probability_processes      = numpy.zeros( self.nmax_limits, 'U24')

                # Calculate total cross sections for all processes
                for i in range(self.types):
                        self.sigma_total_elastic    = self.sigma_total_elastic    + self.sigma_elastic[i]
                        self.sigma_total_ionization = self.sigma_total_ionization + self.sigma_ionization[i]
                        for j in range(self.excitation_types[i]):
                                self.sigma_total_excitation = ( self.sigma_total_excitation
                                                                + self.sigma_excitation[i,j] )
                        for j in range(self.dissociation_types[i]):
                                self.sigma_total_dissociation = ( self.sigma_total_dissociation
                                                                  + self.sigma_dissociation[i,j] )

                # Calculate total global cross-section
                self.sigma_total_global = (self.sigma_total_elastic      + self.sigma_total_ionization
                                           + self.sigma_total_excitation + self.sigma_total_dissociation  )

                # Calculate collision frequencies for all processes
                # f_i(v) = n_i * s_i(v) * v
                for i in range(self.types):
                        self.frequency_elastic[i]    = self.number_density[i] * self.sigma_elastic[i]    * self.sigma_velocities 
                        self.frequency_ionization[i] = self.number_density[i] * self.sigma_ionization[i] * self.sigma_velocities
                        for j in range(self.excitation_types[i]):
                                self.frequency_excitation[i,j] = ( self.number_density[i]
                                                                   * self.sigma_excitation[i,j]
                                                                   * self.sigma_velocities )
                        for j in range(self.dissociation_types[i]):
                                self.frequency_dissociation[i,j] = ( self.number_density[i]
                                                                     * self.sigma_dissociation[i,j]
                                                                     * self.sigma_velocities )

                # Calculate total collision frequency 
                # f_tot(v) = sum{f_i(v), i};  f_i(v) = n_i * s_i(v) * v
                for i in range(self.types):
                        self.frequency_total_global = ( self.frequency_total_global
                                                        + self.frequency_elastic[i]
                                                        + self.frequency_ionization[i] )
                        for j in range(self.excitation_types[i]):
                                self.frequency_total_global = ( self.frequency_total_global
                                                                + self.frequency_excitation[i,j] )                        
                        for j in range(self.dissociation_types[i]):
                                self.frequency_total_global = ( self.frequency_total_global
                                                                + self.frequency_dissociation[i,j] )

                # Calculate the maximum value of the total collision frequency over all velocities
                #   f_max = max{f_tot(v), v}
                self.frequency_total_global_max = self.frequency_total_global.max()

                # Calculate null collision frequency (required to apply the null collision method)
                #   f_null(v) = f_max - f_tot(v)
                self.frequency_null_collision = self.frequency_total_global_max - self.frequency_total_global

                # Calculate relative collision probabilities
                # prel_i(v) = f_i(v) / f_max
                self.p_relative_null_collision = self.frequency_null_collision / self.frequency_total_global_max
                for i in range(self.types):
                        self.p_relative_elastic[i]    = self.frequency_elastic[i]    / self.frequency_total_global_max
                        self.p_relative_ionization[i] = self.frequency_ionization[i] / self.frequency_total_global_max
                        for j in range(self.excitation_types[i]):
                                self.p_relative_excitation[i,j] = ( self.frequency_excitation[i,j]
                                                                    / self.frequency_total_global_max )
                        for j in range(self.dissociation_types[i]):
                                self.p_relative_dissociation[i,j] = ( self.frequency_dissociation[i,j]
                                                                      / self.frequency_total_global_max )
                                
                # Calculate probability boundaries                        
                self.probability_limits[0]    = self.p_relative_elastic[0]
                self.probability_processes[0] = 'el->'+ self.names[0] + ': elastic' 
                for i in range(1, self.types):
                        self.probability_limits[i] = ( self.probability_limits[i-1]
                                                       + self.p_relative_elastic[i] )
                        self.probability_processes[i] = 'el->'+ self.names[i] + ': elastic' 
                for i in range(self.types):
                        self.probability_limits[self.types+i] = ( self.probability_limits[self.types+i-1]
                                                                  + self.p_relative_ionization[i] )
                        self.probability_processes[self.types+i] = 'el->'+ self.names[i] + ': ionization'                       
                self.probability_index_exc = 2 * self.types
                i = self.probability_index_exc  # self.types*2
                
                for j in range(self.types):
                        for k in range(self.excitation_types[j]):
                                self.probability_limits[i] = ( self.probability_limits[i-1]
                                                               + self.p_relative_excitation[j,k] )
                                self.probability_processes[i] = 'el->'+ self.names[j] + ': excitation#' + str(k)                               
                                i = i + 1
                self.probability_index_diss = i
                for j in range(self.types):
                        for k in range(self.dissociation_types[j]):
                                self.probability_limits[i] = ( self.probability_limits[i-1]
                                                               + self.p_relative_dissociation[j,k] )
                                self.probability_processes[i] = 'el->'+ self.names[j] + ': dissociation#' + str(k)
                                i = i + 1                                
                self.n_limits = i
                if debug:
                        print('Probability boundaries for electron scattering')
                        for i in range(self.n_limits):
                                print(i, self.probability_processes[i])
        
                return (status, message)


        def plot_xsec_electrons(self, plot_single=True, plot_total=False,
                                plot_frequencies=False, plot_relative=False, plot_boundaries=False,
                                dot_points=False):
                """ Plots cross sections for electrons collisions and related main statistical parameters, 
                    as a function of electron energy.

                        Plots, for the different processes (elastic scattering, ionization, dissociation types): 
                                - total cross sections
                                - collision frequencies
                                - relative scattering probabilities
                                - probability boundaries

                        Parameters
                        ----------

                        plot_single:      if True, plot the cross section values of single processes as a function of energy
                        plot_total:       if True, plot the total cross section values as a function of energy
                        plot_frequencies: if True, plot the collision frequency values as a function of energy
                        plot_relative:    if True, plot the relative scattering probabilities
                        plot_boundaries:  if True, plot the probability boundaries
                        dot_points:       if True, plot using dots instead of markers

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking

                if (not gnuplot_installed): return

                if dot_points: plot_style = 'dots'
                else:          plot_style = 'points'

                # Plot cross-sections for single processes
                if plot_single:
                        data_list = []
                        # Prepare graph
                        self.cross_section_graph = new_plot(title='Cross sections for electron scattering',
                                                            xlabel='Energy / eV',
                                                            ylabel='Cross section / m**2',
                                                            logx=True,
                                                            logy=True,
                                                            format_x='%g',
                                                            format_y='%g',
                                                            grid=True,
                                                            style=plot_style,
                                                            persistence=XSEC_PLOT_PERSIST,
                                                            redirect_output=True)

                        # Plot cross sections for each process type and gas type
                        for i in range(self.types):
                                data_list.append( [ self.sigma_energies(),
                                                         self.sigma_elastic[i],
                                                         'elastic ' + self.names[i],
                                                         None ] )
                        for i in range(self.types):
                                data_list.append( [ self.sigma_energies(),
                                                         self.sigma_ionization[i],
                                                         'ionization ' + self.names[i],
                                                         None ] )
                        for i in range(self.types):
                                for j in range(self.excitation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.sigma_excitation[i,j],
                                                                 'excitation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )                        
                        for i in range(self.types):
                                for j in range(self.dissociation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.sigma_dissociation[i,j],
                                                                 'dissociation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )

                        plot_curves(self.cross_section_graph, data_list)

                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.cross_section_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.cross_section_graph)

                
                # Plot total cross-sections
                if plot_total:
                        data_list = []
                        # Prepare graph
                        self.cross_section_graph = new_plot(title='Total cross sections for electron scattering',
                                                            xlabel='Energy / eV',
                                                            ylabel='Cross section / m**2',
                                                            logx=True,
                                                            logy=True,
                                                            format_x='%g',
                                                            format_y='%g',
                                                            grid=True,
                                                            style=plot_style,
                                                            persistence=XSEC_PLOT_PERSIST,
                                                            redirect_output=True)
                        # Plot total global cross section
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_global, 
                                                 'total global',
                                                 None ] )
                        # Plot total cross sections for each process type
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_elastic,
                                                 'total elastic',
                                                 None ] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_ionization,
                                                 'total ionization',
                                                 None] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_excitation,
                                                 'total excitation',
                                                 None ] )                        
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_dissociation,
                                                 'total dissociation',
                                                 None ] )

                        plot_curves(self.cross_section_graph, data_list)
                        
                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.cross_section_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.cross_section_graph)
                
                # Plot collision frequencies                
                if plot_frequencies:
                        data_list = []
                        self.coll_freq_graph = new_plot(title='Collision frequencies for electron scattering',
                                                        xlabel='Energy / eV',
                                                        ylabel='Frequency / Hz',
                                                        format_x='%g',
                                                        format_y='%g',
                                                        grid=True,
                                                        style=plot_style,
                                                        persistence=XSEC_PLOT_PERSIST,
                                                        redirect_output=True)                        
                        for i in range(self.types):
                                data_list.append( [ self.sigma_energies(),
                                                         self.frequency_elastic[i],
                                                         'elastic ' + self.names[i],
                                                         None ] )
                                data_list.append( [ self.sigma_energies(),
                                                         self.frequency_ionization[i],
                                                         'ionization ' + self.names[i],
                                                         None ] )
                                for j in range(self.excitation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.frequency_excitation[i,j],
                                                                 'excitation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )                                
                                for j in range(self.dissociation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.frequency_dissociation[i,j],
                                                                 'dissociation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.frequency_total_global,
                                                 'total',
                                                 None ] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.frequency_null_collision,
                                                 'null collision',
                                                 None ] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.frequency_total_global_max*numpy.ones(self.n_sigma),
                                                 'max',
                                                 None ] )

                        plot_curves(self.coll_freq_graph, data_list)
                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.coll_freq_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.coll_freq_graph)
                        
                # Plot relative collision frequencies
                if plot_relative:
                        data_list = []
                        self.rel_coll_freq_graph = new_plot(title='Relative collision frequencies for electron scattering',
                                                            xlabel='Energy / eV',
                                                            ylabel='Relative probability / %',
                                                            logx=True,
                                                            grid=True,
                                                            style=plot_style,
                                                            persistence=XSEC_PLOT_PERSIST,
                                                            redirect_output=True)                        
                        for i in range(self.types):
                                data_list.append( [ self.sigma_energies(),
                                                         self.p_relative_elastic[i]*100,
                                                         'elastic '+self.names[i],
                                                         None ] )
                                data_list.append( [ self.sigma_energies(),
                                                         self.p_relative_ionization[i]*100,
                                                         'ionization '+self.names[i],
                                                         None ] )
                                for j in range(self.excitation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.p_relative_excitation[i,j]*100,
                                                                 'excitation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )                                
                                for j in range(self.dissociation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.p_relative_dissociation[i,j]*100,
                                                                 'dissociation ' + self.names[i] + ' # ' + str(j),
                                                                 None ] )
                        data_list.append( [ self.sigma_energies(),
                                                 self.p_relative_null_collision*100,
                                                 'null collision',
                                                 None ] )

                        plot_curves(self.rel_coll_freq_graph, data_list)
                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.rel_coll_freq_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.rel_coll_freq_graph)
                        
                # Plot probability boundaries
                if plot_boundaries:
                        data_list = []
                        self.boundaries_graph = new_plot(title='Probability boundaries for electron scattering',
                                                         xlabel='Energy / eV',
                                                         ylabel='Probability limit',
                                                         logx=True,
                                                         grid=True,
                                                         style=plot_style,
                                                         persistence=XSEC_PLOT_PERSIST,
                                                         redirect_output=True)
                        for i in range(self.n_limits):
                                data_list.append( [ self.sigma_energies(),
                                                         self.probability_limits[i],
                                                         #'boundary #'+str(i),
                                                         self.probability_processes[i],
                                                         None ] )

                        plot_curves(self.boundaries_graph, data_list)
                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.boundaries_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.boundaries_graph)
                        
                return status, message


        # +-------------------------------------------------------------+
        # | Methods to read cross sections for ions/neutrals scattering |
        # +-------------------------------------------------------------+

        def read_xsec_ions_elastic(self, filename, sep, ion_index, target_index, plot=False):
                """Reads the table of cross sections values for ions-neutrals elastic scattering as a function 
                   of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                        
                        sep:            character used to separate values (usually the tab char)
                        ion_index:      integer value indentifying the type of ion (i.e. the type of neutral that was ionized) 
                                        that elastically scatters with the neutral
                        target_index:   integer value indentifying the type of neutral target
                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for elastic scattering of ' + self.names[ion_index] + '+ on ' + self.names[target_index]
                qm_ratio = -ELECTRON_CHARGE / (self.mass[ion_index] * ATOMIC_UNIT_MASS)

                (status, message) = self._read_xsec(filename, sep, qm_ratio, target_index, ion_index, plot, title)

                if (status == 0): self.sigma_elastic_ions[ion_index, target_index] = self.sigma_return

                return (status, message)


        def read_xsec_ions_charge_exchange(self, filename, sep, ion_index, target_index, plot=False):
                """Reads the table of cross sections values for ions-neutrals charge exchange scattering 
                   as a function of energy from a ascii file.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                        
                        sep:            character used to separate values (usually the tab char)
                        ion_index:      integer value indentifying the type of ion (i.e. the type of neutral that was ionized) 
                                        that scatters with the neutral
                        target_index:   integer value indentifying the type of neutral target


                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                title    = 'Cross section for charge exchange scattering of ' + self.names[ion_index] + '+ on ' + \
                           self.names[target_index] 
                qm_ratio = -ELECTRON_CHARGE / (self.mass[ion_index] * ATOMIC_UNIT_MASS)

                (status, message) = self._read_xsec(filename, sep, qm_ratio, target_index, ion_index, plot, title)

                if (status == 0): self.sigma_charge_exchange_ions[ion_index,target_index] = self.sigma_return

                return (status, message)


        def calculate_total_xsec_ions(self, debug=False):
                """Calculates main statistical parameters from cross sections for ion impact processes.

                        Calculate, for the different processes (elastic scattering, ionization, dissociation types): 
                                - total cross sections
                                - collision frequencies
                                - relative scattering probabilities
                                - probability boundaries
                        (see the module docstring for more information on these quantities and their usage)

                        Parameters
                        ----------

                        Initialized data attributes
                        ---------------------------

                        Calculated total global cross section for each ion, sum of all targets and all processes
                        self.sigma_total_global_ions:          sum of all cross section for all processes

                        Calculated cross sections for each ion/neutral combination, sum of all processes
                        self.sigma_total_ions:                 sum of all cross sections (elastic + charge exchange)
                                                                 for each ion/target combinations

                        Calculated collision frequencies for   each ion/neutral combination
                        self.frequency_elastic_ions:           collision frequencies for elastic scattering for each 
                        self.frequency_charge_exchange_ions:   collision frequencies for charge exchange sc
                        self.frequency_total_ions:             total collision frequency (elastic + charge exchange) 

                        Calculated collision frequencies for   each ion (over all neutral types)
                        self.frequency_total_global_ions:      sum of all collision frequencies, for all possible processes
                                                               and all neutral types 
                        self.frequency_total_global_max_ions:  max value of the total global frequency, over all speeds
                        self.frequency_null_collision_ions:    null collision frequency, defined as the difference between
                                                                 the maximum value (over all velocities) of the 
                                                                 total global collision frequency
                                                                 and the total global collision frequency at that speed

                        Calculated probabily limits for elastic, charge exchange and null collision
                        self-p_null_ions:                      probability to have e null collision for that ion (all targets)
                        self.probability_limits_ions:          once stated which is the target, this are the probability limits
                                                               for the possible collision types: elastic, charge exchange, null

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                self.sigma_total_global_ions            = numpy.zeros( (self.types, self.n_sigma_ions), 'd', order='F')

                self.sigma_total_ions                   = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')

                self.frequency_elastic_ions             = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')
                self.frequency_charge_exchange_ions     = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')
                self.frequency_total_ions               = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')
                
                self.frequency_total_global_ions        = numpy.zeros( (self.types, self.n_sigma_ions), 'd', order='F')
                self.frequency_total_global_max_ions    = numpy.zeros( self.types, 'd')                
                self.frequency_null_collision_ions      = numpy.zeros( (self.types, self.n_sigma_ions), 'd', order='F')

                self.p_relative_null_ions               = numpy.zeros( (self.types, self.n_sigma_ions), 'd', order='F')
                self.p_relative_elastic_ions            = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')
                self.p_relative_charge_exchange_ions    = numpy.zeros( (self.types, self.types, self.n_sigma_ions), 'd', order='F')

                self.probability_limits_ions            = numpy.zeros( (self.types, self.types, 2, self.n_sigma_ions), 'd', order='F')

                # Calculate total cross section (elastic + charge exchange) for each ion/target combination
                for ion in range(self.types):
                        for target in range(self.types):
                                self.sigma_total_ions[ion,target] = (   self.sigma_elastic_ions[ion,target]
                                                                      + self.sigma_charge_exchange_ions[ion,target] )

                # Calculate collision frequencies for all processes and their sum
                # f_(i,t)(v) = n_t * s_(i,t)(v) * v_rel
                for ion in range(self.types):
                        for target in range(self.types):
                                self.frequency_elastic_ions[ion,target]         = (  self.number_density[target]
                                                                                   * self.sigma_elastic_ions[ion,target]
                                                                                   * self.sigma_velocities_ions[ion] )
                                self.frequency_charge_exchange_ions[ion,target] = (  self.number_density[target]
                                                                                   * self.sigma_charge_exchange_ions[ion,target]
                                                                                   * self.sigma_velocities_ions[ion] )
                                self.frequency_total_ions[ion,target] = (   self.frequency_elastic_ions[ion,target]
                                                                          + self.frequency_charge_exchange_ions[ion,target] )
                                                               
                # Calculate total global collision frequency for each ion type, summing over all neutral types
                # f_tot_i(v) = sum{f_(i,t)(v), t}
                for ion in range(self.types):
                        for target in range(self.types):
                                self.frequency_total_global_ions[ion] = (  self.frequency_total_global_ions[ion]
                                                                         + self.frequency_total_ions[ion,target] )

                # Calculate the maximum value of the total global collision frequency over all velocities
                #   f_max = max{f_tot(v), v}
                for ion in range(self.types):
                        self.frequency_total_global_max_ions[ion] = self.frequency_total_global_ions[ion].max()

                # Calculate the null collision frequency
                #   f_null(v) = f_max - f_tot(v)
                for ion in range(self.types):
                        self.frequency_null_collision_ions[ion] = (  self.frequency_total_global_max_ions[ion] 
                                                                   - self.frequency_total_global_ions[ion]    )
                # Calculate relative scattering probabilities
                for ion in range(self.types):
                        # Probability to have a null collision (regardless of the target)
                        self.p_relative_null_ions[ion] = self.frequency_null_collision_ions[ion] / self.frequency_total_global_max_ions[ion]
                        for target in range(self.types):
                                # Relative probability to have elastic or charge exchange, for a specific target,
                                # if the collision is not null
                                self.p_relative_elastic_ions[ion, target] = (  self.frequency_elastic_ions[ion,target]
                                                                             / self.frequency_total_ions[ion,target]   )
                                self.p_relative_charge_exchange_ions[ion, target] = (  self.frequency_charge_exchange_ions[ion,target]
                                                                                     / self.frequency_total_ions[ion,target]   )
                
                # Probability boundaries to have a collision elastic, charge exchange or null
                # once the target has been identified        
                for ion in range(self.types):
                        for target in range(self.types):
                                self.probability_limits_ions[ion, target, 0] = (  (1 - self.p_relative_null_ions[ion])
                                                                                * self.p_relative_elastic_ions[ion,target]  )
                                self.probability_limits_ions[ion, target, 1] = (  self.probability_limits_ions[ion, target, 0]
                                                                                + (1 - self.p_relative_null_ions[ion])
                                                                                * self.p_relative_charge_exchange_ions[ion,target]  )
                return (status, message)


        def plot_xsec_ions(self, ion_type, plot_total=True, plot_frequencies=False, plot_relative=False, plot_boundaries=False,
                           dot_points=False):
                """Plots cross sections for ion collisions and related main statistical parameters, as a function of electron energy.

                        Plots, for the different processes (elastic scattering, ionization, dissociation types): 
                                - total cross sections
                                - collision frequencies
                                - relative probability of elastic scattering

                        Parameters
                        ----------

                        ion_type:         type of ion for which the quantities will be plotted
                        plot_total:       if True, plot the total cross section values as a function of energy
                        plot_frequencies: if True, plot the collision frequency values as a function of energy
                        plot_relative:    if True, plot the relative probability of elastic scattering
                        plot_boundaries:  if True, plot the probavility boundaries

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking

                if (not gnuplot_installed): return

                if dot_points: plot_style = 'dots'
                else:          plot_style = 'points'               

                # Plot total cross-sections                       
                if plot_total:
                        data_list = []
                        # Prepare graph
                        self.cross_section_graph = new_plot(title='Cross sections for ' + self.names[ion_type] + '+ impact',
                                                            xlabel='Energy / eV',
                                                            ylabel='Cross section / m**2',
                                                            logx=True,
                                                            logy=True,
                                                            grid=True,
                                                            style=plot_style,
                                                            persistence=XSEC_PLOT_PERSIST,
                                                            redirect_output=True)
                        # Plot total global cross section (all neutral types and all processes)
                        data_list.append( [ self.sigma_energies_ions(ion_type),
                                                 self.sigma_total_global_ions[ion_type],
                                                 'total global' + self.names[ion_type] + '+',
                                                 None ] )
                        # Plot cross sections for each neutral and process type
                        for target in range(self.types):
                                # Toal cross section (all processes)
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                         self.sigma_total_ions[ion_type,target],
                                                         'total ' + self.names[ion_type] + '+ on ' + self.names[target],
                                                         None ] )
                                # Elastic cross section
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                         self.sigma_elastic_ions[ion_type,target],
                                                         'elastic ' + self.names[ion_type] + '+ on ' + self.names[target],
                                                         None ] )
                                # Charge exchange cross section
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                         self.sigma_charge_exchange_ions[ion_type,target],
                                                         'charge exchange '+self.names[ion_type]+'+ on '+ self.names[target],
                                                         None ] )
                                
                        plot_curves(self.cross_section_graph, data_list)
                        
                        plot_close(self.cross_section_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.cross_section_graph)
                                
                # Plot collision frequencies
                if plot_frequencies:
                        data_list = []
                        self.coll_freq_graph = new_plot(title='Collision frequencies for ' + self.names[ion_type] + '+ impact',
                                                        xlabel='Energy / eV',
                                                        ylabel='Frequency / Hz',
                                                        logx=True,
                                                        logy=True,
                                                        format_x='%g',
                                                        format_y='%g',
                                                        grid=True,
                                                        style=plot_style,
                                                        persistence=XSEC_PLOT_PERSIST,
                                                        redirect_output=True)
                        
                        for target in range(self.types):
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.frequency_elastic_ions[ion_type,target],
                                                    'elastic ' + self.names[ion_type] + '+ on ' + self.names[target],
                                                    None ] )
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.frequency_charge_exchange_ions[ion_type,target],
                                                    'charge exchange ' + self.names[ion_type] + '+ on ' + self.names[target],
                                                    None ] )
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.frequency_total_ions[ion_type, target],
                                                    'total ' + self.names[ion_type] + '+ on ' + self.names[target],
                                                    None ] )
                        data_list.append( [ self.sigma_energies_ions(ion_type),
                                            self.frequency_total_global_ions[ion_type],
                                            'total global ' + self.names[ion_type] + '+',
                                            None ] )                                
                        data_list.append( [ self.sigma_energies_ions(ion_type),
                                            self.frequency_null_collision_ions[ion_type],
                                            'null collision ' + self.names[ion_type] + '+', None ] )
                        data_list.append( [ self.sigma_energies_ions(ion_type),
                                            self.frequency_total_global_max_ions[ion_type]*numpy.ones(self.n_sigma_ions),
                                            'max ' + self.names[ion_type] + '+', None ] )

                        plot_curves(self.coll_freq_graph, data_list)                        

                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.coll_freq_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.coll_freq_graph)

                # Plot relative collision frequencies
                if plot_relative:
                        data_list = []
                        self.cross_sect_rel_graph = new_plot(title='Relative probability of ion scattering for '
                                                             + self.names[ion_type]
                                                             + '+ impact',
                                                             xlabel='Energy / eV',
                                                             ylabel='Relative probability / %',
                                                             logx=True,
                                                             grid=True,
                                                             ymin = 0,
                                                             ymax = 100,
                                                             style=plot_style,
                                                             persistence=XSEC_PLOT_PERSIST,
                                                             redirect_output=True)
                        data_list.append( [ self.sigma_energies_ions(ion_type),
                                            self.p_relative_null_ions[ion_type]*100,
                                            'null '+self.names[ion_type]+'+',
                                            None ] )                        
                        for target in range(self.types):
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.p_relative_elastic_ions[ion_type,target]*100,
                                                    'elastic '+self.names[ion_type]+'+ on '+self.names[target],
                                                    None ] )
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.p_relative_charge_exchange_ions[ion_type,target]*100,
                                                    'charge exchange '+self.names[ion_type]+'+ on '+self.names[target],
                                                    None ] )                                

                        plot_curves(self.cross_sect_rel_graph, data_list)
                        
                        # Shut down gnuplot process, window remain visible due to persist option
                        plot_close(self.cross_sect_rel_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.cross_sect_rel_graph)

                # Plot probability boundaries
                if plot_boundaries:
                        for target in range(self.types):                        
                                data_list = []
                                self.cross_sect_lim_graph = new_plot(title='Probability boundaries for '
                                                                     + self.names[ion_type]
                                                                     + '+ impact on ' + self.names[target],
                                                                     xlabel='Energy / eV',
                                                                     ylabel='Probability limit',
                                                                     logx=True,
                                                                     grid=True,
                                                                     ymin = 0,
                                                                     ymax = 1,
                                                                     style=plot_style,
                                                                     persistence=XSEC_PLOT_PERSIST,
                                                                     redirect_output=True)

                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.probability_limits_ions[ion_type,target,0],
                                                    'limit 0 (elastic) '+self.names[ion_type]+'+ on '+self.names[target],
                                                    None ] )
                                data_list.append( [ self.sigma_energies_ions(ion_type),
                                                    self.probability_limits_ions[ion_type,target,1],
                                                    'limit 1 (charge exchange) '+self.names[ion_type]+'+ on '+self.names[target],
                                                    None ] )                                

                                plot_curves(self.cross_sect_lim_graph, data_list)
                        
                                # Shut down gnuplot process, window remain visible due to persist option
                                plot_close(self.cross_sect_lim_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                                del(self.cross_sect_lim_graph)                        

                return status, message


        # +----------------------------------------------------------------------+
        # | Methods to read cross sections for electrons recombination with ions |
        # +----------------------------------------------------------------------+

        def read_xsec_ele_ion_recomb(self, filename, sep, ion_index, rec_type=None, diss_type=None, plot=False):
                """Reads the table of cross section values for recombination of electrons with ions.

                        Parameters
                        ----------

                        filename:       name of file containing cross section values at energies in the form

                                                #Comment
                                                <energy> <separator> <value>    #Comment
                                                <energy> <separator> <value>    #Comment

                                        where energies are in eV and cross section values are in m**2
                                        
                        sep:            character used to separate values (usually the tab character)
                        ion_index:      integer value indentifying the type of ion (must be a molecule) 
                                        with which the elctron recombine
                        rec_type:       type of recombination process
                                        'dissociative' -> dissociative recombination of e- with a molecular ion
                                        '3body'        -> recombination of e- with ion, with transfer of the released energy 
                                                          to another e-
                                        'radiative'    -> recombination of e- with ion, with radiative release of excess energy

                        diss_type:      integer value indentifying the type of dissociative recombination process to which
                                        the cross section is referred (only for dissociative recombination processes)
                                        (the number of dissociation types for each molecule type is defined 
                                        in the "neutrals.csv" file)
                        plot:           if True, plot the cross section values as a function of energy,
                                        before and after interpolation


                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = negative energy value found
                                         2 = negative cross-section value found
                                         3 = ion index is outside allowed bounds
                                         4 = unspecified recombination type
                                         5 = unknown recombination type
                                         6 = dissociation process index is missing
                                         7 = dissociation process index is out of bounds
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                if ( (ion_index < 0) or (ion_index > self.types-1)):
                        status = 3
                        message = 'ion type ' + str(ion_index) + ' is out of bounds'
                        return(status, message)
                if (rec_type is None):
                        status = 4
                        message = 'unspecified recombination type'
                        return(status, message)
                if (rec_type not in ('dissociative', '3body', 'radiative')):
                        status = 5
                        message = 'unknown recombination type \"' + str(rec_type) + '\"'
                        return(status, message)
                if (rec_type == 'dissociative'):
                        if (diss_type is None):
                                status = 6
                                message = 'dissociation process index is missing'
                                return (status, message)
                        if (diss_type > self.dissociation_types[ion_index]):
                                status = 7
                                message = 'dissociation process index ' + str(diss_type) + ' is out of bounds'
                                return (status, message)

                title    = 'Cross section for ' + rec_type + ' recombination of e- with ' + self.names[ion_index]
                if (rec_type =='dissociative'): 
                        title = title + ' (process #' + str(diss_type) + ')'
                qm_ratio = - ELECTRON_CHARGE / ELECTRON_MASS

                (status, message) = self._read_xsec(filename, sep, qm_ratio, ion_index, None, plot, title)

                if (status == 0): 
                        if (rec_type == 'dissociative'): self.sigma_recombination_diss[ion_index,diss_type] = self.sigma_return
                        if (rec_type == '3body'):        self.sigma_recombination_3body[ion_index]          = self.sigma_return
                        if (rec_type == 'radiative'):    self.sigma_recombination_rad[ion_index]            = self.sigma_return

                return (status, message)


        def calculate_total_xsec_ele_ion_recomb(self):
                """Calculates main statistical parameters from cross sections.

                        Calculate, for the different processes of electron-ion recombination:
                                - total cross sections
                                - recombination rate coefficients (non averaged), i.e. product cross section by electron velocity
                        (see the module docstring for more information on these quantities and their usage)

                        When calculating the recombination process for an electron with a given ion type, 
                        the background density involved in the calculation is the density of that ion type.
                        Since it is not constant [*], it is not possible to pre-calculate 
                        recombination frequencies, nor recombination rates. They must be calculated during simulation
                        at each timestep based on the ion density at that time.

                        [*] Instead, neutral density can be considered appoximately constant, if the ionization degree is very low
                        so collision frequencies with neutrals can be pre-calculated.

                        Parameters
                        ----------

                        Initialized data attributes
                        ---------------------------

                        self.ratecoeff_recomb_diss:     rate coefficients (sigma*v) for dissociative recombination processes 
                                                             of each ion type

                        self.sigma_total_recomb_diss:   sum of dissociative recombination cross sections for all neutral types
                                                             and all dissociation processes for each neutral type

                        self.sigma_total_global_recomb: sum of all cross section for all recombination processes

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking      

                self.sigma_total_recombination_diss    = numpy.zeros( self.n_sigma, 'd' )
#               self.sigma_total_recombination_3body   = numpy.zeros( self.n_sigma, 'd' )
#               self.sigma_total_recombination_rad     = numpy.zeros( self.n_sigma, 'd' )
                self.sigma_total_global_recombination  = numpy.zeros( self.n_sigma, 'd' )
                self.ratecoeff_recomb_diss             = numpy.zeros( (self.types, MAX_DISSOCIATION_TYPES, self.n_sigma), 
                                                                      'd', order = 'F')
                # Calculate total cross sections for each recombination process
                # (at the moment only dissociative recombination is considered)
                for ion_index in range(self.types):
                        for diss_type in range(self.dissociation_types[ion_index]):
                                self.sigma_total_recombination_diss = self.sigma_total_recombination_diss + \
                                                                      self.sigma_recombination_diss[ion_index,diss_type]

                # Calculate total global cross-section, for all the possible recombination types
                # (at the moment only dissociative recombination is considered,
                # so it's equal to total cross section for that type of recombination)
                self.sigma_total_global_recombination  = self.sigma_total_recombination_diss
#                                                       + self.sigma_total_recombination_3body
#                                                       + self.sigma_total_recombination_rad

                # Calculate rate coefficients (not averaged) for each recombination process
                # (at the moment only dissociative recombination is considered)
                # rate coefficient is the product cross section by velocity -> sigma(v) * v
                for ion_index in range(self.types):
                        for diss_type in range(self.dissociation_types[ion_index]):
                                self.ratecoeff_recomb_diss[ion_index,diss_type]= self.sigma_recombination_diss[ion_index,diss_type]*\
                                                                                 self.sigma_velocities
#                        self.ratecoeff_recomb_3body[ion_index] = self.sigma_recombination_3body[ion_index] * self.sigma_velocities     
#                        self.ratecoeff_recomb_rad[ion_index]   = self.sigma_recombination_rad[ion_index]   * self.sigma_velocities

                return (status, message)


        def plot_xsec_ele_ion_recomb(self, plot_total=False, plot_rates=False, dot_points=False):
                """Plots cross sections for electrons collisions and related main statistical parameters, 
                   as a function of electron energy.

                        Plots, for the different processes (elastic scattering, ionization, dissociation types): 
                                - total cross sections
                                - collision rates 

                        Parameters
                        ----------

                        plot_total:       if True, plot the total cross section values as a function of energy
                        plot_rates:       if True, plot the collision rate values as a function of energy

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = error
                                message: a string containing an error message or 'Ok'
                """

                status, message  = 0, 'OK'              # This variable is used for error checking

                if (not gnuplot_installed): return

                if dot_points: plot_style = 'dots'
                else:          plot_style = 'points'                

                # Plot total cross-sections
                if plot_total:
                        data_list = []
                        self.el_ion_recomb_graph = new_plot(title='Cross sections for e-/ion recombination',
                                                            xlabel='Energy / eV',
                                                            ylabel='Cross section / m**2',
                                                            logx=True,
                                                            logy=True,
                                                            format_x='%g',
                                                            format_y='%g',
                                                            grid=True,
                                                            style=plot_style,
                                                            persistence=XSEC_PLOT_PERSIST,
                                                            redirect_output=True)                                                      
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_global_recombination, 
                                                 'total global',
                                                 None ] )

                        # Plot total cross sections for each process type
                        data_list.append( [ self.sigma_energies(),
                                                 self.sigma_total_recombination_diss, 
                                                 'total dissociative recombination',
                                                 None ] )
                        # Plot cross sections for each process type and gas type
                        for i in range(self.types):
                                for j in range(self.dissociation_types[i]):
                                        data_list.append( [ self.sigma_energies(),
                                                                 self.sigma_recombination_diss[i,j],
                                                                 'dissociative recombination ' + self.names[i] + ' # '+str(j),
                                                                 None ] )
                                        
                        plot_curves(self.el_ion_recomb_graph, data_list)
                                         
                        plot_close(self.el_ion_recomb_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.el_ion_recomb_graph)

                # Plot collision rates
                if plot_rates:
                        data_list = []
                        self.coll_rates_graph = new_plot(title='Collision rate coefficients for e-/ion recombination',
                                                         xlabel='Energy / eV',
                                                         ylabel='Collision rate coefficients / m**3 s**-1',
                                                         logx=True,
                                                         logy=True,
                                                         format_x='%g',
                                                         format_y='%g',
                                                         grid=True,
                                                         style=plot_style,
                                                         persistence=XSEC_PLOT_PERSIST,
                                                         redirect_output=True)                                                       
                        for i in range(self.types):
                                for j in range(self.dissociation_types[i]):
                                        self.plot_lsit.append( [ self.sigma_energies(),
                                                                 self.ratecoeff_recomb_diss[i,j],
                                                                 ( 'dissociative recombination '
                                                                   + self.names[i]
                                                                   + ' # '
                                                                   + str(j) ),
                                                                 None ] )

                        plot_curves(self.coll_rates_graph, data_list)

                        plot_close(self.coll_rates_graph, keep_output=True, delay=DEL_DATA_FILES_DELAY)
                        del(self.coll_rates_grap)

                return status, message


        def initialize_savefile(self, filename_neutral, append=False, sep='\t', ext='.csv'):
                """Initializes the file to which simulation data about neutral particles will be saved

                   Parameters
                   ----------
                   filename_neutral:      name of the file to which the simulation quantities will be saved
                   sep:                   character or sequence to separate data columns
                   ext:                   extension to give to filename (deafult is '.csv')

                   Initialized data attributes
                   ---------------------------

                   append:                if True, open the file for appendig new data to the end
                                          header will not be re-written
                   self.sep:              character or sequence to separate data columns
                   self.f_neutral:        file to which values of several quantities will be saved
                """

                # If there are not target molecules, do nothing
                if (self.types_molecules == 0): return
                
                self.sep = sep
                self.filename = str(filename_neutral) + ext
                
                # Open data file where to save time evolution of mean values
                if append:
                        data_file = open(self.filename,'a')
                else:
                        data_file = open(self.filename,'w')
                        # Write header
                        data_file.write('\"t[ns]\"' + self.sep )                
                        for i in range(self.types):
                                for j in range(self.dissociation_types[i]):
                                        data_file.write('\"Diss_rate ' + self.names[i] + '#'+str(j) + \
                                                          ' [m**-3*s**-1]\"' + self.sep   )
                                        data_file.write('\"Diss_rate const ' + self.names[i] + '#'+str(j) + \
                                                          ' [m**3*s**-1]\"' + self.sep   )
                        data_file.write(EOL)
                data_file.close()


        def save_data_to_files(self, time):
                """Saves actual data values of dissociation rates to files """

                # If there are not target molecules, do nothing
                if (self.types_molecules == 0): return
                
                data_file = open(self.filename,'a')
                data_file.write( str(1E9*time) + self.sep )                
                for i in range(self.types):
                        for j in range(self.dissociation_types[i]):
                                data_file.write( str(self.dissociation_rate[i][j]) + self.sep )
                                data_file.write( str(self.dissociation_rate_constant[i][j]) + self.sep )
                data_file.write(EOL)
                data_file.close()
