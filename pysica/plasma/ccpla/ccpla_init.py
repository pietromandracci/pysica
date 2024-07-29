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

""" Simulation of a capacitively coupled plasma discharge: data initialization """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

# Mudules from the standard Python library
import os
import shutil
import math
from random import random
from optparse import OptionParser

# Modules provided by the Python community
import numpy

# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.functions.random_pdf import random_maxwell_velocity
from pysica.plasma.ccpla.ccpla_defaults import *

# Import required modules, classes, and functions
from pysica.managers.io.io_screen import wait_input
from pysica.managers import data_manager, unit_manager
from pysica.functions.physics import pressure_conversion, number_density


def initialize_parameters(parameters, verbose=False, saveonly=False, filename_config='', filename_defaults='', restricted=False):
    """  Read the simulation parameters from a configuration file

         parameters:        an istance of the ConfigurationOptions class, defined in ccpla_defaults.py
         verbose:           if True, write some info to the console
         saveonly:          if True, save default values of the parameters to a file and then return
         restricted:        if True, some parameters are not read 
    """
    

    (status, message) = (0, OK)

    # +-------------------------+
    # | Read configuration file |
    # +-------------------------+        
        
    # Define a dictionary with variable names, their default values, descriptive strings, and allowed ranges
    configuration_data = data_manager.ConfigurationData()        
    configuration_data.d['Nmax_particles']    = [ parameters.Nmax_particles, 'max number of electrons/ions', 'checkrange', 1, 
                                                  NMAXPARTICLES ]
    configuration_data.d['rescale_factor']    = [ parameters.rescale_factor, 'rescale factor', 'checkrange', 1, 
                                                  MAX_RESCALE_FACTOR ]
    configuration_data.d['start_ion_deg']     = [ parameters.start_ion_deg, 'starting ionization degree', 
                                                  'checkrangestrictmin', 0.0, MAX_IONIZATION_DEGREE ]
    configuration_data.d['T_neutrals']        = [ parameters.T_neutrals, 'gas temperature', 'checkminstrict', 0 ]
    configuration_data.d['p_neutrals']        = [ parameters.p_neutrals, 'total pressure', 'checkminstrict', 0 ]
    configuration_data.d['distance']          = [ parameters.distance, 'electrodes distance', 'checkminstrict', 0 ]
    configuration_data.d['length']            = [ parameters.length, 'electrodes lateral length', 'checkminstrict', 0 ]
    configuration_data.d['V_bias']            = [ parameters.V_bias, 'electric bias', 'checkminstrict', 0 ]
    configuration_data.d['frequency']         = [ parameters.frequency, 'electric bias frequency', 'checkmin', 0 ]
    configuration_data.d['phase']             = [ parameters.phase, 'electric bias phase at t=0 (deg)', 'checkrange', 0, 180 ]
    configuration_data.d['N_cells']           = [ parameters.N_cells, 'number of cells in PIC scheme', 'checkrange',
                                                  NMINCELLS, NMAXCELLS ]
    configuration_data.d['duration']          = [ parameters.sim_duration, 'requested duration of the simulation',
                                                  'checkmin', 0.0,  ]
    configuration_data.d['dt_output']         = [ parameters.dt_output, 'time between data outputs', 'checkminstrict', 0 ]
    configuration_data.d['dt']                = [ parameters.dt, 'simulation timestep', 'checkmin', 0 ]
    configuration_data.d['save_delay']        = [ parameters.save_delay, 'data save periodicity', 'checkmin', 0 ]
    configuration_data.d['save_delay_dist']   = [ parameters.save_delay_dist, 'data save periodicity for distributions', 'checkmin', 1 ]
    configuration_data.d['maxcollfreq']       = [ parameters.maxcollfreq, 'maximum allowed collision frequency for variable dt',
                                                  'checkminstrict', 0 ]
    configuration_data.d['min_scattered']     = [ parameters.min_scattered,
                                                  'minimum required scattering events to apply many particles method',
                                                  'checkminstrict', 0 ]
    configuration_data.d['lateral_loss']      = [ 1*parameters.lateral_loss, 'particle lateral loss activation flag', 'boolean' ]
    configuration_data.d['isactive_recomb']   = [ 1*parameters.isactive_recomb, 'recombination processes activation flag',
                                                  'boolean' ]
    configuration_data.d['e_min_sigma']       = [ parameters.e_min_sigma, 'minimum cross section energy for electrons',
                                                  'checkmin', 0 ]
    configuration_data.d['e_max_sigma']       = [ parameters.e_max_sigma, 'maximum cross section energy for electrons',
                                                  'checkmin', 0 ]
    configuration_data.d['e_min_sigma_ions']  = [ parameters.e_min_sigma_ions, 'minimum cross section energy for ions',
                                                  'checkmin', 0 ]
    configuration_data.d['e_max_sigma_ions']  = [ parameters.e_max_sigma_ions, 'maximum cross section energy for ions',
                                                  'checkmin', 0 ]
    configuration_data.d['N_sigma']           = [ parameters.N_sigma,
                                                  'number of electrons cross section values to be interpolated',
                                                  'checkminstrict', 10 ]
    configuration_data.d['N_sigma_ions']      = [ parameters.N_sigma, 'number of ions cross section values to be interpolated',
                                                  'checkminstrict', 10 ]
    configuration_data.d['dot_points']        = [ parameters.dot_points, 'use dots instead of crosses in plots', 'boolean' ]       
    configuration_data.d['n_max_points']      = [ parameters.n_max_points, 'max number of points in historic plots',
                                                  'checkrange', 100, N_MAX_OUTPUT ]
    configuration_data.d['decimation_factor'] = [ parameters.decimation_factor, 'decimation factor for historic plots',
                                                  'checkrange', 2, MAX_DECIMATION_FACTOR ]        
       
    # If requested, only save values on a file and exit
    if saveonly:
        (status, message) = configuration_data.write_file(FILENAME_DEFAULTS)
        return (status, message)

    # Read variable values from configuration file
    if verbose:
        print('\nReading configuration parameters from file \"' + FILENAME_CONFIG + '\" ...')
    (status, message) = configuration_data.read_file(FILENAME_CONFIG, check_values=True)
    if (status != 0):
        message = 'Error in file \"' + FILENAME_CONFIG + '\": ' + message
        return (status, message)

    # +-----------------------------------------------------------------+
    # | Check consistency of configuration values and copy to variables |
    # +-----------------------------------------------------------------+

    if verbose: print('\nChecking parameters ...')

    ERROR = 'Error reading file \"' + FILENAME_CONFIG + '\": '

    # Number of particles
    if not restricted: parameters.Nmax_particles = int(configuration_data.d['Nmax_particles'][0])

    # Rescale factor
    if ( (configuration_data.d['rescale_factor'][0] > 1                           ) and \
         (configuration_data.d['rescale_factor'][0] > parameters.Nmax_particles/10)       ):
            status  = 1
            message = ERROR + 'rescale factor cannot exceed the max number of electrons divided by 10' + EOL
            return (status, message)
    parameters.rescale_factor = configuration_data.d['rescale_factor'][0]        
        
    # Starting ionization degree
    parameters.start_ion_deg = configuration_data.d['start_ion_deg'][0]

    # Temperature and total pressure of gases
    parameters.T_neutrals = configuration_data.d['T_neutrals'][0]
    parameters.p_neutrals = configuration_data.d['p_neutrals'][0]

    # Min and max cross section energies for electrons collisions
    if (configuration_data.d['e_min_sigma'][0] >= configuration_data.d['e_max_sigma'][0]):
        status  = 1
        message = ERROR + 'minimum energy must be lower than maximum energy' + EOL
        return (status, message)
    parameters.e_min_sigma = configuration_data.d['e_min_sigma'][0]
    parameters.e_max_sigma = configuration_data.d['e_max_sigma'][0]

    # Min and max cross section energies for ions collisions
    if (configuration_data.d['e_min_sigma_ions'][0] >= configuration_data.d['e_max_sigma_ions'][0]):
        status  = 1
        message = ERROR + 'minimum energy must be lower than maximum energy' + EOL
        return (status, message)
    parameters.e_min_sigma_ions = configuration_data.d['e_min_sigma_ions'][0]
    parameters.e_max_sigma_ions = configuration_data.d['e_max_sigma_ions'][0]

    # Number of cross-section values to be calculated by interpolation
    parameters.N_sigma      = int(configuration_data.d['N_sigma'][0])
    parameters.N_sigma_ions = int(configuration_data.d['N_sigma_ions'][0])

    # Activation of e/ion recombination processes
    if (configuration_data.d['isactive_recomb'][0] > 0):
        parameters.isactive_recomb = True
    else:
        parameters.isactive_recomb = False

    # Activation of e/ion lateral loss
    if (configuration_data.d['lateral_loss'][0] > 0):
        parameters.lateral_loss = True
    else:
        parameters.lateral_loss = False

    # Electrodes distance and area
    parameters.distance = configuration_data.d['distance'][0]
    parameters.length   = configuration_data.d['length'][0]

    # Electric bias
    parameters.V_bias = configuration_data.d['V_bias'][0]

    # Eletric field frequency
    parameters.frequency = configuration_data.d['frequency'][0]

    # Eletric field phase
    parameters.phase = configuration_data.d['phase'][0] / 180.0 * math.pi

    # Calculate electric field pulsation
    parameters.E_pulsation = 2*math.pi * parameters.frequency        

    # Cells for PIC scheme
    parameters.N_cells = int(configuration_data.d['N_cells'][0])

    # Simulation timestep
    if (configuration_data.d['dt'][0] == 0):
        parameters.dt_var = True
        parameters.dt     = 0
    else:
        # In case of variable electric field, check that dt is much lower than the period of electric field
        if (parameters.frequency != 0): 
            if (configuration_data.d['dt'][0] > 1/(MINTRATIO*parameters.frequency)):
                status  = 1
                message = ERROR + 'iteration time must not exceed 1 / ' \
                        + str(MINTRATIO) + ' of electric field period' + EOL
                return (status, message)
        parameters.dt_var = False                        
        parameters.dt     = configuration_data.d['dt'][0]
        
    # Time between data acquisitions
    if (configuration_data.d['dt_output'][0] < parameters.dt):
        status  = 1
        message = ERROR + 'data acquisition interval must be greater than or equal to timestep' + EOL
        return (status, message)
    parameters.dt_output = configuration_data.d['dt_output'][0]        

    # Requested simulation duration
    if ( (configuration_data.d['duration'][0] > 0 ) and
         (configuration_data.d['duration'][0] <= configuration_data.d['dt_output'][0]) ):
        status  = 1
        message = ERROR + 'requested simulation duration must be greater then data acquisition interval' + EOL
        return (status, message)
    parameters.sim_duration = configuration_data.d['duration'][0]

    # Save data to file periodicity
    parameters.save_delay = int(configuration_data.d['save_delay'][0])

    # Save distributions to file periodicity
    parameters.save_delay_dist = int(configuration_data.d['save_delay_dist'][0])
    
    # Maximum collision frequency allowed in calcolus of dt (for variable dt)
    parameters.maxcollfreq = configuration_data.d['maxcollfreq'][0]

    # Minimum number of collisions required to apply many perticles method
    parameters.min_scattered = configuration_data.d['min_scattered'][0]

    # Use points instead of symbols in plots
    parameters.dot_points = configuration_data.d['dot_points'][0]

    # Decimation parameters
    parameters.n_max_points      = int(configuration_data.d['n_max_points'][0])
    parameters.decimation_factor = int(configuration_data.d['decimation_factor'][0])
        
    # Release memory used by configuration_data instance
    del configuration_data


    # +----------------------------+
    # | Calculate other parameters |
    # +----------------------------+

    # Calculate the number of data values to be registered and plotted
#   parameters.N_output = NMAXOUTPUT
#   parameters.N_output = int(math.ceil(parameters.sim_duration/parameters.dt_output))
#   if (parameters.N_output > NMAXOUTPUT):
#      status  = 1
#      message = 'Number of data acquisitions is too high: reduce simulation duration or increase data acquisition interval'
#      message += EOL
#      return (status, message)

    # Calculate neutral gases total pressure in torr (only for output, not be used in calculations)
    parameters.p_neutrals_torr = pressure_conversion(parameters.p_neutrals,'Pa','torr')

    # Calculate neutral gases total number density (in m**-3)
    parameters.neutrals_density = number_density(parameters.p_neutrals, parameters.T_neutrals)

    # Calculate electrodes area (in m**2)
    parameters.area  = parameters.length * parameters.length

    # Calculate plasma volume (in m**3)
    parameters.volume = parameters.area * parameters.distance

    # Calculate starting electrons number density (in m**-3)
    parameters.start_e_density = parameters.start_ion_deg * parameters.neutrals_density

    # Calculate starting number of electrons
    parameters.N0_electrons = int(parameters.start_e_density * parameters.volume)

    # Initilize particles weight
    parameters.start_weight = START_WEIGHT

    # If the number of starting electrons is greater than the maximum allowed, lower it and increase the weight 
    if (parameters.N0_electrons > parameters.Nmax_particles): 
        parameters.start_weight = 10.0**math.ceil( math.log10( parameters.N0_electrons*1.0 / parameters.Nmax_particles ) ) 
        parameters.N0_electrons = int(parameters.N0_electrons / parameters.start_weight)
        if verbose:
            print("\nWARNING: starting weight was increased to avoid exceeding max number of electrons")

    # If the number of starting electrons is less than 1, set it to 1 and recalculate parameters
    if (parameters.N0_electrons < 1): 
        parameters.N0_electrons = 1
        if verbose:
            print("\nWARNING: starting ionization degree was increased to get at least one electron")

    # Recaculate electron density and ionization degree on the basis of the number of electrons that was obtained
    parameters.start_e_density = parameters.N0_electrons * parameters.start_weight / parameters.volume
    parameters.start_ion_deg   = parameters.start_e_density / parameters.neutrals_density

    return (status, message)


def initialize_cross_sections(neutrals, cl_options, recombination=False):

    (status, message) = (0, OK)

    # +---------------------------------------------+
    # | Load electron-neutral impact cross sections |
    # +---------------------------------------------+

    if (cl_options.verbosity > 0): print('\nReading cross-section tables for e-/neutrals scattering ...')

    # Read elastic scattering cross-sections for electron impact
    for neutral_index in range(neutrals.types):
        filename_sigma = os.path.join(XSECT_DIRECTORY_NAME, SCRIPTNAME + '_sigma_' + neutrals.names[neutral_index] + '_elastic' + EXT)
        if (cl_options.verbosity > 1): print('\n--> file \"' + filename_sigma + '\"')
        ERROR = '\nERROR in file \"' + filename_sigma+'\": '
        (status, message) = neutrals.read_xsec_electrons_elastic(filename_sigma, '\t', neutral_index, 
                                                                 plot=( cl_options.plot_xsec and (cl_options.debug_lev > 1) ) )
        if (status != 0): return (status, ERROR + message )
        if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Read ionization cross-sections for electron impact
    for neutral_index in range(neutrals.types):
        filename_sigma = os.path.join(XSECT_DIRECTORY_NAME,
                                      SCRIPTNAME + '_sigma_' + neutrals.names[neutral_index]
                                      + '_ionization' + EXT ) 
        if (cl_options.verbosity > 1): print('\n--> file \"'+filename_sigma+'\"')
        ERROR = '\nERROR in file \"'+filename_sigma+'\": '
        (status, message) = neutrals.read_xsec_electrons_ionization(filename_sigma, '\t', neutral_index, 
                                                                    check=True,
                                                                    plot=( cl_options.plot_xsec and \
                                                                           (cl_options.debug_lev > 1)
                                                                         )
                                                                   )
        if (status != 0): return (status, ERROR + message)
        if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Read excitation cross-sections for electron impact
    for neutral_index in range(neutrals.types):
        for exc_type in range(neutrals.excitation_types[neutral_index]):
            filename_sigma = os.path.join(XSECT_DIRECTORY_NAME,
                                          SCRIPTNAME+'_sigma_' + neutrals.names[neutral_index]
                                          + '_excitation_' + str(exc_type)+EXT )
            if (cl_options.verbosity > 1): print('\n--> file \"'+filename_sigma+'\"')
            ERROR = '\nERROR in file \"'+filename_sigma+'\": '
            (status, message) = neutrals.read_xsec_electrons_excitation(
                    filename_sigma, '\t', neutral_index, exc_type, check=True,
                    plot=( cl_options.plot_xsec and (cl_options.debug_lev > 1) ) )
            if (status != 0): return (status, ERROR + message )
            if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()
        
    # Read dissociation cross-sections for electron impact
    for neutral_index in range(neutrals.types):
        for diss_type in range(neutrals.dissociation_types[neutral_index]):
            filename_sigma = os.path.join(XSECT_DIRECTORY_NAME,
                                          SCRIPTNAME+'_sigma_' + neutrals.names[neutral_index]
                                          + '_dissociation_' + str(diss_type)+EXT )
            if (cl_options.verbosity > 1): print('\n--> file \"'+filename_sigma+'\"')
            ERROR = '\nERROR in file \"'+filename_sigma+'\": '
            (status, message) = neutrals.read_xsec_electrons_dissociation(
                    filename_sigma, '\t', neutral_index, diss_type, check=True,
                    plot=( cl_options.plot_xsec and (cl_options.debug_lev > 1) ) )
            if (status != 0): return (status, ERROR + message )
            if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Calculate total cross-sections, scattering frequencies and probabilities for electron impact
    (status, message) = neutrals.calculate_total_xsec_electrons()
    if (status != 0): return (status, ERROR + message )

    # Plot total cross-sections, scattering rates and probabilities for electron impact
    (status, message) = neutrals.plot_xsec_electrons( plot_single      = cl_options.plot_xsec,
                                                      plot_total       = cl_options.plot_xsec, 
                                                      plot_frequencies = cl_options.plot_xsec and (cl_options.debug_lev > 0), 
                                                      plot_relative    = cl_options.plot_xsec and (cl_options.debug_lev > 0), 
                                                      plot_boundaries  = cl_options.plot_xsec and (cl_options.debug_lev > 0) 
                                                    )
    if (status != 0):
        return (status, ERROR + message )
    if ( (cl_options.debug_lev > 0) and (not cl_options.batch_mode)): wait_input()


    # +----------------------------------------+
    # | Load ion-neutral impact cross sections |
    # +----------------------------------------+

    if (cl_options.verbosity > 0): print('\nReading cross-section tables for ion/neutrals scattering ...')
        
    # Read elastic scattering cross-sections for ions
    filename_sigma = os.path.join(XSECT_DIRECTORY_NAME, SCRIPTNAME + NAME_ION_ELASTIC + EXT)
    for ion_index in range(neutrals.types):
        for neutral_index in range(neutrals.types):
            if (cl_options.verbosity > 1): print('\n--> file \"'+filename_sigma+'\"')
            ERROR = '\nERROR in file \"'+filename_sigma+'\": '
            (status, message) = neutrals.read_xsec_ions_elastic(
                    filename_sigma, '\t', ion_index, neutral_index,
                    ( cl_options.plot_xsec and (cl_options.debug_lev > 1) ) )
            if (status != 0): return(status,  ERROR + message)
            if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Read charge exchange scattering cross-sections for ions
    filename_sigma = os.path.join(XSECT_DIRECTORY_NAME, SCRIPTNAME + NAME_ION_CHARGE_EX + EXT)
    for ion_index in range(neutrals.types):
        for neutral_index in range(neutrals.types):
            if (cl_options.verbosity > 1):
                print('\n--> file \"' + filename_sigma + '\"')
            ERROR = '\nERROR in file \"' + filename_sigma + '\": '
            (status, message) = neutrals.read_xsec_ions_charge_exchange(
                filename_sigma, '\t', ion_index, neutral_index,
                ( cl_options.plot_xsec and (cl_options.debug_lev > 1) ) )
            if (status != 0): return(status, ERROR + message)
            if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Calculate total cross-sections, scattering rates and probabilities for ions
    (status, message) = neutrals.calculate_total_xsec_ions( )
    if (status != 0):
        return(status, ERROR + message)

    # Plot total cross-sections, scattering rates and probabilities for ions
    for ion_type in range(neutrals.types):
        (status, message) = neutrals.plot_xsec_ions( ion_type         = ion_type, 
                                                     plot_total       = cl_options.plot_xsec, 
                                                     plot_frequencies = cl_options.plot_xsec and (cl_options.debug_lev > 0), 
                                                     plot_relative    = cl_options.plot_xsec and (cl_options.debug_lev > 0)
                                                   )
        if (status != 0): return(status, ERROR + message)
        if ( (cl_options.debug_lev > 0) and (not cl_options.batch_mode)): wait_input()  

                
    # +------------------------------------------------+
    # | Load electron-ion recombination cross sections |
    # +------------------------------------------------+

    # Read cross-sections for electron-ion recombination
    if recombination:
        if (cl_options.verbosity > 0): print('\nReading cross-section tables for e-/ion recombination ...')
        for neutral_index in range(neutrals.types):
            # Check if this neutral is an atom or a molecule
            if (neutrals.molecule_type[neutral_index] == 'a'):
                # For atoms, load the 3body recombination cross section
                filename_sigma = os.path.join(XSECT_DIRECTORY_NAME,
                                              SCRIPTNAME+'_sigma_'+neutrals.names[neutral_index]+'_recombination_3body'+EXT)
                if (cl_options.verbosity > 1): print('\n--> file \"'+filename_sigma+'\"')
                ERROR = '\nERROR reading file \"'+filename_sigma+'\": '
                (status, message) = neutrals.read_xsec_ele_ion_recomb(filename_sigma, '\t', neutral_index,
                                                                      rec_type='3body', 
                                                                      plot=( cl_options.plot_xsec and \
                                                                             (cl_options.debug_lev > 0) ) )
                if (status != 0): return (status,  ERROR + message)
                if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()
            else:
                # For molecules, load the dissociative recombination cross section
                for diss_type in range(neutrals.dissociation_types[neutral_index]):
                    filename_sigma = os.path.join(XSECT_DIRECTORY_NAME,
                                                  SCRIPTNAME
                                                  + '_sigma_'+neutrals.names[neutral_index]
                                                  + '_recombination_diss_'
                                                  + str(diss_type)+EXT )
                    if (cl_options.verbosity > 1):
                        print('\n--> file \"'+filename_sigma+'\"')
                    ERROR = '\nERROR reading file \"'+filename_sigma+'\": '
                    (status, message) = neutrals.read_xsec_ele_ion_recomb(filename_sigma, '\t',
                                                                          neutral_index,
                                                                          rec_type='dissociative',
                                                                          diss_type=diss_type,
                                                                          plot=\
                                                                          ( cl_options.plot_xsec and
                                                                            ( cl_options.debug_lev > 0)
                                                                          )
                                                                         )
                    if (status != 0):  return (status,  ERROR + message)
                    if ( (cl_options.debug_lev > 1) and (not cl_options.batch_mode)): wait_input()

    # Calculate total cross-sections and scattering rates for electron recombination
    # (to do also if recombination is not active)
    (status, message) = neutrals.calculate_total_xsec_ele_ion_recomb()
    if (status != 0): return(status,  ERROR + message)

    # Plot total cross-sections and scattering rates for electron recombination, but only if it is active
    if recombination:
        (status, message) = neutrals.plot_xsec_ele_ion_recomb(
            plot_total = cl_options.plot_xsec,
            plot_rates = cl_options.plot_xsec and (cl_options.debug_lev > 0) )
        if (status != 0): return (status,  ERROR + message)
        if ( (cl_options.debug_lev > 0) and (not cl_options.batch_mode)): wait_input()

    return (status, message)   
    

def initialize_ensambles(charges, neutrals, parameters, cl_options):
    """ Initialize position and velocity of electrons and ions """
    
    (status, message) = (0, OK)

    if (cl_options.verbosity > 0): print('\nInitializing electron ensamble...')

    charges.time = 0
    charges.dt   = parameters.dt

    # Initialize the electron ensamble
    charges.names[0]             = 'e-'
    charges.charge[0]            = ELECTRON_CHARGE
    charges.mass[0]              = ELECTRON_MASS
    charges.cm_ratio[0]          = charges.charge[0] / charges.mass[0]
    charges.weight[0]            = parameters.start_weight

    # Set all particles as inactive, i.e. wipe clear the ensambles before introducing new particles
    # also set all leap-frog flags to true, so that when a particle will be activated,
    # the leap-frog scheme will start correctly
    charges.active.fill(False)
    charges.restart_lf.fill(True)

    # Activate a random set of electrons
    for i in range(parameters.N0_electrons):
        # Set particle active
        charges.active[0][i] = True
        # Set starting position
        charges.x[0][i]  = random() * parameters.length
        charges.y[0][i]  = random() * parameters.length
        charges.z[0][i]  = random() * parameters.distance
        # Set starting velocity components
        mean_v_el = numpy.sqrt(K_BOLTZMANN * neutrals.temperature / charges.mass[0])
        (vx, vy, vz, v) = random_maxwell_velocity(mean_v_el, modulus=True)
        charges.vx[0][i] = vx
        charges.vy[0][i] = vy
        charges.vz[0][i] = vz
        charges.v[0][i]  = v
#        charges.v[0][i]  = math.sqrt(  charges.vx[0][i] * charges.vx[0][i]
#                                     + charges.vy[0][i] * charges.vy[0][i]
#                                     + charges.vz[0][i] * charges.vz[0][i] ) 
        if (charges.v[0][i] != 0): 
            charges.theta[0][i] = math.acos( charges.vz[0][i] / charges.v[0][i] )

    if (cl_options.verbosity > 0): print('\nInitializing ion ensamble...')                        

    # Activate a random set of ions
    for i in range(1, neutrals.types+1):
        charges.names[i]             = neutrals.names[i-1]+'+'
        charges.charge[i]            = - ELECTRON_CHARGE
        charges.mass[i]              = neutrals.mass[i-1] * ATOMIC_UNIT_MASS
        charges.cm_ratio[i]          = charges.charge[i] / charges.mass[i]
        charges.weight[i]            = parameters.start_weight

        # Starting density of ions of type i is proportional to the density of neutral species i 
        parameters.N0_ions = int(parameters.N0_electrons * neutrals.number_density[i-1] / neutrals.number_density.sum())
        for j in range(parameters.N0_ions):
            # Set particle active
            charges.active[i][j] = True
            # Set starting position
            charges.x[i][j]  = random() * parameters.length
            charges.y[i][j]  = random() * parameters.length
            charges.z[i][j]  = random() * parameters.distance
            # Set starting velocity components
            (vx, vy, vz, v) = random_maxwell_velocity(neutrals.mean_v[i-1], modulus=True)
            charges.vx[i][j] = vx
            charges.vy[i][j] = vy
            charges.vz[i][j] = vz
            charges.v[i][j]  = v           
#            charges.v[i][j]  = math.sqrt(  charges.vx[i][j] * charges.vx[i][j]
#                                         + charges.vy[i][j] * charges.vy[i][j]
#                                         + charges.vz[i][j] * charges.vz[i][j] )
            # Set angle between velocity vector and z-axis
            if (charges.v[i][j] != 0): 
                charges.theta[i][j] = math.acos( charges.vz[i][j] / charges.v[i][j] )

    return (status, message)


def initialize_filenames(neutrals, parameters, string=None):
    
    name = ''
    if (string is None):
        for i in range(neutrals.types):
            name += neutrals.names[i]
            if ( (neutrals.types > 1) and (i < neutrals.types-1) ): name += '+'
        name += '_'
        (value, unit) = unit_manager.change_unit(parameters.p_neutrals, 'Pa')
        name += str(value) + unit.strip() + '_'
        (value, unit) = unit_manager.change_unit(parameters.distance, 'm')
        name += str(value) + unit.strip() + '_'
        (value, unit) = unit_manager.change_unit(parameters.V_bias, 'V')
        name += str(value) + unit.strip() + '_'        
        (value, unit) = unit_manager.change_unit(parameters.frequency, 'Hz')
        name += str(value) + unit.strip()
    else:
        name = str(string)
    parameters.basename = name

        
def generate_save_dir_name(parameters, abs_path=True):
    
    # Generate save directory path        
    parameters.path = os.path.abspath('.')

    if abs_path:
        parameters.save_directory = os.path.join(parameters.path, SAVE_DIRECTORY_NAME, parameters.basename)
    else:
        parameters.save_directory = parameters.path

        
def generate_save_file_names(parameters):
    
    # Add save directory path to filenames
    parameters.filename_stat_ele    = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_STAT_ELE)
    parameters.filename_stat_ion    = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_STAT_ION)    
    parameters.filename_stat_neu    = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_STAT_NEU)
    parameters.filename_distrib_ele = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_DISTRIB_ELE)
    parameters.filename_distrib_ion = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_DISTRIB_ION)
    parameters.filename_I           = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_I)
    parameters.filename_V           = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_V)
    parameters.filename_epos_z      = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_Z_ELECTRONS)
    parameters.filename_ipos_z      = os.path.join(parameters.save_directory, SCRIPTNAME + NAME_Z_IONS)
    parameters.filename_config      = os.path.join(parameters.save_directory, FILENAME_CONFIG)
    parameters.filename_neutrals    = os.path.join(parameters.save_directory, FILENAME_NEUTRALS)
        

def create_save_dir(parameters):
        
    generate_save_dir_name(parameters)

    if not os.path.exists(SAVE_DIRECTORY_NAME): os.mkdir(SAVE_DIRECTORY_NAME)

    # If the directory already exists, change the name
    if os.path.exists(parameters.save_directory):
        i = 0
        while True:
            i += 1
            string = parameters.save_directory + ' (' + str(i) + ')'
            if not os.path.exists(string): break
        parameters.save_directory = string
    os.mkdir(parameters.save_directory)

    generate_save_file_names(parameters)

    # Copy configuration and neutrals files to the save directory
    shutil.copy(FILENAME_CONFIG,   parameters.filename_config)
    shutil.copy(FILENAME_NEUTRALS, parameters.filename_neutrals)
 
