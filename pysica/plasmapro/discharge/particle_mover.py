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

    This module contains functions to simulate the motion of particles in plasma discharges, 
    it relies on the fortran module providing an interface with it.

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import sys
import math
import numpy
from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *
from pysica.io.io_screen import wait_input
from .fortran.fmodule import f_main

MIN_DEBUG_LEVEL = 2


def calculate_dt(charges, neutrals, ccp, parameters, debug_level):
        """ Calculate timestep  """
        
        # Let's calculate a convenient dt
        # we need that dt <= maxcollfreq*taumin
        # where taumin is minimum time between collisions
        # taumin = (n*s*vmax)**-1
        # vmax is the maximum velocity (i.e. without collisions)
        # that electrons could reach in the time dt_output
        a       = -charges.cm_ratio[0] * math.fabs(ccp.E_peak)  # Particles acceleration without collisions
        v0      = charges.v_average(0)                          # mean electron velocity [m/s]
        e0      = charges.e_average(0)                          # mean electron kinetic energy [eV]
        vmax    = v0 + a * parameters.dt_output                 # maximum expected velocity in time dt_output [m/s]
        emax    = -0.5/charges.cm_ratio[0]*vmax**2              # maximum expected energy in time dt_output [eV]

        if (debug_level >= MIN_DEBUG_LEVEL):
                print('-> Calculating variable dt ')
                print('a=', a, end=' ')
                print('v0=', v0, end=' ')
                print('e0=', e0, end=' ')
                print('vmax=', vmax, end=' ')
                print('emax=', emax)

        # Calculate the best timestep
        if ( (vmax == 0) or (neutrals.number_density.sum() == 0) ):
                dt1 = dt_output
        else:
                # Find the nearest value of cross section for the actual average velocity
                for j in range(neutrals.n_sigma):
                        if (v0 <= neutrals.sigma_velocities[j]): break 
                sigma   =  neutrals.sigma_total_global[j]       

                # Calculate the minimum expected collision period [s]
                taumin  = (neutrals.number_density.sum() * sigma * vmax)**-1            
                dt1     = parameters.maxcollfreq*taumin                         # dt [s]

        # Get only the order of magnitude of dt
        dt = 10**round(math.log10(dt1))
                
        # Check that timestep is not greater than time between data acquisitions
        if (dt > parameters.dt_output):         dt = parameters.dt_output                       

        # Check that the timestep is much lower than the period of electric field       
        if (ccp.frequency != 0):
                dtmax = 1.0 / (ccp.frequency*MINTRATIO)
                while (dt > dtmax):
                        dt = dt / 10.0

        if (debug_level >= MIN_DEBUG_LEVEL):
                print('E_0= ', e0,' eV', end=' ') #+unit_manager.print_unit(e0,   'eV', 3),
                print('\te_max= ', emax, ' eV', end=' ') #+unit_manager.print_unit(emax, 'eV', 3),  
                print('\t\tdt= ', dt, ' s') # +unit_manager.print_unit(dt,    's', 3),
                print() 

        if ( (dt <= 0) or numpy.isnan(dt) ):
                print("ERROR calculating dt")
                print("dt1     = " + str(dt1))
                print("dt      = " + str(dt))
                print("a       = " + str(a))
                print("v0      = " + str(v0))
                print("e0      = " + str(e0))
                print("vmax    = " + str(vmax))
                print("emax    = " + str(emax))
                wait_input()
                sys.exit()

        return dt




# +-----------------------------------+
# | Move particles in a CCP discharge |
# +-----------------------------------+

def move_particles(charges, neutrals, ccp, parameters, options):
        """Simulate the motion of charged particles (electrons and ions) in a cold plasma.

               Parameters
               ----------

               charges:     instance of the class MovingParticles, a collection of charged particles
               neutrals:    instance of the class TargetParticles, a collection of neutral particles
               ccp:         instance of the class CcpProperties, reactor properties
               duration:    required duration of the particle motion (simulation time)
               dt-vat:      if True, timestep mus be recalculated
               maxcollfreq: maximum allowed vlaue collision frequency

               Initialized data attributes
               ---------------------------
               
               Returns
               -------
        """

        if (options.debug_lev >= MIN_DEBUG_LEVEL): print('\n***** START OF PARTICLE MOVER *****\n')

        if ( parameters.dt_var and (charges.n_active(0) > 0) ):
                charges.dt = calculate_dt(charges, neutrals, ccp, parameters, options.debug_lev)
        
        # Reset collision counters
        neutrals.collisions_null.fill(0.0)
        neutrals.collisions_elastic.fill(0.0)
        neutrals.collisions_ionization.fill(0.0)
        neutrals.collisions_dissociation.fill(0.0)
        neutrals.collisions_recombination.fill(0.0)
        neutrals.dissociation_rate.fill(0.0)
        neutrals.dissociation_rate_constant.fill(0.0)

        # Initialize the variable used to store the average electric current of the previus cycle
        # the variable is declared as a rank-zero numpy array in order to get working f2py automatic conversion 
        # between fortran and python
        ccp.average_current.fill(0)

        # This is needed since f2py seems not able to convert correctly logical arrays between Fortran and Python,
        # while it works well if using integer type arrays
        active  = charges.active.astype('i')
        restart = charges.restart_lf.astype('i')

        if (options.debug_lev >= MIN_DEBUG_LEVEL): print('-> Calling Fortran function')

        f_main.simccp(charges.x, charges.y, charges.z, 
                          charges.vx, charges.vy, charges.vz, charges.v, 
                          active, restart,
                          charges.cm_ratio, 
                          charges.weight,
                          charges.rescale_factor,
                          neutrals.sigma_velocities, 
                          neutrals.frequency_total_global_max, 
                          neutrals.sigma_velocities_ions, 
                          neutrals.frequency_total_global_max_ions, 
                          neutrals.probability_limits, 
                          neutrals.probability_limits_ions,
                          neutrals.n_limits, neutrals.n_limits_ions, 
                          neutrals.ratecoeff_recomb_diss, 
                          neutrals.isactive_recomb, 
                          neutrals.min_scattered,
                          neutrals.v_ratio, 
                          neutrals.v_ratio_ions,
                          neutrals.ionization_energy, neutrals.dissociation_energy, 
                          neutrals.dissociation_types, 
                          neutrals.mean_v,
                          neutrals.secondary_emission,
                          ccp.distance, ccp.length, ccp.V_peak, ccp.pulsation, ccp.phase, ccp.lateral_loss,
                          ccp.charge_density, ccp.potential,
                          ccp.average_current,
                          charges.dt, parameters.dt_output, 
                          neutrals.collisions_null, neutrals.collisions_elastic, neutrals.collisions_ionization, 
                          neutrals.collisions_dissociation, neutrals.collisions_recombination,
                          options.debug_lev_for)

        # This is needed since f2py seems not able to convert correctly logical arrays between Fortran and Python (see above)
        charges.active     = active.astype('bool')
        charges.restart_lf = restart.astype('bool')

        # Update simulation time and phase of electric bias
        charges.time = charges.time + parameters.dt_output
        ccp.phase    = ccp.phase + ccp.pulsation * parameters.dt_output
        if (ccp.pulsation > 0):
                ccp.V  = ccp.V_peak * math.sin(ccp.phase)
        else:
                ccp.V  = ccp.V_peak

        # Calculate the total number of electron collisions, of any type
        neutrals.collisions_total_electron = neutrals.collisions_elastic.sum()      + neutrals.collisions_ionization.sum() + \
                                             neutrals.collisions_dissociation.sum() + neutrals.collisions_recombination.sum()

        # Calculate average time between electron collisions (rough estimation, since charges.n_active(0) changes during the cycle)
        if (neutrals.collisions_total_electron > 0):
                charges.tau_mis = parameters.dt_output / neutrals.collisions_total_electron * charges.n_active(0) * charges.weight[0]
        else:
                charges.tau_mis = 0

        # Calculate average electron collision probability (rough estimation, since charges.n_active(0) changes during the cycle)
        N_iterations = int(parameters.dt_output/charges.dt)
        if (charges.n_active(0) > 0):
                charges.p_coll =  100.0 * neutrals.collisions_total_electron / ( N_iterations * charges.n_active(0) * charges.weight[0] )
        else:
                charges.p_coll = 0
        
        # Calculate charges number densities
        for i in range(charges.types):
                charges.number_density[i] = charges.n_active(i) * charges.weight[i] / ccp.volume

        # Calculate Debye length and plasma frequency
        for i in range(charges.types):
                if (charges.number_density[i] > 0):
                        charges.debye_length[i]     = numpy.sqrt( - EPSILONZERO * ELECTRON_CHARGE * charges.kT(i) / \
                                                                  ( charges.charge[i]**2 * charges.number_density[i] ) \
                                                                )
                else:
                        charges.debye_length[i] = 0                        
                charges.plasma_frequency[i] = numpy.sqrt( charges.charge[i]**2 * charges.number_density[i] / \
                                                          ( EPSILONZERO * charges.mass[i] ) \
                                                        )

                if (options.debug_lev >= MIN_DEBUG_LEVEL): print('-> e-/ion, kTe, Dl, f =', i, \
                                                  charges.kT(i), charges.debye_length[i], charges.plasma_frequency[i])
                
        # Calculate neutrals dissociation rates and rate constants
        neutrals.dissociation_rate = neutrals.collisions_dissociation / ccp.volume / parameters.dt_output
        for i in range(neutrals.types):
                if ( (neutrals.number_density[i] > 0) and (charges.number_density[0] > 0) ):
                        neutrals.dissociation_rate_constant[i] = neutrals.dissociation_rate[i] / \
                                                                 (neutrals.number_density[i]* charges.number_density[0])
                else:
                        neutrals.dissociation_rate_constant[i] = 0
                if (options.debug_lev >= MIN_DEBUG_LEVEL):                      
                        print('-> neutral, n0, ne, R_diss, k_diss =', i, neutrals.number_density[i], charges.number_density[0],\
                                                                           neutrals.number_density[i] * charges.number_density[0], \
                                                                           neutrals.dissociation_rate[i], \
                                                                           neutrals.dissociation_rate_constant[i])

        # Calculate angles between electron velocity and E-field, avoiding division by zero for particles having null velocity
        for i in range(charges.types):
                for k in range(charges.n):
                        if charges.active[i][k]:
                                # if an active particle has zero speed, avoid division, and set angle to zero
                                if (charges.v[i][k] == 0.0): 
                                        charges.theta[i][k] = 0.0
                                else:
                                        charges.theta[i][k] = numpy.arccos(charges.vz[i][k] / charges.v[i][k])

        if (options.debug_lev >= MIN_DEBUG_LEVEL): print('\n***** END OF PARTICLE MOVER *****\n')
        
        if (options.debug_lev >= MIN_DEBUG_LEVEL+1): wait_input()
        
