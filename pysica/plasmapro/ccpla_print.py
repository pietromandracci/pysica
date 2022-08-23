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

""" Simulation of a capacitively coupled plasma discharge. """

import math
from numpy import isnan, isinf

from pysica.managers import unit_manager
from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *

# +-----------+
# | Functions |
# + ----------+

def print_simulation_information(neutrals, ccp, parameters, options,
                                 print_physical=True, print_simulation=True,
                                 print_output=True, print_filenames=True):
    string = ''
    # Print information about physical parameters
    if print_physical:
        string += 'PHYSICAL PARAMETERS\n'
        string +=('Total pressure of the gas mixture:   ' + unit_manager.print_unit(parameters.p_neutrals,'Pa',4)
                                                          + ' ('
                                                          + unit_manager.print_unit(parameters.p_neutrals_torr,'Torr',4)
                                                          + ')' + '\n')
        string +=('Temperature of the gas mixture:      ' + unit_manager.print_unit(parameters.T_neutrals, 'K',4)
                                                          + ' ('
                                                          + unit_manager.print_unit(parameters.T_neutrals-ZERO_CELSIUS,'C',4)
                                                          + ')' + '\n')
        string += 'Total number density of molecules:   '  + unit_manager.print_exp(parameters.neutrals_density, 2) + ' m**-3' + '\n'
        string += 'Distance between electrodes:         '  + unit_manager.print_unit(parameters.distance,'m') + '\n'
        string += 'Length of electrodes:                '  + unit_manager.print_unit(parameters.length,'m') + '\n'
        string += 'Area of electrodes                   '  + unit_manager.print_exp(ccp.area, 2) + ' m**2' + '\n'
        string += 'Plasma volume:                       '  + unit_manager.print_exp(ccp.volume, 2) + ' m**3' + '\n'
        string += 'Starting electrons number density:   '  + unit_manager.print_exp(parameters.start_e_density, 2) + ' m**-3' + '\n'
        string += 'Starting ionization degree:          '  + unit_manager.print_exp(parameters.start_ion_deg, 2) + '\n'
        string += 'Electric bias peak value:            '  + unit_manager.print_unit(parameters.V_bias,'V', 4)  + '\n'
        string += 'Mean electric field intensity:       '  + unit_manager.print_unit(ccp.E_peak,'V/m', 4) + '\n'
        string += 'Electric bias frequency:             ' 
        if (ccp.frequency == 0.0):
            string += 'STATIC' + '\n'
        else:
            string += unit_manager.print_unit(ccp.frequency,'Hz', 4) + '\n'
            string += 'Electric bias period:                ' + unit_manager.print_unit(ccp.period,'s', 4) + '\n'
            string += 'Electric bias phase at t=0:          ' + str(parameters.phase * 180.0 / math.pi) + ' deg' + '\n'

    # Print information about simulation parameters
    if print_simulation:
        if not options.gui_mode: string += '\n'
        string += 'SIMULATION PARAMETERS\n'
        string += 'Maximum  number of particles:        ' + str(parameters.Nmax_particles) + '\n'
        string += 'Starting number of electrons:        ' + str(parameters.N0_electrons) + '\n'
        string += 'Starting computational weight:       ' + unit_manager.print_exp(parameters.start_weight,3) + '\n'
        string += 'Rescaling factor                     ' + str(parameters.rescale_factor) + '\n'
        string += 'Number of cells in PIC scheme:       ' + str(parameters.N_cells) + '\n'
        string += 'Cell dimension:                      ' + unit_manager.print_unit(ccp.delta_grid,'m', 4) + '\n'
        string += 'Collisions required for fast method: ' + str(neutrals.min_scattered) + '\n'
        string += 'Simulation timestep:                 '
        if parameters.dt_var:
            string += 'variable' + '\n'
            string += 'Max allowed collision frequency:     ' + str(parameters.maxcollfreq*100) + ' %' + '\n'
        else:
            string += unit_manager.print_unit(parameters.dt,'s') + '\n'
        string += 'Data output interval:                ' + unit_manager.print_unit(parameters.dt_output,'s', 4) + '\n'
        string += 'Required simulation time:            '
        if (parameters.sim_duration > 0):
            string += unit_manager.print_unit(parameters.sim_duration,'s', 4)  + '\n'
        else:
            string += 'unlimited\n'
        string += 'Electron and ion lateral loss:       '
        if ccp.lateral_loss:
            string += 'ON' + '\n'
        else:
            string += 'OFF' + '\n'
        string += 'Electron-ion recombination:          ' 
        if neutrals.isactive_recomb:
            string += 'ON' + '\n'
        else:
            string += 'OFF' + '\n'
        string += 'Electron impact cross section        ' + '\n'
        string += '- Number of cross section values:    ' + str(parameters.N_sigma) + '\n'
        string += '- Minimum cross section energy:      ' + unit_manager.print_unit(parameters.e_min_sigma, 'eV', 4) + '\n'
        string += '- Maximum cross section energy:      ' + unit_manager.print_unit(parameters.e_max_sigma, 'eV', 4) + '\n'
        string += 'Ion impact cross section             ' + '\n'
        string += '- Number of cross section values:    ' + str(parameters.N_sigma_ions) + '\n'
        string += '- Minimum cross section energy:      ' + unit_manager.print_unit(parameters.e_min_sigma_ions, 'eV', 4) + '\n'
        string += '- Maximum cross section energy:      ' + unit_manager.print_unit(parameters.e_max_sigma_ions, 'eV', 4) + '\n'

    # Print output information
    if print_output:
        if not options.gui_mode: string += '\n'
        string += 'DATA OUTPUT PARAMETERS\n'
        string += 'Verbosity level [0..3]:              ' + str(options.verbosity) + '\n'
        string += 'Python debug level [0..2]:           ' + str(options.debug_lev) + '\n'
        string += 'Fortran debug level [0..2]:          ' + str(options.debug_lev_for) + '\n'
        string += 'Plot cross sections:                 ' 
        if options.plot_xsec: string += 'YES' + '\n'
        else:                 string += 'NO' + '\n'
        string += 'Save data to file:                   '
        if (parameters.save_delay == 0):
            string += 'never' + '\n'
        elif (parameters.save_delay == 1):
            string += 'every cycle' + '\n'
        else:
            string += 'every ' + str(parameters.save_delay) + ' cycles' + '\n'
        if (print_filenames and (parameters.save_delay > 0) ):
            string += '\n'
            string += 'Filename electron mean data:         ' + '\"' + parameters.filename_stat_ele    + EXT +'\"\n'
            string += 'Filename ion mean data:              ' + '\"' + parameters.filename_stat_neu    + EXT + '\"\n'
            string += 'Filename eedf:                       ' + '\"' + parameters.filename_distrib_ele + EXT + '\"\n'
            string += 'Filename iedf:                       ' + '\"' + parameters.filename_distrib_ion + EXT + '\"\n'
            string += 'Filename z positions of electrons:   ' + '\"' + parameters.filename_epos_z      + EXT + '\"\n'
            string += 'Filename z positions of ions         ' + '\"' + parameters.filename_ipos_z      + EXT + '\"\n'            
            string += 'Filename electric current:           ' + '\"' + parameters.filename_I           + EXT + '\"\n'
            string += 'Filename electric potential:         ' + '\"' + parameters.filename_V           + EXT + '\"\n'
            
        if options.gui_mode:
            string += 'Max number of points in plots:       ' + str(parameters.n_max_points) + '\n'
            string += 'Decimation factor:                   ' + str(parameters.decimation_factor) + '\n'
                        
    return string


def print_gas_information(neutrals):
    string = ''
    for i in range(neutrals.types):
        string += '\n' 
        string += 'Name:                                     ' + str(neutrals.names[i]) + '\n'
        if   (neutrals.molecule_type[i]=='a'):
            string += 'Type:                                     ' + 'atom' + '\n'
        elif (neutrals.molecule_type[i]=='m'):
            string += 'Type:                                     ' + 'non polar molecule' + '\n'
        elif (neutrals.molecule_type[i]=='p'):
            string += 'Type:                                     ' + 'polar molecule' + '\n'
        string +=('Flow rate percentage:                     '
                  + str( unit_manager.fix_digits(neutrals.gas_flow[i]/neutrals.gas_flow.sum()*100.0, 3) )
                  + ' %' + '\n')
        string +=('Partial pressure [Pa]:                    '
                  + str( unit_manager.fix_digits(neutrals.partial_pressure[i], 3) ) + '\n')
        string +=('Number density [m**-3]:                   ' + str( unit_manager.fix_digits(neutrals.number_density[i], 4) )
                                                               + '\n')
        string += 'Mass [u]:                                 ' + str( neutrals.mass[i] ) + '\n'
        string += 'Mean molecule speed [m s**-1]             ' + str( unit_manager.fix_digits(neutrals.mean_v[i], 4) ) + '\n'   
        string +=('Electron energy loss (elastic):           ' + str( unit_manager.fix_digits(neutrals.energy_loss[i], 6) )
                                                               + '\n')
        string += 'Secondary emission coefficient            ' + str( neutrals.secondary_emission[i] ) + '\n'
        string += 'Ionization energy [eV]:                   ' + str( neutrals.ionization_energy[i] ) + '\n'
        string += 'Number of excitation processes:           ' + str( neutrals.excitation_types[i] ) + '\n'
        string += 'Excitation energies [eV]                  '
        for j in range(neutrals.excitation_types[i]):
            string += str(neutrals.excitation_energy[i][j])
        string += '\n'

        if (neutrals.molecule_type[i] != 'a'):
            string += 'Number of dissociation processes:         ' + str( neutrals.dissociation_types[i] ) + '\n'
            string += 'Dissociation energies [eV]                ' 
            for j in range(neutrals.dissociation_types[i]):
                string += str(neutrals.dissociation_energy[i][j])
            string += '\n'
            
        #string += '\nTotal number density relative error = ' + \
        #           str( (neutrals.number_density.sum() - neutrals_density) / neutrals_density ) + '\n'

    return string


def print_runtime_info(charges, neutrals, ccp, parameters, options, time_before, electric_bias_before, n_active_el_before):

    # N. of char allocated to place gas names
    gas_length = 10 
    
    text = ''
#   text += 'Elapsed time before iteration            = ' + unit_manager.print_unit(time_before, 's', 4) + '\n'
    text += 'Iteration time                           = ' + unit_manager.print_unit(parameters.dt_output, 's', 3) + '\n'
    text += 'Timestep                                 = ' + unit_manager.print_unit(charges.dt, 's', 3) + '\n'
    text += 'Computational electrons before cycle     = ' + str(n_active_el_before)
    text += ' (max ' + str((parameters.Nmax_particles))+')' + '\n'
    text += 'Electric bias intensity before cycle     = ' + unit_manager.print_unit(electric_bias_before,'V',3)
    text += ' (max = ' + unit_manager.print_unit(ccp.V_peak, 'V', 3) + ')' + '\n'
    text += 'Reduced electric field before cycle      = ' + unit_manager.print_unit(electric_bias_before
                                                                                   /ccp.distance
                                                                                    /parameters.neutrals_density*1.0e21,'Td',3) 
    text += ' (max = ' + unit_manager.print_unit(ccp.V_peak/ccp.distance/parameters.neutrals_density*1.0e21,'Td',3) + ')\n'
    text += '...' + '\n'
    text += 'Elapsed time after cycle                 = ' + unit_manager.print_unit(charges.time, 's', 4)
    if (parameters.sim_duration > 0):
        text +=  ' of ' + unit_manager.print_unit(parameters.sim_duration, 's', 4)                
        text += ' (' + str('%3.1f' % (100.0*(charges.time)/parameters.sim_duration)) + '%) '
    else:
        text += ' (unlimited)'
    text += '\n'

    text += 'Average time between collisions          = ' + unit_manager.print_unit(charges.tau_mis,'s',3) + '\n'
    text += 'Average current density                  = ' + unit_manager.print_unit(ccp.average_current[0]/ccp.area,'A/m**2',3)\
                                                      + '\n'
    charge_density = 0.0
    for i in range(charges.types):
        if (i==0): charge_density = charge_density + charges.number_density[i] 
        else:      charge_density = charge_density - charges.number_density[i]
    charge_density = charge_density * ELECTRON_CHARGE #ELECTRON_CHARGE is negative
    text += 'Net charge density                       = ' + unit_manager.print_unit(charge_density,'C*m**-3',3) + '\n'
    text += 'Real electron collisions                 = ' + unit_manager.print_exp(neutrals.collisions_total_electron, 3) + \
                                                            str(' (%2.2f' % charges.p_coll) + '%) ' + '\n'
    text += 'Null electron collisions                 = ' + unit_manager.print_exp(neutrals.collisions_null[0], 3) + '\n'
    if True: #(options.verbosity > 2):
        text += '\n'
        text += 'Number of computational particles (max '+str(parameters.Nmax_particles) + \
                ') and physical parameters for electrons and ions' + '\n'
        for i in range(charges.types):
            string = ''
            string += charges.names[i].ljust(gas_length) + ':'
            string += '  N = '        + str(charges.n_active(i)).rjust(8) + ';'
            #string += '  w = '        + unit_manager.print_exp(charges.weight[i],2, align=True).rjust(8) + ';'
            string += '  w = '        + unit_manager.print_exp(charges.weight[i],2, align=True).rjust(5) + ';' 
            string += '  n = '        + unit_manager.print_exp(charges.number_density[i],
                                                               2, align=True).rjust(8) + ' m**-3'  + ';'
            string += '  alpha = '    + unit_manager.print_exp(charges.number_density[i]/parameters.neutrals_density,
                                                               2, align=True).rjust(8) + ';'
            string += '  <E> = '      + unit_manager.print_unit(charges.e_average(i), 'eV',
                                                                3, align=True).ljust(9) + ';'
#            string += '  E/N = '      + unit_manager.print_exp(ccp.V/ccp.distance/parameters.neutrals_density,
#                                                               2, align=True).rjust(8) + ' Td'  + ';'       
            string += '  lD = '       + unit_manager.print_unit(charges.debye_length[i],     'm',
                                                                3, align=True).ljust(9) + ';' 
            string += '  ND = '       + unit_manager.print_exp( 4.0/3.0**math.pi*charges.debye_length[i]**3*\
                                                                            charges.number_density[i], 2,
                                                                            align=True ).rjust(8) + ';'            
            if (charges.debye_length[i] >0):
                    string += '  f = '        + unit_manager.print_unit(charges.v_average(i)/charges.debye_length[i],
                                                                        'Hz', 3, align=True).rjust(9)
            else:
                    string += '  f = '        + 'N/A'.rjust(9)
            if (charges.v_average(i) > 0):
                    string += ' ('            + unit_manager.print_unit(charges.plasma_frequency[i], 'Hz',
                                                                        3, align=True) + ');'
            else:
                    string += '  t = '        + 'N/A'
            text += string + '\n'     

        # Print table of collision rates
            
        text += '\n'
        text += 'Electron collisions' + '\n'
        string = ''
        length = 14
        string += ( 'Gas'.ljust(gas_length)
                    + '| Total'.ljust(length+1)
                    + '| Elastic'.ljust(length+1)
                    + '| Ionization'.ljust(length+1)
                    + '| Excitation'.ljust(length+1)
                    + '| Dissociation'.ljust(length+1)
                    + '| Recombination'.ljust(length+1)                    
                   )

#        for j in range(MAX_DISSOCIATION_TYPES):
#                string += '| Type' + str(j).rjust(2) + '  '
        text += string + '|' + '\n'   
#        string_line = ''.ljust(gas_length+1 + (length+1)*6 + 10*MAX_DISSOCIATION_TYPES, '-')
        string_line = ''.ljust(gas_length+1 + (length+1)*6, '-')
        text += string_line + '\n'
        text += ( 'TOT'.ljust(gas_length) + '|'
                 + unit_manager.print_exp(neutrals.collisions_total_electron,2).rjust(length)      + '|'
                 + unit_manager.print_exp(neutrals.collisions_elastic.sum(),2).rjust(length)       + '|'
                 + unit_manager.print_exp(neutrals.collisions_ionization.sum(),2).rjust(length)    + '|'
                 + unit_manager.print_exp(neutrals.collisions_excitation.sum(),2).rjust(length)    + '|'
                 + unit_manager.print_exp(neutrals.collisions_dissociation.sum(),2).rjust(length)  + '|'
                 + unit_manager.print_exp(neutrals.collisions_recombination.sum(),2).rjust(length) + '|'                  
                 + '\n' )
        text += string_line + '\n'

        for i in range(neutrals.types):
            string = ''
            string += str(neutrals.names[i]).ljust(gas_length) + '|'
            string += unit_manager.print_exp( ( neutrals.collisions_elastic[i]
                                                + neutrals.collisions_ionization[i]
                                                + neutrals.collisions_excitation[i].sum()
                                                + neutrals.collisions_dissociation[i].sum()
                                                + neutrals.collisions_recombination[i].sum()) ,
                                                2 ).rjust(length) + '|'
            string += unit_manager.print_exp( neutrals.collisions_elastic[i],              2 ).rjust(length) + '|'
            string += unit_manager.print_exp( neutrals.collisions_ionization[i],           2 ).rjust(length) + '|'
            string += unit_manager.print_exp( neutrals.collisions_excitation[i].sum(),     2 ).rjust(length) + '|'
            if (neutrals.molecule_type[i] == 'a'):
                string += '-'.rjust(length) + '|'
            else:
                string += unit_manager.print_exp(neutrals.collisions_dissociation[i].sum(), 2).rjust(length)  + '|'
#            for k in range(MAX_DISSOCIATION_TYPES):
#                if ( (neutrals.molecule_type[i] == 'a') or (k >= neutrals.dissociation_types[i]) ):
#                    string += '-'.rjust(9) + '|'
#                else:
#                    string += unit_manager.print_exp(neutrals.collisions_dissociation[i][k],2).rjust(9) + '|'
            string += unit_manager.print_exp( neutrals.collisions_recombination[i].sum(),  2 ).rjust(length) + '|'                    
            text += string + '\n'

        # Print table of dissociation rates and rate constants
        
        text += '\n'
        text += 'Dissociation rate parameters' + '\n'
        string = ''
        length = 20
        string += 'Gas'.ljust(gas_length)
        for j in range(MAX_DISSOCIATION_TYPES):
            string += (  ('| Diss. rate '+str(j)).ljust(length+1)
                       + ('| Diss. rate const. '+str(j)).ljust(length+1) )
        text += string + '|' + '\n'   
        string_line = ''.ljust(gas_length+1 + (length+1)*2*MAX_DISSOCIATION_TYPES , '-')
        text += string_line + '\n'        
        string = 'Unit'.ljust(gas_length)
        for j in range(MAX_DISSOCIATION_TYPES):
            string += (   '| m**-3 * s**-1'.ljust(length+1)
                        + '| m**3  * s**-1'.ljust(length+1)  )
        text += string + '|' + '\n'
        text += string_line + '\n'

        for i in range(neutrals.types):
            string = ''
            # If this gas is atomic, no dissociation is possible, skip it
            if (neutrals.molecule_type[i] != 'a'):
                string += str(neutrals.names[i]).ljust(gas_length) + '|'
                for j in range(MAX_DISSOCIATION_TYPES):
                    if (j >= neutrals.dissociation_types[i]):
                        string += '-'.rjust(length) + '|'
                        string += '-'.rjust(length) + '|'                        
                    else:
                        string += unit_manager.print_exp(neutrals.dissociation_rate[i][j],2).rjust(length) + '|'
                        k = neutrals.dissociation_rate_constant[i][j]
                        if ( isnan(k) or isinf(k) ): 
                            string += 'N/A'
                        else:
                            string += unit_manager.print_exp(k,2).rjust(length) + '|'
                #string += '|'
                text += string + '\n'                    

    return text
