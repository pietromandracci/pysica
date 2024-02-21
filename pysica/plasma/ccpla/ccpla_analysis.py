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

""" Read and show data saved by the ccpla script during a simulation run """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

# Mudules from the standard Python library
import math

# Modules provided by the Python community
import numpy
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.inset_locator as inset_locator
import Gnuplot

# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.plasma.ccpla.ccpla_defaults import *

# Modules provided by plasma.ccpla package
from pysica.managers import data_manager
from pysica.managers.unit_manager import print_unit, print_exp
from pysica.analysis import univariate
from pysica.functions.pdf import pdf_maxwell_energy
from pysica.plasma.ccpla.discharge.reactors import CcpProperties
from pysica.plasma.ccpla.discharge.target_particles import TargetParticles
from pysica.plasma.ccpla.ccpla_init import *


class CcplaSavedData:
    """Class of data saved by ccpla program"""

    def __init__(self, verbose=False, read_potential=False, read_edf=False, read_position=False,
                 size_default=26, size_title=30, size_axes=26, size_legend=26):
        """Initialize the collection of data to be analyzed"""

        self.tick_length_major = 10
        self.tick_length_minor = 5

        # Set the font dimensions
        plt.rcParams["mathtext.default"] = 'regular'
        plt.rcParams["font.size"]        = 20
        plt.rc('font',   size      = size_default) #controls default text size
        plt.rc('axes',   titlesize = size_title)   #fontsize of the title
        plt.rc('axes',   labelsize = size_axes)    #fontsize of the x and y labels
        plt.rc('xtick',  labelsize = size_axes)    #fontsize of the x tick labels 
        plt.rc('ytick',  labelsize = size_axes)    #fontsize of the y tick labels
        plt.rc('legend', fontsize  = size_legend)  #fontsize of the legend

        # Reading EEDF and IEDF as well as potential distribution can be suppressed to save memory and execution time
        self.read_V        = read_potential
        self.read_edf      = read_edf
        self.read_z        = read_position

        # Define some constants
        self.DEBYE_CONST_A = 1.51E13                                # (J m)**(-1/2)
        self.DEBYE_CONST_B = 4/3 * numpy.pi * self.DEBYE_CONST_A**3 # (J m)**(-3/2)
        self.DEBYE_CONST_C = 56.41                                  # m**(3/2) s**-1         

        # Read data from configuration files
        self.parameters   = ConfigurationOptions()
        (status, message) = initialize_parameters(self.parameters, verbose=verbose, saveonly=False,
                                                  filename_config=FILENAME_CONFIG)
        if (status != 0):
            self.error = (status, message)
            return

        if verbose: print('\nLoading reactor properties\n')
        self.ccp = CcpProperties(self.parameters.distance, self.parameters.length,
                                 self.parameters.V_bias,   self.parameters.frequency, self.parameters.phase,
                                 self.parameters.N_cells,  self.parameters.lateral_loss)

        if verbose: print('Reading neutral properties from file \"'+ FILENAME_NEUTRALS + '\"\n')
        self.neutrals = TargetParticles(self.parameters.N_sigma, self.parameters.N_sigma_ions,
                                        self.parameters.T_neutrals, self.parameters.p_neutrals,
                                        self.parameters.min_scattered, self.parameters.isactive_recomb,
                                        filename=FILENAME_NEUTRALS)
        (status, message) = self.neutrals.read_error
        if (status != 0):
                self.error = (status, message)
                return
        (status, message) = self.neutrals.read_properties(FILENAME_NEUTRALS, '\t',
                                                          self.parameters.e_min_sigma, self.parameters.e_max_sigma,
                                                          self.parameters.e_min_sigma_ions, self.parameters.e_max_sigma_ions,
                                                          debug=False)
        if (status != 0):
                self.error = (status, message)
                return

        # If there are not molecular gases, data about dissociation have not been saved and cannot be read
        if (self.neutrals.types_molecules > 0):
            self.read_neu_stat = True
        else:
            self.read_neu_stat = False

            
        # Read data saved during simulation
        generate_save_dir_name(self.parameters, abs_path=False)
        generate_save_file_names(self.parameters)

        self.parameters.filename_stat_ele += EXT
        self.parameters.filename_stat_neu += EXT
        self.parameters.filename_epos_z   += EXT
        self.parameters.filename_I        += EXT
        if self.read_edf:
            self.parameters.filename_distrib_ele += EXT
        if self.read_V:
            self.parameters.filename_V += EXT                        

        if verbose: print('Reading electron mean data from file \''+ self.parameters.filename_stat_ele +'\'\n')
        self.means_ele = data_manager.DataGrid()
        (status, message) = self.means_ele.read_file(self.parameters.filename_stat_ele, sep=SEP, transpose=True, skip=1)
        if (status != 0):
            string = 'Error reading file \"' + self.parameters.filename_stat_ele + '\": '+message
            self.error =  (status, string)
            return
        self.means_ion = []
        for i in range(self.neutrals.types):
            filename = self.parameters.filename_stat_ion + '_' + self.neutrals.names[i] + '+' + EXT
            if verbose: print('Reading ' + self.neutrals.names[i] + '+  mean data from file \''+ filename +'\'\n')
            self.means_ion.append(data_manager.DataGrid())
            (status, message) = self.means_ion[i].read_file(filename, sep=SEP, transpose=True, skip=1)
            if (status != 0):
                string = 'Error reading file \"' + filename + '\": '+message
                self.error =  (status, string)
                return
            
        if self.read_neu_stat:
            if verbose: print('Reading neutrals data from file \''+ self.parameters.filename_stat_neu +'\'\n')
            self.means_neu = data_manager.DataGrid()
            (status, message) = self.means_neu.read_file(self.parameters.filename_stat_neu, sep=SEP, transpose=True, skip=1)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_stat_neu + '\": '+message
                self.error =  (status, string)
                return
        
        if self.read_V:
            if verbose: print('Reading electric current data from file \''+ self.parameters.filename_I +'\'\n')
            self.means_I = data_manager.DataGrid()
            (status, message) = self.means_I.read_file(self.parameters.filename_I, sep=SEP, transpose=True, skip=1)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_I + '\": '+message
                self.error =  (status, string)
                return               
            if verbose: print('Reading electric potential data from file \''+ self.parameters.filename_V +'\'\n')
            self.V = data_manager.DataGrid()
            (status, message) = self.V.read_file(self.parameters.filename_V, sep=SEP)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_V + '\": '+message
                self.error =  (status, string)
                return

        if self.read_edf:
            if verbose: print('Reading EEDF data from file \''+ self.parameters.filename_distrib_ele +'\'\n')
            self.eedf = data_manager.DataGrid()
            (status, message) = self.eedf.read_file(self.parameters.filename_distrib_ele, sep=SEP, pad_value=numpy.nan)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_distrib_ele + '\": '+message
                self.error =  (status, string)
                return
            self.iedf = []
            for i in range(self.neutrals.types):
                filename = self.parameters.filename_distrib_ion + '_' + self.neutrals.names[i] + '+' + EXT
                if verbose: print('Reading IEDF data from file \''+ filename +'\'\n')
                #iedf = data_manager.DataGrid()
                self.iedf.append( data_manager.DataGrid() )
                (status, message) = self.iedf[i].read_file(filename, sep=SEP, pad_value=numpy.nan)
                if (status != 0):
                    string = 'Error reading file \"' + filename + '\": '+message
                    self.error =  (status, string)
                    return

        if self.read_z:
            if verbose: print('Reading electron z positions from file \''+ self.parameters.filename_epos_z +'\'\n')
            self.epos_z = data_manager.DataGrid()
            (status, message) = self.epos_z.read_file(self.parameters.filename_epos_z, sep=SEP, pad_value=numpy.nan)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_epos_z + '\": '+message
                self.error =  (status, string)
                return
            self.ipos_z = []
            for i in range(self.neutrals.types):
                filename = self.parameters.filename_ipos_z + '_' + self.neutrals.names[i] + '+' + EXT
                if verbose: print('Reading ion z positions file \''+ filename +'\'\n')
                self.ipos_z.append( data_manager.DataGrid() )
                (status, message) = self.ipos_z[i].read_file(filename, sep=SEP, pad_value=numpy.nan)
                if (status != 0):
                    string = 'Error reading file \"' + filename + '\": '+message
                    self.error =  (status, string)
                    return
                
        self.n_rows           = len(self.means_ele.data_array[0])
        self.n_rows_dist      = int( self.n_rows / self.parameters.save_delay_dist ) 
        self.max_time         = self.means_ele.data_array[0, self.n_rows-1]
        self.sim_duration     = self.means_ele.data_array[0, self.n_rows-1] * 1E-9
        self.output_timestep  = self.means_ele.data_array[0, self.n_rows-1] * 1E-9 / self.n_rows
        if self.read_neu_stat:
            self.n_cols           = len(self.means_neu.data_array)
            self.mol_types        = int((self.n_cols - 1) / 2)

        # lambda_D = A * sqrt(<E_e>**3 / n_e)
        # NOTE: ELECTRON_CHARGE is negative, multiplying by it is necessary to convert energy from eV to J
        self.debye_length     = (  self.DEBYE_CONST_A
                                 * numpy.sqrt(- ELECTRON_CHARGE
                                              * self.means_ele.data_array[3]
                                              / self.means_ele.data_array[14] )
                                )
        # N_D = 4/3 * pi * lamnda_D**3 * n_e
        self.debye_number     = 4/3 * numpy.pi * self.debye_length**3 * self.means_ele.data_array[14]
        # f_pla ~ sqrt(n_e)
        self.plasma_frequency = self.DEBYE_CONST_C * numpy.sqrt(self.means_ele.data_array[14])

        if verbose: print('Data loading completed\n')                

        self.error = (status, message)

        
    # +----------------------------------------------+
    # | Functions that return information about data |
    # +----------------------------------------------+
                
    def get_row(self, time):
        for row in range(self.n_rows):
            if (self.means_ele.data_array[0][row] >= time): break
        return row

    def get_row_dist(self, time):
        if (time < self.means_ele.data_array[0][0] * self.parameters.save_delay_dist):
            return None
        else:
            row = self.get_row(time)
            return round( (row+1) / self.parameters.save_delay_dist - 1 ) 

    def get_time(self, row):
        if ( (row < 0) or (row > self.n_rows-1) ):
            return None
        else:
            return self.means_ele.data_array[0][row]

    def get_time_dist(self, row_dist):
        if ( (row_dist < 0) or (row_dist > self.n_rows_dist-1) ):
            return None
        else:
            row_means = (row_dist+1) * self.parameters.save_delay_dist - 1 
            return self.get_time(row_means)
   
    def is_molecule(self, neutral_index):
        neutral_index = int(neutral_index)
        if ( (neutral_index < 0) or (neutral_index > self.neutrals.types) ):
            return None            
        if (self.neutrals.molecule_type[neutral_index] == 'a'):
            return False
        else:
            return True
    
    def get_molecule_index(self, neutral_index):
        if ( self.is_molecule(neutral_index) is True ):
            molecule_index = 0
            for i in range(neutral_index):
                if (self.neutrals.molecule_type[i] != 'a'):
                    molecule_index += 1
            return molecule_index
        elif ( self.is_molecule(neutral_index) is False):
            print(self.neutrals.names[neutral_index] + ' is not a molecule')
            return None
        else:
            print('Invalid index')
            return None            

    # +---------------------------------------------+
    # | Functions that print information about data |
    # +---------------------------------------------+
        
    def print_parameters(self):
        print('')
        print('Electrodes distance                 : ' + print_unit(self.ccp.distance,'m'))
        print('Electrodes lateral length           : ' + print_unit(self.ccp.length,'m'))
        print('Plasma volume                       : ' + str(self.ccp.volume) + ' m**3')
        print('Number of PIC cells                 : ' + str(self.parameters.N_cells))
        print('Cell dimension                      : ' + print_unit(self.ccp.delta_grid, 'm', 4 ))
        print('Electric bias                       : ' + print_unit(self.ccp.V_peak,'V'))
        print('Electric bias frequency             : ' + print_unit(self.ccp.frequency,'Hz'))
        print('Electric bias phase (at t=0)        : ' + str(self.ccp.phase))
        print('')
        print('Gas temperature                     : ' + str(self.neutrals.temperature) + ' K')
        print('Total pressure                      : ' + str(self.neutrals.total_pressure) + ' Pa')
        print('Number of tabulated e- xsec values  : ' + str(self.neutrals.n_sigma))
        print('Number of tabulated ion xsec values : ' + str(self.neutrals.n_sigma_ions))
        print('')
        print('Maximum number of electrons         : ' + str(self.parameters.Nmax_particles))
        print('Required starting ionization degree : ' + print_exp(self.parameters.start_ion_deg,4))
        print('')
        print('Timestep [0=automatic]              : ' + print_unit(self.parameters.dt, 's'))
        print('Mean time between data acquisitions : ' + print_unit(self.output_timestep, 's'))
        print('Number of data acquisitions         : ' + str(self.n_rows-1))
        print('Overall simulated time              : ' + print_unit(self.sim_duration, 's'))
        print('Data saved every                    : ' + str(self.parameters.save_delay) + ' data acquisitions')
        print('Distributions saved every           : ' + str(self.parameters.save_delay_dist) + ' data saves')

        return (0, OK)

    def print_ions(self):
        for i in range(self.neutrals.types):
            print('neutral index: ' + str(i))
            print('ion index:     ' + str(i+1))
            print('name:          ' + self.neutrals.names[i])
            print('type:          ' + self.neutrals.molecule_type[i] + '\n')
                
        return (0, OK)   
        
    def print_molecules(self):
        for i in range(self.neutrals.types):            
            if self.is_molecule(i):
                print('neutral  index: ' + str(i))
                print('molecule index: ' + str( self.get_molecule_index(i)))
                print('name:           ' + self.neutrals.names[i] + '\n')
                
        return (0, OK)

    
    # +------------------------------------------+
    # | Electron parameters time evolution plots |
    # +------------------------------------------+
        
    def plot_electron_number(self, real=True, computational=True, line='None',
                             symbol_real='x', symbol_comp='+',
                             color_real='red', color_comp='blue'):
        """Plot the number of electrons as a function of the simulation time"""
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)        
        plt.title('Number of electrons')
        plt.xlabel('Time / ns')
        plt.ylabel('Number of electrons')
        if computational:
            plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[1], 
                         marker=symbol_comp, linestyle=line, color=color_comp, label=r'Computational $e^-$')
        if real:
            plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[1] * self.means_ele.data_array[2], 
                         marker=symbol_real, linestyle=line, color=color_real, label=r'Real $e^-$')
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)
        

    def plot_electron_weight(self,line='None', symbol='.', color='red'):
        """ Plot the electron weight as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Electron weight')
        plt.xlabel('Time / ns')
        plt.ylabel('Weight')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[2], 
                     marker=symbol, linestyle=line, color=color, label = r'Weight $e^-$')
        plt.legend()
        plt.grid()
        plt.show()

        return (0, OK)

    def plot_electron_density(self, line='None', symbol='.', color='red'):
        """ Plot the number density of electrons as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Electron density')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$n_e$ / $m^{-3}$')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[14], 
                     marker=symbol, linestyle=line, color=color, label = r'$n_e$')
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)

    def plot_plasma_frequency(self, line='None', symbol='.', color='red'):
        """ Plot the plasma frequency as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)        
        plt.title('Plasma frequency')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$\nu_P$ / $s^{-1}$')
        plt.semilogy(self.means_ele.data_array[0], self.plasma_frequency, 
                     marker=symbol, linestyle=line, color=color, label = r'$\nu_P$')
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)

    def plot_debye_length(self, line='None', symbol='.', color='red'):
        """ Plot the Debye length as a function of the simulation time """        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Debye length')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$\lambda_D$ / m')
        plt.semilogy(self.means_ele.data_array[0], self.debye_length, 
                     marker=symbol, linestyle=line, color=color, label = r'$\lambda_D$')
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)
    
        
    def plot_debye_number(self,line='None', symbol='.', color='red'):
        """ Plot the Debye number as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Debye number')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$N_D$')
        plt.semilogy(self.means_ele.data_array[0], self.debye_number, 
                     marker=symbol, linestyle=line, color=color, label = r'$N_D$')
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)
    
               
    def plot_electron_mean_energy(self, plot_sigma=False, plot_min=False, plot_max=False, semilog=False,
                                  line='None', symbol='.', color='red', ecolor='orange',
                                  maxcolor='cyan', mincolor='blue'):
        """ Plot the mean energy of electrons as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()      
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Electron mean energy')
        plt.xlabel('Time / ns')
        plt.ylabel(r'<$E_e$> / eV')
        if semilog:
                plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[3],
                             marker=symbol, linestyle=line, color=color, label=r'<$E_e$>')
        else:
            if plot_sigma:
                plt.errorbar(self.means_ele.data_array[0], self.means_ele.data_array[3], self.means_ele.data_array[4], 
                             marker=symbol, linestyle=line, color=color, ecolor=ecolor, label = r'<$E_e$>')
            else:
                plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[3], 
                         marker=symbol, linestyle=line, color=color, label = r'<$E_e$>')

            if plot_min: plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[5], 
                                  marker=symbol, linestyle=line, color=mincolor, label = r'$E_e^{min}$')
            if plot_max: plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[6], 
                                  marker=symbol, linestyle=line, color=maxcolor, label = r'$E_e^{max}$')
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)

        
    def plot_electron_angle(self, plot_sigma=True, plot_min=True, plot_max=True,
                            line='None', symbol='.', color='red', ecolor='orange',
                            maxcolor='cyan', mincolor='blue'):
        """ Plot the angle between z direction for electrons as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Mean angle between electron velocity and electric field')
        plt.xlabel('Time / ns')
        plt.ylabel(r'<$\theta$> / °')
        if plot_sigma:                        
            plt.errorbar(self.means_ele.data_array[0], self.means_ele.data_array[7], self.means_ele.data_array[8], 
                         marker=symbol, linestyle=line, color=color, ecolor=ecolor, label=r'<$\theta$>')
        else:
            plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[7], 
                     marker=symbol, linestyle=line, color=color, label = r'<$\theta$>')

        if plot_min: plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[9], 
                              marker=symbol, linestyle=line, color=mincolor, label = r'$\theta_{min}$')
        if plot_max: plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[10], 
                              marker=symbol, linestyle=line, color=maxcolor, label = r'$\theta_{max}$')
        plt.legend()
        plt.grid()                
        plt.show()
        
        return (0, OK)

        
    def plot_tau(self, plot_dt=False, line='None', symbol_dt='+', symbol_tau='x', color_dt='blue', color_tau='red'):
        """ Plot the mean time between collisions and the timestep as a function of the simulation time """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        if plot_dt:
            plt.title('Timestep and mean time between collisions')
        else:
            plt.title('Mean time between collisions')
        plt.xlabel('Time / ns')
        plt.ylabel('Time / fs')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[12], 
                     marker=symbol_tau, linestyle=line, color=color_tau,  label = r'$\tau$')
        if plot_dt:
            plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[11], 
                         marker=symbol_dt, linestyle=line, color=color_dt, label = 'dt')
        
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)

        
    def plot_collision_frequency(self, line='None', symbol='.', color='red'):
        """ Plot the mean collision frequency for electron collisions """
        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Collision frequency')
        plt.xlabel('Time / ns')
        plt.ylabel('f %')
        plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[13], 
                 marker=symbol, linestyle=line, color=color, label = 'f')
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)

    # +-------------------------------------+
    # | Ion parameters time evolution plots |
    # +-------------------------------------+

    def plot_ion_number(self, ion_type, ion_string=None, real=True, computational=True, line='None',
                        symbol_real='x', symbol_comp='+',
                        color_real='red', color_comp='blue'):
        """ Plot the number of ions as a function of the simulation time """
        
        ion_type = int(ion_type)
        if ((ion_type < 1) or (ion_type > self.neutrals.types)):
            status  = 1
            message = 'Ion type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            ion_type -= 1            
        if (ion_string is None): ion_string = self.neutrals.names[ion_type] + r'$^{+}$'
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Number of ' + ion_string +' ions')
        plt.xlabel('Time / ns')
        plt.ylabel('Number of ions')
        if computational:
            plt.semilogy(self.means_ion[ion_type].data_array[0],
                         self.means_ion[ion_type].data_array[1], 
                         marker=symbol_comp, linestyle=line, color=color_comp,
                         label='Computational ' + ion_string)
        if real:
            plt.semilogy( self.means_ion[ion_type].data_array[0],
                          self.means_ion[ion_type].data_array[1]
                        * self.means_ion[ion_type].data_array[2], 
                         marker=symbol_real, linestyle=line, color=color_real,
                         label='Real ' + ion_string)
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)
    

    def plot_ion_weight(self, ion_type, ion_string=None, line='None', symbol='.', color='red'):
        """ Plot the ion weight as a function of the simulation time """
        
        ion_type = int(ion_type)
        if ((ion_type < 1) or (ion_type > self.neutrals.types)):
            status  = 1
            message = 'Ion type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            ion_type -= 1
        if (ion_string is None): ion_string = self.neutrals.names[ion_type] + r'$^{+}$'
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)        
        plt.title(ion_string + ' weight')
        plt.xlabel('Time / ns')
        plt.ylabel('Weight')
        plt.semilogy(self.means_ion[ion_type].data_array[0],
                     self.means_ion[ion_type].data_array[2], 
                     marker=symbol, linestyle=line, color=color,
                     label = 'Weight ' + ion_string)
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)

    def plot_ion_density(self, ion_type, ion_string=None, line='None', symbol='.', color='red'):
        """ Plot the ion number density as a function of the simulation time """

        
        ion_type = int(ion_type)
        if ((ion_type < 1) or (ion_type > self.neutrals.types)):
            status  = 1
            message = 'Ion type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            ion_type -= 1
        if (ion_string is None): ion_string = self.neutrals.names[ion_type] + r'$^{+}$'        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title(ion_string + ' density')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$n_i$ / $m^{-3}$')
        plt.semilogy(self.means_ion[ion_type].data_array[0],
                     self.means_ion[ion_type].data_array[11], 
                     marker=symbol, linestyle=line, color=color, label=ion_string)
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)

        
    def plot_ion_mean_energy(self, ion_type, ion_string=None,
                             plot_sigma=False, plot_min=False, plot_max=False, semilog=False,
                             line='None', symbol='.', color='red', ecolor='orange',
                             maxcolor='cyan', mincolor='blue'):
        """ Plot the ion mean energy as a function of the simulation time """
        
        ion_type = int(ion_type)
        if ((ion_type < 1) or (ion_type > self.neutrals.types)):
            status  = 1
            message = 'Ion type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            ion_type -= 1
        if (ion_string is None): ion_string = self.neutrals.names[ion_type] + r'$^{+}$'        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title(ion_string + ' mean energy')
        plt.xlabel('Time / ns')
        plt.ylabel(r'<$E_i$> / eV')
        if semilog:
                plt.semilogy(self.means_ion[ion_type].data_array[0],
                             self.means_ion[ion_type].data_array[3],
                             marker=symbol, linestyle=line, color=color, label=r'<$E_i$>')
        else:
            if plot_sigma:
                plt.errorbar(self.means_ion[ion_type].data_array[0],
                             self.means_ion[ion_type].data_array[3],
                             self.means_ion[ion_type].data_array[4], 
                             marker=symbol, linestyle=line, color=color, ecolor=ecolor, label = r'<$E_i$>')
            else:
                plt.plot(self.means_ion[ion_type].data_array[0],
                         self.means_ion[ion_type].data_array[3], 
                         marker=symbol, linestyle=line, color=color, label = r'<$E_i$>')

            if plot_min: plt.plot(self.means_ion[ion_type].data_array[0],
                                  self.means_ion[ion_type].data_array[5], 
                                  marker=symbol, linestyle=line, color=mincolor, label = r'$E_i^{min}$')
            if plot_max: plt.plot(self.means_ion[ion_type].data_array[0],
                                  self.means_ion[ion_type].data_array[6], 
                                  marker=symbol, linestyle=line, color=maxcolor, label = r'$E_i^{max}$')
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)
        

    def plot_ion_angle(self, ion_type, ion_string=None,
                       plot_sigma=True, plot_min=True, plot_max=True,
                       line='None', symbol='.', color='red', ecolor='orange',
                       maxcolor='cyan', mincolor='blue'):
        """ Plot the angle with z direction for ions as a function of the simulation time """
        
        ion_type = int(ion_type)
        if ((ion_type < 1) or (ion_type > self.neutrals.types)):
            status  = 1
            message = 'Ion type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            ion_type -= 1
        if (ion_string is None): ion_string = self.neutrals.names[ion_type] + r'$^{+}$'                
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Mean angle between ' + ion_string + ' velocity and electric field')
        plt.xlabel('Time / ns')
        plt.ylabel(r'<$\theta$> / °')
        if plot_sigma:                        
            plt.errorbar(self.means_ion[ion_type].data_array[0],
                         self.means_ion[ion_type].data_array[7],
                         self.means_ion[ion_type].data_array[8], 
                         marker=symbol, linestyle=line, color=color, ecolor=ecolor, label=r'<$\theta$>')
        else:
            plt.plot(self.means_ion[ion_type].data_array[0],
                     self.means_ion[ion_type].data_array[7], 
                     marker=symbol, linestyle=line, color=color, label = r'<$\theta$>')

        if plot_min: plt.plot(self.means_ion[ion_type].data_array[0],
                              self.means_ion[ion_type].data_array[9], 
                              marker=symbol, linestyle=line, color=mincolor, label = r'$\theta_{min}$')
        if plot_max: plt.plot(self.means_ion[ion_type].data_array[0],
                              self.means_ion[ion_type].data_array[10], 
                              marker=symbol, linestyle=line, color=maxcolor, label = r'$\theta_{max}$')
        plt.legend()
        plt.grid()                
        plt.show()
        
        return (0, OK)

        
    # +---------------------------+
    # | Electric properties plots |
    # +---------------------------+
    
    def plot_current_density(self, absolute=False, line='None', symbol='.', color='red'):
        """ Plot the current density as a function of the simulation time """        

        if (not self.read_V):
            status  = 1
            message = 'Electric potential and current are not available'
            return (status, message)        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Current density')
        plt.xlabel('Time / ns')
        plt.ylabel(r'j / $A m^{-2}$')
        if absolute:
            plt.plot(self.means_I.data_array[0] / (self.ccp.length * self.ccp.length),
                     numpy.abs(self.means_I.data_array[1]), 
                     marker=symbol, linestyle=line, color=color, label = 'j')            

        else:
            plt.plot(self.means_I.data_array[0] / (self.ccp.length * self.ccp.length),
                     self.means_I.data_array[1], 
                     marker=symbol, linestyle=line, color=color, label = 'j')
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)


    def plot_potential(self, time, line='-', symbol='.', color='red'):
        """ Plot the electric potential as a function of z coordinate """        
        
        if (not self.read_V):
            status  = 1
            message = 'Electric potential and current are not available'
            return (status, message)        
        if ((time < 0) or (time > self.max_time)):
            status  = 2
            message = 'Time must be in range [0,' + str(self.max_time) + '] ns'
            return (status, message)
        row = self.get_row(time)
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Electric potential (t = ' + str(time) + ' ns)')
        plt.xlabel('z / mm')
        plt.ylabel(r'$\Delta V$ / V')
        plt.plot(self.ccp.grid_points*1.0E3, self.V.data_array[row], 
                 marker=symbol, linestyle=line, color=color)
        plt.grid()
        plt.show()
        
        return (0, OK)
        
    # +---------------------------------------+
    # | Dissociation rate/rate constant plots |
    # +---------------------------------------+
    
    def plot_dissociation_rate(self, neutral_index, neutral_string=None, line='None', symbol='.', color='red'):
        """ Plot the dissociation rate as a function of the simulation time """

        if ((neutral_index < 1) or (neutral_index > self.neutrals.types)):
            status  = 1
            message = 'Neutral type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            neutral_index -= 1        
        
        if (not self.is_molecule(neutral_index)):
            status  = 1
            message = self.neutrals.names[neutral_index] + ' is not a molecule'
            return (status, message)
        
        if (neutral_string is None): neutral_string = self.neutrals.names[neutral_index]
        
        molecule_index = self.get_molecule_index(neutral_index)
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Dissociation rate')
        plt.xlabel('Time / ns')
        plt.ylabel(r'R / $m^3$ $s^{-1}$')
        col = 1 + molecule_index * 2
        plt.semilogy(self.means_neu.data_array[0], self.means_neu.data_array[col],
                     marker=symbol, linestyle=line, color=color, label=neutral_string)        
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)
        

    def plot_dissociation_rate_const(self, neutral_index, neutral_string=None, line='None', symbol='.', color='red'):
        """ Plot the dissociation rate constant as a function of the simulation time """

        if ((neutral_index < 1) or (neutral_index > self.neutrals.types)):
            status  = 1
            message = 'Neutral type must be in range [1,' + str(self.neutrals.types) + ']'
            return (status, message)
        else:
            neutral_index -= 1     
        
        if (not self.is_molecule(neutral_index)):
            status  = 1
            message = self.neutrals.names[neutral_index] + ' is not a molecule'
            return (status, message)
        if (neutral_string is None): neutral_string = self.neutrals.names[neutral_index] + r'$^{+}$'
        molecule_index = self.get_molecule_index(neutral_index)
        plt.ioff()        
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Dissociation rate constant')
        plt.xlabel('Time / ns')
        plt.ylabel(r'k / $m^3$ $s^{-1}$')
        col = 2 + molecule_index * 2
        plt.plot(self.means_neu.data_array[0],self.means_neu.data_array[col],
                 marker=symbol, linestyle=line, color=color, label=neutral_string)         
        plt.grid()
        plt.legend()
        plt.show()
        
        return (0, OK)
            

    # +------------------+
    # | Z position plots |
    # +------------------+
                
    def plot_position(self, index, name_string=None, time=None, row=None, log=False, color='red', symbol='.'):
        """ Plot the particle z coordinate as a function of the simulation time """        
        
        if (not self.read_z):
            status  = 1
            message = 'Positions are not available'
            return (status, message)
        elif ((time is None) and (row is None)):
            status  = 2
            message = 'Give either simulation time or data row number !'
            return (status, message)
        elif (time is None):
            time = self.means_ele.data_array[0][row]            
        elif ((time < 0) or (time > self.max_time)):
            status  = 3
            message = 'Time must be in range [0,' + str(self.max_time) + '] ns'
            return (status, message)
        else:
            row = self.get_row(time)        
        string = ': row = ' + str(row) + "; time = " + str(time) + ' ns'        
        if (index == 0):
            if (name_string is None): name = r'$e^{-}$'
            ensamble = self.epos_z.data_array[row] * 1.0e3
        else:
            if (name_string is None): name = self.neutrals.names[index-1] + r'$^{+}$ '
            ensamble = self.ipos_z[index-1].data_array[row] * 1.0e3
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Ensamble of ' + name)
        plt.xlabel('Particle index')
        plt.ylabel('z / mm')
        if log:
            plt.semilogy(ensamble, marker=symbol, linestyle='None', color=color, label=name+string)
        else:    
            plt.plot(ensamble, marker=symbol, linestyle='None', color=color, label=name+string)
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)


    def plot_z_distribution(self, particle_index, time, name_string=None, bins=None, density=True, hist_type='step', color='red'):
        """ Plot the particle distribution along the z axis """
        
        if ((particle_index < 0) or (particle_index > self.neutrals.types)):
             status  = 1
             message = 'Particle index must be in range [0,' + str(self.neutrals.types) + ']'
             return (status, message)
        if (not self.read_z):
            status  = 2
            message = 'Positions are not available'
            return (status, message)        
        if ((time < 0) or (time > self.max_time)):
            status  = 3
            message = 'Time must be in range [0,' + str(self.max_time) + '] ns'
            return (status, message)
        row = self.get_row(time)
        if (particle_index == 0):
            if (name_string is None): name = r'$e^{-}$'
            ensamble = self.epos_z.data_array[row] * 1.0e3
        else:
            if (name_string is None): name = self.neutrals.names[particle_index-1] + r'$^{+}$ '
            ensamble = self.ipos_z[particle_index-1].data_array[row] * 1.0e3
        if (bins is None):
            bins = self.ccp.grid_points * 1.0E3
            #bins = self.parameters.N_cells
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title(name + ' spatial distribution (t = ' + str(time) + ' ns)')
        plt.xlabel('z / mm')
        if density:
            plt.ylabel(r'Particle density / $mm^{-3}$')
        else:
            plt.ylabel(r'Number of particles')
        plt.hist(ensamble, bins=bins, density=density, histtype=hist_type, color=color)
        plt.grid()
        plt.show()
        
        return (0, OK)
    
        
    # +---------------------------+
    # | Energy distribution plots |
    # +---------------------------+
                
    def plot_ensamble(self, index, time=None, row=None, log=False, color='red', symbol=','):
        """ Plot the particle energy """
        
        if (not self.read_edf):
            status  = 1
            message = 'Energy distribution functions are not available'
            return (status, message)
        elif ((time is None) and (row is None)):
            status  = 2
            message = 'Give either simulation time or data row number !'
            return (status, message)
        elif (time is None):
            time = self.means_ele.data_array[0][row]
        elif ((time < 0) or (time > self.max_time)):
            status  = 3
            message = 'Time must be in range [0,' + str(self.max_time) + '] ns'
            return (status, message)
        else:
            row = self.get_row(time)                        
        string = 'row #' + str(row) + "; time = " + str(time) + ' ns'
        if (index == 0):
            name = 'e- '
            ensamble = self.eedf.data_array[row]
        else:
            name     = self.neutrals.names[index-1] + '+ '
            ensamble = self.iedf[index-1].data_array[row]            
        plt.ioff()
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        plt.title('Ensamble of ' + name)
        plt.xlabel('Particle index')
        plt.ylabel('Energy / eV')
        if log:
            plt.semilogy(ensamble, marker=symbol, linestyle='None', color=color, label=name+string)
        else:    
            plt.plot(ensamble, marker=symbol, linestyle='None', color=color, label=name+string)
        plt.legend()
        plt.grid()
        plt.show()
        
        return (0, OK)
    

    def plot_edf(self, particle_index, time, name_string=None, bins='fd',
                 e_min=None, e_max=None, plot_inset=False, hist_type='step', color='red',
                 major_grid=False, minor_grid=False):
        """ Plot EEDF or IEDF """
        
        if not self.read_edf:
            status  = 1 
            message = 'Energy distribution functions are not available'
            return (status, message)
        particle_index=int(particle_index)
        if ( (particle_index < 0) or (particle_index > self.neutrals.types) ):
            status  = 2
            message = 'Particle index must be in range [0,' + str(self.neutrals.types) + ']'
            return (status, message)
        if ((time < 0) or (time > self.max_time)):
            status  = 3
            message = 'Time must be in range [0,' + str(self.max_time) + '] ns'
            return (status, message)        
        row = self.get_row(time) 
        #print('\nRow # ' + str(row) + ' of ' + str(self.n_rows-1))        
        if (particle_index==0):
            if (name_string is None): name_string= r'$e^{-}$'
            title_string = 'EEDF (time = ' + str(time) + ' ns)'
            x = self.eedf.data_array[row]
        else:
            ion_type = particle_index-1
            if (name_string is None): name_string = self.neutrals.names[particle_index-1]
            title_string = 'IEDF ' + name_string + r'$^{+}$ (time = ' + str(time) + ' ns)'            
            x = self.iedf[ion_type].data_array[row]
        edf = x[~numpy.isnan(x)]

        # Prepare main figure
        plt.figure()
        plt.title(title_string)
        plt.xlabel('Energy / eV')
        plt.ylabel(r'Probability density / $eV^{-1}$')

        # Plot main histogram
        plt.hist(edf, bins=bins, density=True, histtype=hist_type, color=color)
        plt.minorticks_on()
        plt.tick_params(which='major', direction='in', length=self.tick_length_major)
        plt.tick_params(which='minor', direction='in', length=self.tick_length_minor)
        if major_grid: plt.tick_params(axis='x', which='major', grid_linestyle='-',  grid_color='grey')
        if minor_grid: plt.tick_params(axis='x', which='minor', grid_linestyle='--', grid_color='grey')
        plt.ylim(-0.0025)
        if (major_grid and minor_grid):
            plt.grid(axis='x', which='both')
        elif major_grid:
            plt.grid(axis='x', which='major')
        elif minor_grid:
            plt.grid(axis='x', which='minor')
        if plot_inset:
            if (e_min is None): e_min = - edf.min()
            if (e_max is None): e_max = edf.max() / 4
            # Prepare inset
            ax2 = inset_locator.inset_axes(plt.gca(), width='50%', height='60%')
            # Plot secondary histogram in inset
            ax2.hist(edf, bins=bins, density=True, histtype=hist_type, color=color)
            ax2.minorticks_on()
            ax2.tick_params(which='major', direction='in', length=self.tick_length_major)
            ax2.tick_params(which='minor', direction='in', length=self.tick_length_minor)
            ax2.tick_params(axis='x', which='major', grid_linestyle='-',  grid_color='grey')
            ax2.tick_params(axis='x', which='minor', grid_linestyle='--', grid_color='grey')
            ax2.set_xlim(left=e_min, right=e_max)
            ax2.set_ylim(bottom=-0.01)
            ax2.grid(axis='x', which='both')
        plt.show()             

        return (0, OK)

        
    def plot_eedf_un(self, time=None, row=None, pdf='Maxwell', intervals=0, int_method='sqrt', method='fixed', plot_interface='pylab'):
        """ Compare the particle distribution with a distribution fuction. 

            A histogram will be created, the number of bins can be given or automatically calculated

            Parameters
            ----------

            time:                   simulation time at which the eedf/iedf is requested
            row:                    data row containint the eedf/iedf requested
                                    must be given if time is not
            pdf:                    probability distribution function describing the model
                                   'Maxwell' -> Maxwell pdf for particle kinetic energy

            intervals:              maximum allowed number of intervals (0 means no limit given)
            int_method:             method used to calculate the maximum number of intervals, used only if intervals=0
                                    'sqrt':  n_bins = sqrt(n_points)
                                    'log2':  n_bins = 1 + log2(n_points)
                                    'root3': n_bins = n_points**(1/3)
                                    'norm':  n_bins = n_data / (3.5 * stdv / n_points**(1/3) )
            method:                 method used to calculate the intervals
                                    'fixed': a fixed number of equiparted intervals is used
                                    'tails': width of each interval is a multiple of self.width/sqrt(self.n_data)
                                    expected frequency for each interval (except the last one) 
                                    is not less than 5                      
                                    'iterate': all intervals have the same width
                                    none of the intervals has zero frequency 
                                    not more than 20% of the intervals has a frequency lower than 5
            plot_interface:         graphic interface to be used to plot the histogram
                                    'pylab'   -> use matplotlib
                                    'gnuplot' -> use gnuplot

        """
        
        if not self.read_edf:
            status  = 1
            message = 'Energy distribution functions are not available'
            return (status, message)
        
        if (time is None):
                if (row is None):
                        status =  2
                        message = 'Give either simulation time or data row number !'
                        return (status, message)
        else:
                row = self.get_row(time) 
        print('\nRow # ' + str(row) + ' of ' + str(self.n_rows-1))
        self.n_electrons = int(self.means_ele.data_array[1, row])
        self.h= univariate.DataSet(self.eedf.data_array[row][0:self.n_electrons])
        if (pdf == 'Maxwell'):
                self.kt = 2.0 * self.h.mean / 3.0
                self.pdf = lambda x: pdf_maxwell_energy(x, self.kt)
        else:
                print('Unknown distribution type !')
                return
        print(self.h.expected_frequencies(self.pdf, intervals=intervals, int_method=int_method, method=method))
        print(self.h.observed_frequencies())
        print(self.h.chisquare_estimation())
        print('Time        = ' + print_unit(self.means_ele.data_array[0, row]*1E-9,'s'))
        print('Electrons   = ' + str(self.n_electrons))
        print('Weight      = ' + str(self.means_ele.data_array[2, row]))
        print('Mean energy = ' + print_unit(self.means_ele.data_array[3, row],  'eV', 3))
        print('Stdv energy = ' + print_unit(self.means_ele.data_array[4, row],  'eV', 3))
        print('Min energy  = ' + print_unit(self.means_ele.data_array[5, row],  'eV', 3))
        print('Max energy  = ' + print_unit(self.means_ele.data_array[6, row],  'eV', 3))
        print('Mean angle  = ' + print_unit(self.means_ele.data_array[7, row],  'deg', 3))
        print('Stdv angle  = ' + print_unit(self.means_ele.data_array[8, row],  'deg', 3))
        print('Min angle   = ' + print_unit(self.means_ele.data_array[9, row],  'deg', 3))
        print('Max angle   = ' + print_unit(self.means_ele.data_array[10, row],         'deg', 3))
        print('dt          = ' + print_unit(self.means_ele.data_array[11,row]*1E-15, 's', 3))
        print('tau         = ' + print_unit(self.means_ele.data_array[12,row]*1E-15, 's', 3))
        print('p           = ' + str(self.means_ele.data_array[13,row]) + ' %')
#       print 'E_mean      = ' + print_unit(h.mean,             'eV', 3)
#       print 'E_max       = ' + print_unit(h.max,              'eV', 3)
#       print 'E_min       = ' + print_unit(h.min,              'eV', 3)
        print('k*T         = ' + print_unit(self.kt,            'eV', 3))
        print('N           = ' + str(self.h.n_data))
        print('DF          = ' + str(self.h.freedom_degrees))
        print('Chisquare   = ' + str(self.h.chisquare))
        print('P-value     = ' + str(self.h.p_value))
        print(self.h.plot_histogram(interface=plot_interface))
        
        return (0, OK)
        
