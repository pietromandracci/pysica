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

""" Read and show data saved by ccpla script """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

# Mudules from the standard Python library
import math

# Modules provided by the Python community
import numpy
import matplotlib.pyplot as plt
import Gnuplot

# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *
#from pysica.plasmapro.plot.plot_parameters import *

# Modules provided by plasmapro package
from pysica.managers import data_manager
from pysica.managers.unit_manager import print_unit, print_exp
from pysica.analysis import univariate
from pysica.functions.pdf import pdf_maxwell_energy
from pysica.plasmapro.discharge.reactors import CcpProperties
from pysica.plasmapro.discharge.target_particles import TargetParticles
from pysica.plasmapro.ccpla_init import *

plt.rcParams["mathtext.default"] = 'regular'
plt.rcParams["font.size"]        = 20


class CcplaSavedData:
    """Class of data saved by ccpla program"""

    def __init__(self, name=None, verbose=False, read_edf=False, read_potential=False):
        """Initialize the collection of data to be analyzed"""

        # Reading EEDF and IEDF as well as potential distribution can be suppressed to save memory and execution time
        self.read_edf = read_edf
        self.read_V   = read_potential

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

        # Read data saved during simulation
        generate_save_dir_name(self.parameters, abs_path=False)
        generate_save_file_names(self.parameters)

        self.parameters.filename_stat_ele +=  EXT
        self.parameters.filename_stat_neu +=  EXT
        self.parameters.filename_I        +=  EXT
        if self.read_edf:
            self.parameters.filename_distrib_ele +=  EXT
        if self.read_V:
            self.parameters.filename_V +=  EXT                        

        if verbose: print('Reading electron mean data from file \''+ self.parameters.filename_stat_ele +'\'\n')
        self.means_ele = data_manager.DataGrid()
        (status, message) = self.means_ele.read_file(self.parameters.filename_stat_ele, sep=SEP, transpose=True, skip=1)
        if (status != 0):
            string = 'Error reading file \"' + self.parameters.filename_stat_ele + '\": '+message
            self.error =  (status, string)
            return

        if verbose: print('Reading electron neutrals data from file \''+ self.parameters.filename_stat_neu +'\'\n')
        self.means_neu = data_manager.DataGrid()
        (status, message) = self.means_neu.read_file(self.parameters.filename_stat_neu, sep=SEP, transpose=True, skip=1)
        if (status != 0):
            string = 'Error reading file \"' + self.parameters.filename_stat_neu + '\": '+message
            self.error =  (status, string)
            return

        if verbose: print('Reading electric current data from file \''+ self.parameters.filename_I +'\'\n')
        self.means_I = data_manager.DataGrid()
        (status, message) = self.means_I.read_file(self.parameters.filename_I, sep=SEP, transpose=True, skip=1)
        if (status != 0):
            string = 'Error reading file \"' + self.parameters.filename_I + '\": '+message
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

        if self.read_V:
            if verbose: print('Reading electric potential data from file \''+ self.parameters.filename_V +'\'\n')
            self.V = data_manager.DataGrid()
            (status, message) = self.V.read_file(self.parameters.filename_V, sep=SEP)
            if (status != 0):
                string = 'Error reading file \"' + self.parameters.filename_V + '\": '+message
                self.error =  (status, string)
                return
                

        self.n_rows           = len(self.means_ele.data_array[0])
        self.sim_duration     = self.means_ele.data_array[0, self.n_rows-1] * 1E-9
        self.output_timestep  = self.means_ele.data_array[0, self.n_rows-1] * 1E-9 / self.n_rows
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

    def get_row(self, time):
        for row in range(self.n_rows):
            if (self.means_ele.data_array[0][row] >= time): break
        return row
    
    # +----------------------------------------+
    # | Plasma parameters time evolution plots |
    # +----------------------------------------+
        
    def plot_electron_number(self, real=True, computational=True, line='None', symbol='.',
                             color_real='red', color_comp='blue'):
        plt.ioff()
        plt.title('Number of electrons')
        plt.xlabel('Time / ns')
        plt.ylabel('Number of electrons')
        if computational:
            plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[1], 
                         marker=symbol, linestyle=line, color=color_comp, label=r'Computational $e^-$')
        if real:
            plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[1] * self.means_ele.data_array[2], 
                         marker=symbol, linestyle=line, color=color_real, label=r'Real $e^-$')
        plt.legend()
        plt.grid()
        plt.show()

    def plot_electron_weight(self,line='None', symbol='.', color='red'):
        plt.ioff()
        plt.title('Electron weight')
        plt.xlabel('Time / ns')
        plt.ylabel('Weight')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[2], 
                     marker=symbol, linestyle=line, color=color, label = r'Weight $e^-$')
        plt.legend()
        plt.grid()
        plt.show()     

    def plot_electron_density(self, line='None', symbol='.', color='red'):
        plt.title('Electron density')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$n_e$ / $m^{-3}$')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[14], 
                     marker=symbol, linestyle=line, color=color, label = r'$n_e$')
        plt.grid()
        plt.legend()
        plt.show()

    def plot_plasma_frequency(self,line='None', symbol='.', color='red'):
        plt.title('Plasma frequency')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$\nu_P$ / $s^{-1}$')
        plt.semilogy(self.means_ele.data_array[0], self.plasma_frequency, 
                     marker=symbol, linestyle=line, color=color, label = r'$\nu_P$')
        plt.grid()
        plt.legend()
        plt.show()        

    def plot_debye_length(self,line='None', symbol='.', color='red'):
        plt.title('Debye length')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$\lambda_D$ / m')
        plt.semilogy(self.means_ele.data_array[0], self.debye_length, 
                     marker=symbol, linestyle=line, color=color, label = r'$\lambda_D$')
        plt.grid()
        plt.legend()
        plt.show()
        
    def plot_debye_number(self,line='None', symbol='.', color='red'):
        plt.title('Debye number')
        plt.xlabel('Time / ns')
        plt.ylabel(r'$N_D$')
        plt.semilogy(self.means_ele.data_array[0], self.debye_number, 
                     marker=symbol, linestyle=line, color=color, label = r'$N_D$')
        plt.grid()
        plt.legend()
        plt.show()    
               
    def plot_electron_mean_energy(self, plot_sigma=False, plot_min=False, plot_max=False, semilog=False,
                                  line='None', symbol='.', color='red', ecolor='orange',
                                  maxcolor='cyan', mincolor='blue'):
        plt.ioff()
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

        
    def plot_electron_angle(self, plot_sigma=True, plot_min=True, plot_max=True,
                            line='None', symbol='.', color='red', ecolor='orange',
                            maxcolor='cyan', mincolor='blue'):
        plt.ioff()
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

        
    def plot_tau(self, line='None', symbol='.', color='red'):
        plt.ioff()
        plt.title('Timestep and time between collisions')
        plt.xlabel('Time / ns')
        plt.ylabel('Time / fs')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[11], 
                     marker=symbol, linestyle=line, color=color, label = 'dt')
        plt.semilogy(self.means_ele.data_array[0], self.means_ele.data_array[12], 
                     marker=symbol, linestyle=line, color=color,  label = r'$\tau$')
        plt.legend()
        plt.grid()
        plt.show()

        
    def plot_collision_frequency(self, line='None', symbol='.', color='red'):
        plt.ioff()
        plt.title('Collision frequency')
        plt.xlabel('Time / ns')
        plt.ylabel('f %')
        plt.plot(self.means_ele.data_array[0], self.means_ele.data_array[13], 
                 marker=symbol, linestyle=line, color=color, label = 'f')
        plt.legend()
        plt.grid()
        plt.show()
        
    # +---------------------------+
    # | Electric properties plots |
    # +---------------------------+
    
    def plot_current_density(self, absolute=False, line='None', symbol='.', color='red'):
        plt.ioff()
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


    def plot_potential(self,time, line='-', symbol='.', color='red'):
        row = self.get_row(time)
        plt.ioff()
        plt.title('Electric potential')
        plt.xlabel('z / mm')
        plt.ylabel(r'$\Delta V$ / V')
        plt.plot(self.ccp.grid_points*1.0E3, self.V.data_array[row], 
                 marker=symbol, linestyle=line, color=color)
        plt.grid()
        plt.show()           
        
    # +---------------------------------------+
    # | Dissociation rate/rate constant plots |
    # +---------------------------------------+
        
    def plot_dissocation_rates(self, line='None', symbol='.', color='red'):
        plt.title('Dissociation rates')
        plt.xlabel('Time / ns')
        plt.ylabel('R / m**3 s**-1')
        i_name = 0
        for i in range(self.mol_types):
            col = 1 + i*2
            # search the name of the molecule
            for j in range(i_name, self.neutrals.types):
                if (self.neutrals.molecule_type[j] != 'a'):
                    i_name = j
                    break
            name = self.neutrals.names[i_name]
            i_name += 1
            plt.semilogy(self.means_neu.data_array[0], self.means_neu.data_array[col],
                         marker=symbol, linestyle=line, color=color, label = name)
        plt.grid()
        plt.legend()
        plt.show()

        
    def plot_dissocation_rate_const(self, line='None', symbol='.', color='red'):
        plt.title('Dissociation rate constants')
        plt.xlabel('Time / ns')
        plt.ylabel('k / m**3 s**-1')
        i_name = 0                
        for i in range(self.mol_types):
            col = 1 + i*2 +1
            # search the name of the molecule
            for j in range(i_name, self.neutrals.types):
                if (self.neutrals.molecule_type[j] != 'a'):
                    i_name = j
                    break
            name = self.neutrals.names[i_name]
            i_name += 1                        
            plt.semilogy(self.means_neu.data_array[0],self.means_neu.data_array[col],
                         marker=symbol, linestyle=line, color=color, label = name)
        plt.grid()
        plt.legend()
        plt.show()

    # +---------------------------+
    # | Energy distribution plots |
    # +---------------------------+
                
    def plot_ensamble(self, index, time=None, row=None, log=False, color='red', symbol=','):
        if not self.read_edf:
            print('Distribution functions are not available')
            return
        if (time is None):
            if (row is None):
                print('Give either simulation time or data row number !')
                return
        else:
            row = self.get_row(time)
        time = self.means_ele.data_array[0][row]                
        string = 'row #' + str(row) + "; time = " + str(time) + ' ns'
        if (index == 0):
            name = 'e- '
            ensamble = self.eedf.data_array[row]
        else:
            name     = self.neutrals.names[index-1] + '+ '
            ensamble = self.iedf[index-1].data_array[row]
        plt.ioff()
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
                
                
    def plot_eedf(self, time=None, row=None, pdf='Maxwell', intervals=0, int_method='sqrt', method='fixed', plot_interface='pylab'):
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
            print('Distribution functions are not available')
            return
        
        if (time is None):
                if (row is None):
                        print('Give either simulation time or data row number !')
                        return
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
        
