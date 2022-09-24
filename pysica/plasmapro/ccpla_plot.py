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

""" Simulation of a capacitively coupled plasma discharge: functions used to plot data in real time """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

import math
import numpy


# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *
from pysica.managers.gnuplot_manager import *


#+---------------------------+
#| General purpose functions |
#+---------------------------+

def _decimate(a, k=MAX_DECIMATION_FACTOR):
    """ Keeps only an item every k ones from an array 

        Parameters
        ----------

        a:  original array
        k : decimation factor        

        Returns
        -------
        decimated array
    """
    
    if (k < 2): k=2
    d = numpy.zeros(len(a))
    n = int(len(a)/k)
    for i in range(n):
        d[i] = a[int(i*k)]
    return d,n


def decimate(a, n=MAX_DECIMATION_FACTOR):
    """ Removes an item every n ones from an array

        Parameters
        ----------

        a: original array
        n: decimation factor

        Returns
        -------
        decimated array
    """
    
    if (n < 2)       : n = 2
    if (n > len(a)-1): n = len(a)-1
    d = numpy.zeros(len(a))
    j = 0 # index on d
    k = 0 # decimation counter
    for i in range(1, len(a)):
        k += 1
        if (k < n):
            d[j] = a[i]
            j += 1                        
        else:
            k =  0
    return d, j



#+---------------+
#| Runtime plots |
#+---------------+
    
class CcplaPlots:
        
    def __init__(self, parameters, charges, ccp, screen_width=SCREEEN_WIDTH, screen_height=SCREEEN_HEIGHT,
                 vertical_shift=VERTICAL_SHIFT, hborder=HBORDER, vborder=VBORDER):
        """ Defines an istance of the class CcplaPlots, which contains all the runtime plots

            The ensamle of plots is initialized (dimensions, positions etc...)
        """
        
        # Define dimensions of plot windows
        self.window_width          = int( screen_width / 4 - hborder)
        self.window_height         = int( (screen_height - 2*vertical_shift) / 3 - vborder )
        self.vertical_shift        = vertical_shift
        self.hborder               = hborder
        self.vborder               = vborder
        #self.dot_points            = parameters.dot_points
        self.ccp                   = ccp
        if parameters.dot_points: self.plot_style = 'dots'
        else:                     self.plot_style = 'points'

        # Type of gnuplot terminal
        self.str_terminal          = DEFAULT_TERMINAL
        
        self.initialize_graphs()

        # Length of history arrays after which they are decimated
        self.history_limit           = parameters.n_max_points + 1
        # Decimation factor
        self.history_decimation      = parameters.decimation_factor

        # Define numpy arrays that will be used to plot average electron energies during simulation
        self.history_time            = numpy.zeros(self.history_limit,'d')      # Elapsed simulation time
        self.history_e_average       = numpy.zeros(self.history_limit,'d')      # Mean electrons kinetic energy
        self.history_n_electrons     = numpy.zeros(self.history_limit,'i')      # Number of electrons
        self.history_nw_electrons    = numpy.zeros(self.history_limit,'d')      # Number of electrons, weighted

        self.reset_history()

        
    def plot_time(self, graph, time):
        """ Prints simulation time on a graph 
             
            Parameters
            ----------
            graph:  graph on which the time must be printed
                    it must have been created using the new_plot function
            time:   value of the simulation time

            Returns
            -------
            graph_x: absolute position of the graph window in screen coordinates

        """        
        string_time = 't = ' + str(round(time*1E9,3)) + ' ns'
        plot_label(graph, x=LABEL_XPOS, y=LABEL_YPOS, label=string_time, erase=True)    

        
    def graph_x(self, x):
        """ Trasforms logical x coordinate in the absolute coordinate to use for graph positioning 
             
            Parameters
            ----------
            x:  x logical coordinate: 
                x=0 means first graph from the left, 
                x=1 means second graph from the left ....

            Returns
            -------
            graph_x: absolute position of the graph window in screen coordinates

        """
        return x * (self.window_width + self.hborder)

    
    def graph_y(self, y):
        """ Trasforms logical y coordinate in the absolute coordinate to use for graph positioning 
             
            Parameters
            ----------
            y:  y logical coordinate: 
                y=0 means first graph from the top, 
                y=1 means second graph from the top ....

            Returns
            -------
            graph_y: absolute position of the graph window in screen coordinates

        """        
        return self.vertical_shift + y *(self.window_height + self.vborder)
                
               
    def initialize_graphs(self, reset=False):
        """ Creates and initializes all the graphs 
              
            Parameters
            ----------

            reset: if set to True, 
                   - all open graphs are closed before opening the new ones
                   - the history is reset
        """

        # If requested, clears all open graphs and reset history
        if reset:
            # Close all open graphs
            plot_close_all()
            # Reset numpy arrays that will be used to plot average electron energies during simulation          
            self.reset_history()
           
        # Number of electrons vs time
        self.n_electrons_graph              = new_plot(xpos=self.graph_x(0), ypos=self.graph_y(0),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       title='Number of electrons',
                                                       xlabel='Time / ns', ylabel='Electrons',
                                                       xmin=0,
                                                       logy=True,
                                                       redirect_output=True,
                                                       purge=False)
        # Mean energy of electrons vs time
        self.electrons_average_energy_graph = new_plot(xpos=self.graph_x(1), ypos=self.graph_y(0),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       title='Electron mean energy',
                                                       xlabel='Time / ns', ylabel='Energy / eV',
                                                       redirect_output=True,
                                                       purge=False)
        # Charge density vs z-position
        self.charge_density_graph           = new_plot(xpos=self.graph_x(2), ypos=self.graph_y(0),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0,
                                                       xmax=self.ccp.distance*1E3,
                                                       format_y='%.2g',
                                                       title='Charge density',
                                                       xlabel='z / mm', ylabel='Rho / C  m**-3',
                                                       redirect_output=True,
                                                       purge=False)
        # Electric potential vs z-position
        self.potential_graph                = new_plot(xpos=self.graph_x(3), ypos=self.graph_y(0),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0,
                                                       xmax=self.ccp.distance*1E3,
                                                       title='Electric potential',
                                                       xlabel='z / mm', ylabel='Potential / V',
                                                       redirect_output=True,
                                                       purge=False)
        # Electron energy vs angle
        self.electrons_phase_space_graph    = new_plot(xpos=self.graph_x(0), ypos=self.graph_y(1),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=-5,
                                                       xmax=185,
                                                       format_y='%.2g',
                                                       style=self.plot_style,
                                                       title='Electron energy vs angle',
                                                       xlabel='Angle / deg', ylabel='Energy / eV',
                                                       redirect_output=True,
                                                       purge=False)
        # Ion energy vs angle
        self.ions_phase_space_graph         = new_plot(xpos=self.graph_x(1), ypos=self.graph_y(1),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=-5,
                                                       xmax=185,
                                                       format_y='%.2g',
                                                       style=self.plot_style,
                                                       title='Ion energy vs angle',
                                                       xlabel='Angle / deg', ylabel='Energy / eV',
                                                       redirect_output=True,
                                                       purge=False)
        # Electron energy vs z-position
        self.electrons_phase_space_graph2   = new_plot(xpos=self.graph_x(2), ypos=self.graph_y(1),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0,
                                                       xmax=self.ccp.distance*1E3,
                                                       format_y='%#.2g',
                                                       style=self.plot_style,
                                                       title='Electron energy vs z-position',
                                                       xlabel='z / mm', ylabel='Energy / eV',
                                                       redirect_output=True,
                                                       purge=False)
        # Ion energy vs z-position
        self.ions_phase_space_graph2        = new_plot(xpos=self.graph_x(3), ypos=self.graph_y(1),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0,
                                                       xmax=self.ccp.distance*1E3,
                                                       format_y='%#.2g',
                                                       style=self.plot_style,
                                                       title='Ion energy vs z-position',
                                                       xlabel='z / mm', ylabel='Energy / eV',
                                                       redirect_output=True,
                                                       purge=False)
        # Electron energy distribution
        self.eedf_graph                     = new_plot(xpos=self.graph_x(0), ypos=self.graph_y(2),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       style='histeps',
                                                       title='Electron energy distribution',
                                                       xlabel='Energy / eV', ylabel='Number of electrons / %',
                                                       redirect_output=True,
                                                       purge=False)
        # Ion energy distribution
        self.iedf_graph                     = new_plot(xpos=self.graph_x(1), ypos=self.graph_y(2),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       style='histeps',
                                                       title='Ion energy distribution',
                                                       xlabel='Energy / eV', ylabel='Number of ions / %',
                                                       redirect_output=True,
                                                       purge=False)
        # Position of electons in 3D
        self.electrons_position_graph       = new_plot(xpos=self.graph_x(2), ypos=self.graph_y(2),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0, xmax=self.ccp.length*1E3,
                                                       ymin=0, ymax=self.ccp.length*1E3,
                                                       zmin=0, zmax=self.ccp.distance*1E3,
                                                       style=self.plot_style,
                                                       plot_type='3D',
                                                       title='Electron location',
                                                       xlabel='x / mm',
                                                       ylabel='y / mm',
                                                       zlabel='z / mm',
                                                       redirect_output=True,
                                                       purge=False)
        # Position of ions in 3D
        self.ions_position_graph            = new_plot(xpos=self.graph_x(3), ypos=self.graph_y(2),
                                                       terminal=self.str_terminal,
                                                       width=self.window_width,
                                                       height=self.window_height,
                                                       xmin=0, xmax=self.ccp.length*1E3,
                                                       ymin=0, ymax=self.ccp.length*1E3,
                                                       zmin=0, zmax=self.ccp.distance*1E3,
                                                       style=self.plot_style,
                                                       plot_type='3D',
                                                       title='Ion location',
                                                       xlabel='x / mm',
                                                       ylabel='y / mm',
                                                       zlabel='z / mm',
                                                       redirect_output=True,
                                                       purge=False)
        

    def graph_list(self, expanded=False):
        """ Return a list of all the plots """
        
        return plot_list(expanded=expanded, printout=False, getstring=True)

    def print_graphs(self, terminal, dirname=None, options=None):
        """ Export all the graphs to files """

        return plot_print_all(terminal, dirname, options)

    def reset_graphs(self, ccp, parameters):
        """ Reset graphs dimensions and other parameters """

        self.ccp = ccp
        if parameters.dot_points: self.plot_style = 'dots'
        else:                     self.plot_style = 'points'        
        
        # Charge density vs z-position
        plot_set(self.charge_density_graph,         xmin=0, xmax=self.ccp.distance*1E3)
                 
        # Electric potential vs z-position
        plot_set(self.potential_graph,              xmin=0, xmax=self.ccp.distance*1E3)

        # Electron energy vs angle
        plot_set(self.electrons_phase_space_graph,  style=self.plot_style)
        
        # Ion energy vs angle
        plot_set(self.ions_phase_space_graph,       style=self.plot_style)       

        # Electron energy vs z-position
        plot_set(self.electrons_phase_space_graph2, xmin=0, xmax=self.ccp.distance*1E3, style=self.plot_style)
        
        # Ion energy vs z-position
        plot_set(self.ions_phase_space_graph2,      xmin=0, xmax=self.ccp.distance*1E3, style=self.plot_style)

        # Position of electons in 3D
        plot_set(self.electrons_position_graph,     xmin=0, xmax=self.ccp.length*1E3,
                                                    ymin=0, ymax=self.ccp.length*1E3,
                                                    zmin=0, zmax=self.ccp.distance*1E3,
                                                    style=self.plot_style)
        # Position of ions in 3D
        plot_set(self.ions_position_graph,          xmin=0, xmax=self.ccp.length*1E3,
                                                    ymin=0, ymax=self.ccp.length*1E3,
                                                    zmin=0, zmax=self.ccp.distance*1E3,
                                                    style=self.plot_style)
        

    def clear_graphs(self):
        """ Clear all the graph windows """
        
        #plot_clear_all()
        plot_reset_all()        

    def raise_graphs(self):
        """ Raises all graph windows """

        plot_raise_all()


    def lower_graphs(self):
        """ Lowers all graph windows """

        plot_lower_all()
                                                  
    def reset_history(self):
        """ Reset to zero all values of history """
        
        self.i_history = 0
        self.history_time.fill(0)
        self.history_e_average.fill(0)
        self.history_n_electrons.fill(0)
        self.history_nw_electrons.fill(0)
                
    def update_history(self, charges, time):
        """ Updates the history array """
        
        if (self.i_history >= self.history_limit):
            (self.history_time,new_i) = decimate(self.history_time,         self.history_decimation)
            self.history_e_average    = decimate(self.history_e_average,    self.history_decimation)[0]
            self.history_n_electrons  = decimate(self.history_n_electrons,  self.history_decimation)[0]
            self.history_nw_electrons = decimate(self.history_nw_electrons, self.history_decimation)[0]
            self.i_history = new_i

        self.history_time[self.i_history]         = time
        self.history_e_average[self.i_history]    = charges.e_average(0)
        self.history_n_electrons[self.i_history]  = charges.n_active(0)
        self.history_nw_electrons[self.i_history] = charges.n_active(0) * charges.weight[0]
        self.i_history += 1

                                                       
    def plot_positions(self, time, charges):
        """ Plots positions of electrons and ions in 3D """

        # Update time on graphs
        self.plot_time(self.electrons_position_graph, time)                
        self.plot_time(self.ions_position_graph, time)
        
        # Plot electrons graphs
        plot3d(self.electrons_position_graph,
               charges.x[0][charges.active[0]]*1.0E3,
               charges.y[0][charges.active[0]]*1.0E3, 
               charges.z[0][charges.active[0]]*1.0E3,
               label=''+charges.names[0])

        # Plot ions graphs
        self.plot_list = []
        for i in range(1, charges.types):
            #print "Plotting ion type:"+str(i)+" (" + charges.names[i] + ") active n. "+str(charges.active[i].sum())
            # Plot only if there are active ions
            if (charges.n_active(i) > 0):
                self.plot_list.append( [ charges.x[i][charges.active[i]]*1.0E3, 
                                         charges.y[i][charges.active[i]]*1.0E3,
                                         charges.z[i][charges.active[i]]*1.0E3,
                                         charges.names[i], None ] )                
        plot_curves(self.ions_position_graph,self.plot_list)          
                

    def plot_phase_space_data(self, time, charges):
        """ Plots electron and ion energies vz angle and z-position """

        # Update time on graphs
        self.plot_time(self.electrons_phase_space_graph ,time)
        self.plot_time(self.electrons_phase_space_graph2 ,time)
        self.plot_time(self.ions_phase_space_graph, time)
        self.plot_time(self.ions_phase_space_graph2, time)
        
        # Plot electron energies vs angle and mean energy vs mean angle
        self.plot_list = []
        self.plot_list.append( [ numpy.degrees(charges.theta[0][charges.active[0]]),
                                 charges.energies(0),
                                 charges.names[0] + ' ', None ] )
        self.plot_list.append( [ numpy.array([numpy.degrees(charges.theta_average(0))]),
                                 numpy.array([charges.e_average(0)]),
                                 '<' + charges.names[0] + '>', None ] )
        plot_curves(self.electrons_phase_space_graph, self.plot_list)
        
        # Plot electron energies vs z-position
        plot2d(self.electrons_phase_space_graph2,
               charges.z[0][charges.active[0]]*1.0E3,
               charges.energies(0),
               label=charges.names[0])
        
        # Plot ion energy vs angle
        self.plot_list = []
        for i in range(1,charges.types):
            # Plot only if there are active ions
            if (charges.n_active(i) > 0):              
                self.plot_list.append( [ numpy.degrees(charges.theta[i][charges.active[i]]),
                                         charges.energies(i),
                                         charges.names[i], None ] )
        plot_curves(self.ions_phase_space_graph, self.plot_list)
            
        # Plot ion energy vs z-position
        self.plot_list = []
        for i in range(1,charges.types):
            # Plot only if there are active ions
            if (charges.n_active(i) > 0):              
                self.plot_list.append( [ charges.z[i][charges.active[i]]*1.0E3,
                                         charges.energies(i),
                                         charges.names[i], None ] )
        plot_curves(self.ions_phase_space_graph2, self.plot_list)

                                        
    def plot_distributions(self, time, charges):
        """ Plots electron and ion energy distributions """
        
        self.plot_time(self.eedf_graph, time)
        self.plot_time(self.iedf_graph, time)

        # Calculate the number of bins for electron distribution, as the square root of the number of electrons
        self.nbins_ele = int( math.sqrt( charges.n_active(0) ) )
        if (self.nbins_ele == 0): self.nbins_ele = 1
#        if (self.nbins_ele > N_BINS_LOG): plot_set(self.eedf_graph, logx=True)
#        else:                             plot_set(self.eedf_graph, logx=False) 

        # Plot EEDF
        if (charges.n_active(0) > 0):
            # Calculate bins and frequencies
            (h,b) = numpy.histogram( charges.energies(0), bins=self.nbins_ele,
                                    range=(0.0, charges.e_max(0)), density=False )
            # Remove the first item of the bins array (left value of the first bin)
            b = b[1:]
            # Calculate the bins width
            d = ( b.max() - b.min() ) / len(b)
            # Create an array with the center point of each bin
            b = b - d / 2.0
            # Transform from absolute frequencies (number of particles in each bin) to relative frequencies
            # (percentage of particles in each bin)
            h = h*100.0 / ( charges.n_active(0) * 1.0 )
            plot2d(self.eedf_graph, b, h, label=charges.names[0])
        
        # Calculate the number of bins for ion distribution, as the square root of the minimum number of ions
        # (all ion types)
        self.n = charges.n
        for i in range(1, charges.types):
            if (charges.n_active(i) < self.n): self.n = charges.n_active(i)
        self.nbins_ion = int(math.sqrt(self.n))
        if (self.nbins_ion == 0): self.nbins_ion = 1
#        if (self.nbins_ion > N_BINS_LOG): plot_set(self.iedf_graph, logx=True)
#        else:                             plot_set(self.iedf_graph, logx=False)
        
        # Set x-axis limit for ion graphs to max ion energy (all ion types)
        self.ion_emax = 0.0
        for i in range(1,charges.types):
            if (charges.n_active(i) > 0):
                if (charges.e_max(i) > self.ion_emax): self.ion_emax = charges.e_max(i)
        if (self.ion_emax > 0.0): plot_set(self.iedf_graph, xmax=self.ion_emax)

        # Plot IEDF for each ion type
        self.plot_list = []
        for i in range(1, charges.types):
            if (charges.n_active(i) > 0):
                # Calculate bins and frequencies
                (h,b) = numpy.histogram( charges.energies(i), bins=self.nbins_ion,
                                         range=(0.0, self.ion_emax), density=False )
                # Remove the first item of the bins array (left value of the first bin)
                b = b[1:]
                # Calculate the bins width
                d = ( b.max() - b.min() ) / len(b)
                # Create an array with the center point of each bin
                b = b - d / 2.0
                # Transform from absolute frequencies (number of particles in each bin) to relative frequencies
                # (percentage of particles in each bin)
                h = h*100.0 / ( charges.n_active(i) * 1.0 )
                self.plot_list.append( [b, h, charges.names[i], None] )
        plot_curves(self.iedf_graph, self.plot_list)
                    
                                        
    def plot_n_and_energy_data(self, time):
        """ Plot mean electron enery and number of electrons vs simulation time  """
        
        self.plot_time(self.electrons_average_energy_graph, time)
        self.plot_time(self.n_electrons_graph, time)

        # Rescale x-scale to start from min value left after decimation
        self.xmin = self.history_time[0]*1.0E9
        self.xmax = self.history_time[self.i_history-1]*1.0E9

        # Plot mean energy vs time
        plot_set(self.electrons_average_energy_graph, xmin=self.xmin, xmax=self.xmax)
        plot2d(self.electrons_average_energy_graph, self.history_time*1.0E9, self.history_e_average) 
                                                       
        #Plot number of computational and real (weighted) electrons                                                       
        plot_set(self.n_electrons_graph, xmin=self.xmin, xmax=self.xmax)
        self.plot_list = [ [self.history_time*1.0E9, self.history_n_electrons,  'comp e-', None],
                           [self.history_time*1.0E9, self.history_nw_electrons, 'real e-', None] ]
        plot_curves(self.n_electrons_graph,self.plot_list) 

        

    def plot_potential(self, time, ccp):
        """ Plots charge densiy and elecgric potential """

        self.plot_time(self.charge_density_graph, time)
        self.plot_time(self.potential_graph, time)
        plot2d(self.potential_graph,      ccp.grid_points*1.0E3, ccp.potential)
        plot2d(self.charge_density_graph, ccp.grid_points*1.0E3, ccp.charge_density)


#+----------------------+
#| Cross sections plots |
#+----------------------+

def plot_cross_sections(neutrals, electrons=True, ions=True, recombination=False, plot_all=False, dot_points=False):
    """ Plots cross sections """

    (status,message) = (0, OK)
    
    if electrons:
        (status, message) = neutrals.plot_xsec_electrons( plot_single      = not plot_all,
                                                          plot_total       = not plot_all, 
                                                          plot_frequencies = plot_all, 
                                                          plot_relative    = plot_all, 
                                                          plot_boundaries  = plot_all,
                                                          dot_points       = dot_points )
        if (status != 0): return (status, message)
                
    if ions:
        for ion_type in range(neutrals.types):
            (status, message) = neutrals.plot_xsec_ions( ion_type         = ion_type, 
                                                         plot_total       = not plot_all, 
                                                         plot_frequencies = plot_all, 
                                                         plot_relative    = plot_all, 
                                                         dot_points       = dot_points )
            if (status != 0): return (status, message)

    if recombination:                 
            (status, message) = neutrals.plot_xsec_ele_ion_recomb( plot_total = True, 
                                                                   plot_rates = plot_all,
                                                                   dot_points = dot_points )
   
    return (status, message)
