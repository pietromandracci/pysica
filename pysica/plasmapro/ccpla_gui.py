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

""" Simulation of a capacitively coupled plasma discharge: GUI """

# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

# Modules from the standard Python library
import math
from subprocess import run
from multiprocessing import Process, Pipe
import tkinter as tk
import tkinter.messagebox
import tkinter.scrolledtext

# Import required constants and parameters
from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *

# Import required modules, classes, and functions
from pysica.plasmapro.ccpla_init import *
from pysica.plasmapro.ccpla_print import *
from pysica.plasmapro.ccpla_plot import *
from pysica.plasmapro.ccpla_kernel import kernel
                
class CcplaWindow(tk.Frame):
    """ Create the main window """
        
    def __init__(self, charges, neutrals, ccp, parameters, options,
                 screen_width=SCREEEN_WIDTH, screen_height=SCREEEN_HEIGHT, master=None):

        # Class instances
        self.charges       = charges
        self.neutrals      = neutrals
        self.ccp           = ccp
        self.parameters    = parameters
        self.options       = options
        self.debug         = (self.options.debug_lev > 0)        

        # Create an ensamble of runtime plots and initialize them
        self.screen_width  = screen_width
        self.screen_height = screen_height
        self.runtime_plots = CcplaPlots(self.parameters, self.charges,
                                        self.ccp, self.screen_width,
                                        self.screen_height)

        # Normal variables
        self.sim_status         = False #True=running; False=stopped; None=paused
        self.error              = False
        #self.plot_counter       = 0
        #self.save_counter       = 0

        # Tk control variables
        self.dt_exp             = tk.IntVar()
        self.plot_delay         = tk.IntVar()
        self.plot_history       = tk.BooleanVar()
        self.plot_phase_space   = tk.BooleanVar()
        self.plot_field         = tk.BooleanVar()
        self.plot_distributions = tk.BooleanVar()
        self.plot_3D_positions  = tk.BooleanVar()

        # Activate all the plots by default
        self.set_all_plots(True)

        # Initialize Tk control variables 
        self.dt_exp.set(int(math.log10(self.parameters.dt_output)))
        self.plot_delay.set(DEF_PLOT_DELAY)

        # Create main window
        tk.Frame.__init__(self, master)

        # Add menu item "File"
        self.menubutton_file           = tk.Menubutton(self, text='File', justify=tk.LEFT)
        self.menubutton_file.grid(row=0, column=0, columnspan=1, sticky=tk.NW)

        self.menubutton_file.menu_file = tk.Menu(self.menubutton_file, tearoff=0)
        self.menubutton_file['menu']   = self.menubutton_file.menu_file
        #self.menubutton_file.menu_file.add_command(label='Change save directory name',
        #                                           command=self.change_basename)
        self.menubutton_file.menu_file.add_command(label='Reload configuration files',
                                                   command=self.reload_parameters)
        self.menubutton_file.menu_file.add_command(label='Edit configuration files',
                                                   command=self.edit_parameters)
        self.menubutton_file.menu_file.add_command(label='Quit' , command=self.exit_gui)

        # Add menu item "Parameters"
        self.menubutton_parameters = tk.Menubutton(self, text="Parameters", justify=tk.LEFT)
        self.menubutton_parameters.grid(row=0, column=1, columnspan=4, sticky=tk.NW)

        self.menubutton_parameters.menu_parameters = tk.Menu(self.menubutton_parameters,
                                                             tearoff=0)
        self.menubutton_parameters['menu'] = self.menubutton_parameters.menu_parameters
        self.menubutton_parameters.menu_parameters.add_command(
                label='Show physical parameters',
                command=lambda: self.show_parameters(physical=True))
        self.menubutton_parameters.menu_parameters.add_command(label='Show simulation parameters',
                                                               command=lambda: self.show_parameters(simulation=True))
        self.menubutton_parameters.menu_parameters.add_command(label='Show output parameters',
                                                               command=lambda: self.show_parameters(output=True))                
        self.menubutton_parameters.menu_parameters.add_command(label='Show gas properties', command=self.show_gases)
        self.menubutton_parameters.menu_parameters.add_separator()
        self.menubutton_parameters.menu_parameters.add_command(label='Show e-/neutral impact cross sections',
                                                               command=lambda: self.draw_cross_sections(electrons=True))
        self.menubutton_parameters.menu_parameters.add_command(label='Show ion/neutral impact cross sections',
                                                               command=lambda: self.draw_cross_sections(ions=True))  
        self.menubutton_parameters.menu_parameters.add_command(label='Show e-/ion recomb cross sections',
                                                               command=lambda: self.draw_cross_sections(recombination=True))
        self.menubutton_parameters.menu_parameters.add_separator()                
        self.menubutton_parameters.menu_parameters.add_command(label='Show e-/neutral impact parameters',
                                                               command=lambda: self.draw_cross_sections(electrons=True,
                                                                                                        plot_all=True))  
        self.menubutton_parameters.menu_parameters.add_command(label='Show ion/neutral impact parameters',
                                                               command=lambda: self.draw_cross_sections(ions=True, plot_all=True))
        self.menubutton_parameters.menu_parameters.add_command(label='Show e-/ion recombination parameters',
                                                               command=lambda: self.draw_cross_sections(recombination=True,
                                                                                                        plot_all=True))

        # Add menu item "Runtime plots"
        self.menubutton_runtime_plot = tk.Menubutton(self, text='Runtime Plots', justify=tk.LEFT)
        self.menubutton_runtime_plot.grid(row=0, column=5, columnspan=5, sticky=tk.NW)

        self.menubutton_runtime_plot.menu_runtime_plot = tk.Menu(self.menubutton_runtime_plot, tearoff=0)
        self.menubutton_runtime_plot['menu']   = self.menubutton_runtime_plot.menu_runtime_plot
        self.menubutton_runtime_plot.menu_runtime_plot.add_command(label='Select all',
                                                                    command = lambda: self.set_all_plots(True))
        self.menubutton_runtime_plot.menu_runtime_plot.add_command(label='Unselect all',
                                                                    command = lambda: self.set_all_plots(False))
        self.menubutton_runtime_plot.menu_runtime_plot.add_separator()
        self.menubutton_runtime_plot.menu_runtime_plot.add_checkbutton(label='Mean el energy and number vs time',
                                                                       variable=self.plot_history)
        self.menubutton_runtime_plot.menu_runtime_plot.add_checkbutton(label='Phase space plots',
                                                                       variable=self.plot_phase_space)
        self.menubutton_runtime_plot.menu_runtime_plot.add_checkbutton(label='Eletric potential and charge',
                                                                       variable=self.plot_field)
        self.menubutton_runtime_plot.menu_runtime_plot.add_checkbutton(label='EEDF and IEDF',
                                                                       variable=self.plot_distributions)
        self.menubutton_runtime_plot.menu_runtime_plot.add_checkbutton(label='3D e- and ion positions',
                                                                       variable=self.plot_3D_positions)
        self.menubutton_runtime_plot.menu_runtime_plot.add_separator()        
        self.menubutton_runtime_plot.menu_runtime_plot.add_command(label='Plot list', command=self.show_plot_info)
        self.menubutton_runtime_plot.menu_runtime_plot.add_command(label='Plot export', command=self.print_plots)        

        # Add menu item "Help"
        self.menubutton_help           = tk.Menubutton(self, text='Help', justify=tk.LEFT)
        self.menubutton_help.grid(row=0, column=10, columnspan=1, sticky=tk.NW)
        self.menubutton_help.menu_help = tk.Menu(self.menubutton_help, tearoff=0)
        self.menubutton_help['menu']   = self.menubutton_help.menu_help
        self.menubutton_help.menu_help.add_command(label='Online documentation (open in browser)', command=self.show_doc)
        self.menubutton_help.menu_help.add_command(label='About',         command=self.show_about)

        # Add text region
        self.text_output = tk.Text(self,
                                   height=options.text_window_height,
                                   width=options.text_window_width,
                                   font=options.text_window_font,
                                   background=TEXT_WINDOW_BACK,
                                   foreground=TEXT_WINDOW_FORE)
        self.text_output.grid(row=1, column=0, columnspan=50)       

        # Add dt output scale
        self.dt_out_set = tk.Scale(self, from_=MIN_DT_OUTPUT_EXP, to=MAX_DT_OUTPUT_EXP, showvalue=False,
                                   orient=tk.HORIZONTAL, length=300, 
                                   variable=self.dt_exp, command=self.set_dt_output)
        self.dt_out_set.grid(row=2, column=5, columnspan=10)

        # Add plot delay scale
        self.dt_out_set = tk.Scale(self, from_=MIN_PLOT_DELAY, to=MAX_PLOT_DELAY, resolution=RES_PLOT_DELAY, showvalue=False,
                                   orient=tk.HORIZONTAL, length=300, 
                                   variable=self.plot_delay, command=self.set_plot_delay)
        self.dt_out_set.grid(row=3, column=5, columnspan=10)


        # Add label to show dt_output
        self.dt_out_label = tk.Label(self, font = ('mono', 8 ), width=22, justify=tk.LEFT, anchor=tk.W)
        self.show_dt_out_label()

        # Add label to show paused/running
        self.status_label = tk.Label(self)
        self.show_status_label(text='Press RESET to start', color='red')

        # Add buttons
        self.button_reset = tk.Button(self, text='RESET', foreground='red', background='light yellow',
                                      command=self.reset_simulation)
        self.button_reset.grid(row=2, column=0, columnspan=2, sticky=tk.W)
        self.button_start = tk.Button(self, text='START', foreground='red', command=self.start_simulation)
        self.button_start.grid(row=2, column=35, columnspan=5)
        self.button_pause = tk.Button(self, text='Pause', foreground='blue', command=self.pause_restore_simulation)
        self.button_pause.grid(row=2, column=40, columnspan=5)
        self.button_stop = tk.Button(self, text='STOP', foreground='red', command=self.stop_simulation)
        self.button_stop.grid(row=2, column=45, columnspan=5)

        # Set status of buttons at program start
        self.button_reset["state"] = tk.NORMAL
        self.button_start["state"] = tk.DISABLED
        self.button_pause["state"] = tk.DISABLED
        self.button_stop["state"]  = tk.DISABLED                

        self.set_menu()

        self.grid()

        self.after(IDLE_TIME_START, self.show_about)

        self.initialize_kernel()

        
    # +-----------------------------------------------------+    
    # | Method to quit the main window and exit the program |
    # +-----------------------------------------------------+
    
    def exit_gui(self):
        """ Terminate the kernel and close the main window """
        self.quit_kernel()
        self.f_kernel.join()
        self.master.destroy() #self.quit()

        
    # +--------------------------------------+
    # | Methods to set the window properties |
    # +--------------------------------------+  
        
    def set_menu(self, change_names=True, edit_parameters=True, reload_parameters=True):
        index = self.menubutton_parameters.menu_parameters.index(tk.END)
        if self.parameters.isactive_recomb: self.state = tk.NORMAL
        else:                               self.state = tk.DISABLED
        self.menubutton_parameters.menu_parameters.entryconfigure(index,   state=self.state)
        self.menubutton_parameters.menu_parameters.entryconfigure(index-4, state=self.state)
        if change_names:      self.state = tk.NORMAL
        else:                 self.state = tk.DISABLED
        self.menubutton_file.menu_file.entryconfigure(0, state=self.state)         
        if edit_parameters:   self.state = tk.NORMAL
        else:                 self.state = tk.DISABLED
        self.menubutton_file.menu_file.entryconfigure(1, state=self.state)
        if reload_parameters: self.state = tk.NORMAL
        else:                 self.state = tk.DISABLED
        self.menubutton_file.menu_file.entryconfigure(2, state=self.state)
        

    # +-------------------------------------------+
    # | Methods to interact with kernel and pipes | 
    # +-------------------------------------------+    

    def initialize_kernel(self):
        """ Initialize the kernel process and the pipes """

        if self.debug: print('[GUI]   -> INITIALIZING KERNEL')
        
        (self.parent_signal, self.child_signal) = Pipe(duplex=True)
        (self.parent_time,   self.child_time)   = Pipe(duplex=True)        
        (self.parent_data,   self.child_data)   = Pipe(duplex=True)        
        self.f_kernel = Process(target=kernel,
                                args=(self.child_signal,
                                      self.child_time,
                                      self.child_data,
                                      (self.options.debug_lev >= KERNEL_DEBUG_LEV)
                                     )
                               )
        self.f_kernel.daemon = True
        self.f_kernel.start()
        
               
    def start_kernel(self):
        """ Tell the kernel to start iteration and send initializerd data to it """
        
        if self.debug: print('[GUI]   -> Sending start signal to kernel')
        self.parent_signal.send('start')
        if self.debug: print('[GUI]   -> Sending initialized data to kernel')
        self.parent_data.send((self.charges, self.neutrals, self.ccp, self.parameters, self.options))
        if self.debug: print('[GUI]   -> Sending requested simulation time to kernel')
        self.parent_time.send(self.parameters.dt_output)

    def stop_kernel(self):
        """ Flush pipes and tell the kernel to stop iteration """
        
        self.flush_pipes()
        self.parent_signal.send('stop')
        
    def quit_kernel(self):
        """ Tell kernel to quit, if doesn't, kill it """
        
        if self.debug: print('[GUI]   -> quitting kernel')
        self.show_status_label(text='Waiting for kernel to shutdown', color='red')
        self.parent_signal.send('quit')        
        self.after(IDLE_TIME_SHUTDOWN)
        if self.f_kernel.is_alive():
            self.show_status_label(text='Forcing kernel to shut down...', color='red')
            if self.debug: print('[GUI]   -> terminating kernel')            
            self.f_kernel.terminate()
            self.after(IDLE_TIME_SHUTDOWN)
            if self.f_kernel.is_alive():
                if self.debug: print('[GUI]   -> killing kernel')
                self.f_kernel.kill()
        if self.debug: print('[GUI]   -> joining kernel')
        self.f_kernel.join()        


    def flush_pipes(self, flush_signal=True, flush_time=True, flush_data=True):
        """ Read any residual data from the pipes """
        
        if self.debug: print('[GUI]   -> CHECKING PIPES')
        if flush_signal:
            while self.parent_signal.poll():
                if self.debug: print('[GUI]   -> FLUSHING SIGNAL PIPE')
                signal = self.parent_signal.recv()
        if flush_time:    
            while self.parent_time.poll():
                if self.debug: print('[GUI]   -> FLUSHING TIME PIPE')                
                time = self.parent_time.recv()
        if flush_data:
            while self.parent_data.poll():
                if self.debug: print('[GUI]   -> FLUSHING DATA PIPE')                
                (self.charges, self.neutrals, self.ccp, self.parameters, self.options) = self.parent_data.recv()

        
    def close_pipes(self):
        """ Close the pipes """
        if self.debug: print('[GUI]   -> CLOSING PIPES')
        
        self.parent_signal.close()
        self.parent_time.close()
        self.parent_data.close()

        
    # +------------------------------------------------+
    # | Main iteration and methods to interact with it |
    # +------------------------------------------------+    
        
    def iteration(self):
        """ Main simulation iteration """

        # If self.sim_status if False, the simulation must end: exit from simulation loop
        if (self.sim_status is False):
            return        
        # If self.sim_status if True or None, the simulation is running or paused
        else:            
            # Check if the kernel has new calculated data
            self.data_ready = self.parent_signal.poll()
            if self.data_ready:
                self.show_status_label(text='Kernel is ready  ', color='green')
                if (self.sim_status is None): self.button_stop["text"] = 'STOP'          
            else:
                self.show_status_label(text='Kernel is busy ...', color='red')
                if (self.sim_status is None): self.button_stop["text"] = 'KILL'
            self.update_idletasks()
            # If there is new data and simulation is not paused
            if (self.data_ready and (self.sim_status is True)):
                message = self.parent_signal.recv()
                if self.debug: print('[GUI]   -> Message from kernel: ' + message)# + '\n')
                if self.debug: print('[GUI]   -> Reading data from kernel')
                (self.charges, self.neutrals, self.ccp, self.parameters, self.options) = self.parent_data.recv()
                if self.debug: print('[GUI]   -> Data read from kernel')                
                # Syncronize requested iteration duration with window selector
                self.parameters.dt_output = 10.0**int(self.dt_exp.get())                
                # Update data on text window
                self.show_runtime_info()                
                # Update plots 
                self.draw_plots()                
                # Save data to files
                if (self.parameters.save_delay > 0):
                    self.save_counter += 1
                    if (self.save_counter >= self.parameters.save_delay):
                        if self.debug: print('[GUI]   -> Saving data')                                        
                        self.charges.save_data_to_files()
                        self.neutrals.save_data_to_files(self.charges.time)
                        self.ccp.save_data_to_files(self.charges.time)
                        if self.debug: print('[GUI]   -> Data saved')
                        self.save_counter = 0                    
                # Save actual time to be shown in next iteration
                self.time_before          = self.charges.time
                self.n_active_el_before   = self.charges.n_active(0)
                self.electric_bias_before = self.ccp.V
                if (self.charges.n_active(0) <= 0):
                    self.show_end_simulation(text ='\n\nSimulation interrupted \n (no more electrons)\n\n',
                                             color='red')
                    self.stop_simulation(confirm_stop=False)
                elif ( (self.parameters.sim_duration > 0) and (self.charges.time >= self.parameters.sim_duration) ):
                    self.show_end_simulation(text='\n\nSimulation completed\n\n', color='blue')
                    self.stop_simulation(confirm_stop=False)
                else:
                    self.parent_signal.send('continue')
                    self.parent_time.send(self.parameters.dt_output)                    
        # Show the simulation output time
        self.show_dt_out_label()       
        # Before repeating the loop, update widgets and wait some time to allow responsivity of the window
        self.update_idletasks()
        self.after(IDLE_TIME, self.iteration)
        
                
    def start_simulation(self):
        """ Begin a new simulation run """

        self.plot_counter = 0
        self.time_before          = self.charges.time
        self.n_active_el_before   = self.charges.n_active(0)
        self.electric_bias_before = self.ccp.V      
        
        self.draw_plots(force_plot=True)       

        self.set_menu(change_names=False, edit_parameters=False, reload_parameters=False)
        self.button_start["state"] = tk.DISABLED
        self.button_pause["state"] = tk.NORMAL
        self.button_stop["state"]  = tk.DISABLED         
        self.show_runtime_info()
        
        self.update_idletasks()
         
        self.sim_status = True
        if self.debug: print('[GUI]   -> start kernel')
        self.start_kernel()
        if self.debug: print('[GUI]   -> start iteration')
        self.after(IDLE_TIME, self.iteration)

        
    def pause_restore_simulation(self):
        """ Pause the simulation or restore it from paused state """
        
        if (self.sim_status == True):
            self.sim_status = None
            self.button_stop["state"] = tk.NORMAL
            self.button_pause["text"] = 'Continue'
        elif (self.sim_status == None):
            self.sim_status = True
            self.button_stop["state"] = tk.DISABLED
            self.button_pause["text"] = 'Pause'
        if self.debug: print('[GUI]   -> sim_status = ' + str(self.sim_status))
                

    def stop_simulation(self, confirm_stop=True):
        """ End the running simulation """
        
        if self.data_ready:
            if confirm_stop:
                self.answer = tkinter.messagebox.askyesno( title='Stop simulation',
                                                           message=('Stop the simulation ?'),
                                                           icon='question')
                if not self.answer: return
            self.stop_kernel()            
        else:
            self.answer = tkinter.messagebox.askyesno( title='WARNING',
                                                       message=('WARNING! Killing the running kernel may produce data loss \n'
                                                                + 'Do you want to proceed ?'),
                                                       icon='warning')
            if not self.answer: return
            else:
                self.quit_kernel()
                self.close_pipes()
                self.initialize_kernel()
        #self.flush_pipes()
        self.sim_status = False        
        self.button_pause["state"]  = tk.DISABLED
        self.button_stop["text"]    = 'STOP'
        self.button_stop["state"]   = tk.DISABLED
        self.show_status_label(text='Press RESET to start over', color='blue')
        self.button_reset["state"]  = tk.NORMAL
        self.set_menu()       


    def reset_simulation(self):
        """ Reset all data for a new simulation run """
        
        # Initialize particle data
        if self.debug: print('[GUI]   -> Initializing particle data')        
        (status, message) = initialize_ensambles(self.charges, self.neutrals, self.parameters, self.options)
        if (status != 0):
            self.show_message(message, fore='red', back='light yellow')
            return
        self.parameters.dt_output = 10.0**int(self.dt_exp.get())

        # Initialize plots
        if self.debug: print('[GUI]   -> Initializing plots')        
        #self.runtime_plots.initialize_graphs(reset=True)
        self.runtime_plots.reset_graphs(self.ccp, self.parameters)
        self.runtime_plots.clear_graphs()
        self.runtime_plots.reset_history()
        self.draw_plots(force_plot=True)                    
        self.after(IDLE_TIME_GRAPHS)        
        self.runtime_plots.lower_graphs()
        self.plot_counter = 0

        # Create save directory and initialize savefiles
        if (self.parameters.save_delay > 0):
            if self.debug: print('[GUI]   -> Initializing save directory')
            initialize_filenames(self.neutrals, self.parameters)                        
            self.answer = tkinter.messagebox.askyesno( title='Check',
                                                       message=('Use this name for the save directory? \"'
                                                                + str(self.parameters.basename)
                                                                + '\"') )
            if not self.answer: self.change_basename()
            create_save_dir(self.parameters)
            if self.debug: print('[GUI]   -> Initializing savefiles')
            self.charges.initialize_savefiles(self.parameters.filename_stat_ele,
                                              self.parameters.filename_distrib_ele,
                                              self.parameters.filename_distrib_ion,
                                              self.parameters.filename_epos_z,
                                              self.parameters.filename_ipos_z,                                              
                                              append=False, sep=SEP, ext=EXT)
            self.neutrals.initialize_savefile(self.parameters.filename_stat_neu,
                                              append=False, sep=SEP, ext=EXT)
            self.ccp.initialize_savefiles(self.parameters.filename_I, self.parameters.filename_V,
                                          append=False, sep=SEP, ext=EXT)            
            self.save_counter = 0          
            
        # Arrange flags, windows and buttons status
        if self.debug: print('[GUI]   -> Initializing status')  
        self.sim_status = False
        self.show_runtime_info(erase_only=True)
        self.show_status_label(text='Press START to begin simulation', color='red')
        self.button_reset["state"] = tk.DISABLED
        self.button_pause["text"]  = 'Pause'
        self.button_pause["state"] = tk.DISABLED
        self.button_stop["text"]   = 'STOP'
        self.button_stop["state"]  = tk.DISABLED
        self.button_start["state"] = tk.NORMAL

        self.update_idletasks()
        


    # +-----------------------------------------------+
    # | Metods to interact with simulation parameters |
    # +-----------------------------------------------+

    def edit_parameters(self):
        try:
            #editor_process = run([EDITOR_NAME, FILENAME_CONFIG, FILENAME_NEUTRALS])
            editor_process = run([EDITOR_NAME, FILENAME_CONFIG])
        except FileNotFoundError:
            string = ('Error while trying to execute command: \"'
                      + EDITOR_NAME
                      + ' ' 
                      + FILENAME_CONFIG
#                      + ' '
#                      + FILENAME_NEUTRALS
                      + '\" ')
            self.show_message(message=string, fore='red', back='light yellow')
            return
        self.reload_parameters()
        
                        
    def reload_parameters(self): 
        self.error   = False
        self.warning = False
        ERRORS = ''
        # Reload parameters from configuration file
        (status, message) = initialize_parameters(self.parameters,
                                                  verbose=False,
                                                  saveonly=False,
                                                  filename_config=FILENAME_CONFIG,
                                                  filename_defaults=FILENAME_DEFAULTS,
                                                  restricted=True)
        if (status != 0):
            self.error= True
            ERRORS += message + '\n' 
        # Recofigure reactor parameters
        self.ccp.__init__(self.parameters.distance,
                          self.parameters.length,
                          self.parameters.V_bias,
                          self.parameters.frequency,
                          self.parameters.phase,
                          self.parameters.N_cells,
                          self.parameters.lateral_loss)
        
        # Reconfigure neutral particle properties
#        self.neutrals.__init__(self.parameters.N_sigma,
#                               self.parameters.N_sigma_ions,
#                               self.parameters.T_neutrals,
#                               self.parameters.p_neutrals,
#                               self.parameters.min_scattered,
#                               self.parameters.isactive_recomb,
#                               filename=FILENAME_NEUTRALS)
#        (status, message) = self.neutrals.read_error
#        if (status !=0):
#            self.error= True
#            ERRORS += 'Error in file \"' + FILENAME_NEUTRALS + '\": ' + message + '\n'
#        # Read neutral properties from file
#        (status, message) = self.neutrals.read_properties(FILENAME_NEUTRALS,
#                                                          SEP, 
#                                                          self.parameters.e_min_sigma,
#                                                          self.parameters.e_max_sigma,
#                                                          self.parameters.e_min_sigma_ions,
#                                                          self.parameters.e_max_sigma_ions)
#        if (status !=0):
#            self.error= True
#            ERRORS += 'Error in file \"' + FILENAME_NEUTRALS + '\": ' + message + '\n'

#        if not self.error:
#            # Reload cross sections and calculate other impact parameters
#            (status, message) = initialize_cross_sections(self.neutrals, self.options, self.parameters.isactive_recomb)
#            if (status !=0):
#                self.error= True
#                ERRORS += message + '\n'
        if self.error:
            self.button_reset["state"] = tk.DISABLED
            self.button_start["state"] = tk.DISABLED
            self.show_message('Errors encountered while loading parameters \n\n' + ERRORS,
                              fore='red', back='light yellow')
            self.show_status_label(text='Errors: check configuration files', color='red')
        else:
            # Reconfigure electron and ion properties
            #self.charges.__init__(self.neutrals.types+1, self.parameters.Nmax_particles, START_WEIGHT,
            #                      self.parameters.rescale_factor)            
            initialize_ensambles(self.charges, self.neutrals, self.parameters, self.options)            
            self.show_message('Parameters reloaded successfully', fore='blue')            
            self.button_reset["state"] = tk.NORMAL
            self.button_start["state"] = tk.DISABLED
            self.set_menu()
            self.set_dt_output()
            self.show_status_label(text='Press RESET to start over', color='blue')
                        
    def set_dt_output(self, set_value=None):
        # If set_value is None, self.parameters.dt_output is assumed to have been already set
        if (set_value is not None): self.parameters.dt_output = 10.0**int(set_value)
        # Check that output periodicity is not smaller than timestep
        if ( (not self.parameters.dt_var) and (self.parameters.dt_output < self.charges.dt) ):
            self.parameters.dt_output = self.charges.dt
        self.dt_exp.set(int(math.log10(self.parameters.dt_output)))
        self.show_dt_out_label()

    def set_plot_delay(self, set_value):
        self.plot_counter = 0
        self.show_dt_out_label()                

    def set_all_plots(self, state=True):
        self.plot_history.set(state)
        self.plot_phase_space.set(state)
        self.plot_field.set(state)
        self.plot_distributions.set(state)
        self.plot_3D_positions.set(state)


    # +--------------------------------------+
    # | Methods to deal with data file names |
    # +--------------------------------------+  
        
    def change_basename(self):
        self.filenames_window = tk.Toplevel(self)
        self.filenames_window.title('Enter save directory name')
        self.filenames_window.resizable(False,False)
        self.filenames_window.message = tk.Entry(self.filenames_window, 
                                                 font=FONT_MESSAGE, foreground='blue',takefocus=True,
                                                 width=40,
                                                 justify=tk.LEFT, relief=tk.RIDGE)
        self.filenames_window.message.insert('0', self.parameters.basename)
        self.filenames_window.message.grid()
        self.filenames_window.button_ok = tk.Button(self.filenames_window, text="Ok", command=self.get_basename)
        self.filenames_window.button_ok.grid()
        self.filenames_window.message.focus_set()
        self.filenames_window.grab_set()
        self.wait_window(self.filenames_window)

    def get_basename(self):
        self.parameters.basename = self.filenames_window.message.get()
        initialize_filenames(self.neutrals, self.parameters, self.parameters.basename)
        self.filenames_window.destroy()

        
    # +------------------------------------------------+
    # | Methods to show the information on the windows |
    # +------------------------------------------------+        
    
    def show_dt_out_label(self):
        if (self.charges.dt == 0): self.dt_print = 'variable'
        else:                      self.dt_print = unit_manager.print_unit(self.charges.dt,'s', 1)
        # This is necessay since scale widget with resolution > 1  and minimum 1 will give zero (not 1) as minimum value
        if (self.plot_delay.get() == 0): self.dt_plot = self.parameters.dt_output
        else:                            self.dt_plot = self.plot_delay.get() * self.parameters.dt_output                
        self.dt_out_label["text"] = 'Output every: ' + unit_manager.print_unit(self.parameters.dt_output,'s', 1) + '\n' +\
                                    'Plot every:   ' + unit_manager.print_unit(self.dt_plot, 's', 1) + '\n' +\
                                    'Timestep:     ' + self.dt_print
        self.dt_out_label.grid(row=2, column=15, columnspan=8) 
              
    def show_status_label(self, text, color):
        self.status_label["text"]       = text
        self.status_label["foreground"] = color
        self.status_label.grid(row=2, column=22, columnspan=12)

    def show_doc(self):
        try:
            editor_process = run([BROWSER_NAME, DOC_URL])
        except FileNotFoundError:
            string = ('Error while trying to execute command: \"'
                      + BROWSER_NAME
                      + ' '
                      + DOC_URL
                      + '\"')
            self.show_message(message=string, fore='red', back='light yellow')
        
    def show_about(self):
        self.about_window = tk.Toplevel(self)
        self.about_window.title('Copyright Info')
        self.about_window.resizable(False,False)
        self.about_window.message = tk.Message(self.about_window, text=GPL_MESSAGE, justify=tk.CENTER, relief=tk.RIDGE)
        self.about_window.message.grid()
        self.about_window.button  = tk.Button(self.about_window, text="Dismiss", command=self.about_window.destroy)
        self.about_window.button.grid()
        self.about_window.attributes("-topmost", True)         
                        
    def show_gases(self):
        self.gases_window = tk.Toplevel(self)
        self.gases_window.title('Gas Properties')
        self.gases_window.resizable(False,False)
        self.gases_window.message = tk.Message(self.gases_window, text=print_gas_information(self.neutrals),
                                               font=FONT_MESSAGE_SMALL, foreground='blue', width=1500,
                                               justify=tk.LEFT, relief=tk.RIDGE)
        self.gases_window.message.grid()
        self.gases_window.button  = tk.Button(self.gases_window, text="Dismiss", command=self.gases_window.destroy)
        self.gases_window.button.grid()
                
    def show_parameters(self, physical=False, simulation=False, output=False, filenames=False):
        self.parameters_window = tk.Toplevel(self)
        self.parameters_window.title('Parameters')
        self.parameters_window.resizable(False,False)
        self.parameters_window.message = tk.Message(self.parameters_window,
                                                    text=print_simulation_information(self.neutrals, self.ccp,
                                                                                      self.parameters, self.options,
                                                                                      print_physical=physical,
                                                                                      print_simulation=simulation,
                                                                                      print_output=output,
                                                                                      print_filenames=filenames),
                                                    font=FONT_MESSAGE_SMALL, foreground='blue', width=1500,
                                                    justify=tk.LEFT, relief=tk.RIDGE)
        self.parameters_window.message.grid()
        self.parameters_window.button  = tk.Button(self.parameters_window, text="Dismiss",
                                                   command=self.parameters_window.destroy)
        self.parameters_window.button.grid()
        
    def show_plot_info(self):
        try:
            string = self.runtime_plots.graph_list(expanded=True)
            if (string is None): string = 'No plots'
        except AttributeError:
            string = 'No plots'
        self.plotinfo_window = tk.Toplevel(self)
        self.plotinfo_window.title('Plot Windows Properties')
        self.plotinfo_window.resizable(False,False)
        self.plotinfo_window.text = tk.scrolledtext.ScrolledText(self.plotinfo_window,
                                                                 font=FONT_MESSAGE_SMALL,
                                                                 background='light grey',
                                                                 foreground='blue',
                                                                 width=90,
                                                                 height=30)
        self.plotinfo_window.text.grid()
        self.plotinfo_window.text.insert(tk.INSERT, string)
        self.plotinfo_window.text.configure(state='disabled')
        self.plotinfo_window.button  = tk.Button(self.plotinfo_window, text="Dismiss",
                                                   command=self.plotinfo_window.destroy)
        self.plotinfo_window.button.grid()
            
    def show_runtime_info(self, erase_only=False):
        if erase_only:
            self.text_output.delete('0.0',tk.END)
        else:
            self.text_output.delete('0.0',tk.END)
            self.text_output.insert('0.0',print_runtime_info(self.charges,
                                                             self.neutrals,
                                                             self.ccp,
                                                             self.parameters,
                                                             self.options,
                                                             self.time_before,
                                                             self.electric_bias_before,
                                                             self.n_active_el_before)
                                        )
        self.text_output.grid(row=1, column=0, columnspan=50)
                
    def show_end_simulation(self, text='End simulation', color='black'):
        self.end_simulation_window = tk.Toplevel(self) 
        self.end_simulation_window.title('Ccpla End')
        self.end_simulation_window.resizable(False,False)
        self.end_simulation_window.message = tk.Message(self.end_simulation_window,
                                                        text=text, font=FONT_MESSAGE_BIG,
                                                        foreground=color,
                                                        width=1500,
                                                        justify=tk.CENTER, relief=tk.RIDGE)
        self.end_simulation_window.message.grid()
        self.end_simulation_window.button = tk.Button(self.end_simulation_window, text="Dismiss",
                                                      command=self.end_simulation_window.destroy)
        self.end_simulation_window.button.grid()
        self.end_simulation_window.attributes("-topmost", True)       
                        
    def show_message(self, message, fore='black', back='light grey'):
        self.message_window = tk.Toplevel(self)
        self.message_window.title('Message')
        self.message_window.resizable(False,False)
        self.message_window.message = tk.Message(self.message_window,
                                                 text=message,
                                                 background=back,
                                                 foreground=fore,
                                                 width=1500,
                                                 font=self.options.text_window_font,
                                                 justify=tk.CENTER,
                                                 relief=tk.RIDGE)
        self.message_window.message.grid()
        self.message_window.button  = tk.Button(self.message_window, text="Dismiss", command=self.message_window.destroy)
        self.message_window.button.grid()
        self.message_window.attributes("-topmost", True)

    def draw_cross_sections(self, electrons=False, ions=False, recombination=False, plot_all=False):
        (status, message) = plot_cross_sections(self.neutrals,
                                                electrons=electrons,
                                                ions=ions,
                                                recombination=recombination,
                                                plot_all=plot_all,
                                                dot_points=self.parameters.dot_points) 
        if (status !=0): self.show_message(message, fore='red', back='light yellow')

    def draw_plots(self, force_plot=False):
        self.runtime_plots.update_history(self.charges, self.charges.time)
        if ( force_plot or (self.plot_counter >= self.plot_delay.get()) ):
            self.plot_counter = 0
            if self.plot_history.get():
                self.runtime_plots.plot_n_and_energy_data(self.charges.time)
                self.update_idletasks()
            if self.plot_phase_space.get():
                self.runtime_plots.plot_phase_space_data(self.charges.time, self.charges)
                self.update_idletasks()
            if self.plot_field.get():
                self.runtime_plots.plot_potential(self.charges.time, self.ccp)
                self.update_idletasks()
            if self.plot_distributions.get():
                self.runtime_plots.plot_distributions(self.charges.time, self.charges)
                self.update_idletasks()
            if self.plot_3D_positions.get():
                self.runtime_plots.plot_positions(self.charges.time, self.charges)
                self.update_idletasks()
        self.plot_counter += 1

    def print_plots(self):
        self.runtime_plots.print_graphs('png', SAVE_DIRECTORY_NAME)
    

def main(charges, neutrals, ccp, parameters, options):
    """ Graphical interface for ccpla script """

    root = tk.Tk()
    root.title("Ccpla GUI")
    root.resizable(False, False)
        
    # Get the screen dimensions
    screen_width  = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
       
    # Main window
    main_window = CcplaWindow(charges, neutrals, ccp, parameters, options, screen_width, screen_height, master=root)

    root.mainloop()

              
