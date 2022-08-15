# COPYRIGHT 2020 Pietro Mandracci

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

""" Simulation of a capacitively coupled plasma discharge: numerical kernel """

from pysica.plasmapro.ccpla_defaults import *
from pysica.plasmapro.discharge import particle_mover


def kernel(signal_connection, time_connection, data_connection, debug=False):

    """ This function runs parallel to the GUI and calls the particle mover """

    if debug: print('[Kernel]-> ***INSTALLED***\n')

    stay_alive = True
    while stay_alive:
        # Wait until a message arrives on the message pipe
        if debug: print('[Kernel]-> waiting for message [start/quit]')
        message = signal_connection.recv()
        if debug: print('[Kernel]-> message is \"' + message + '\"')        
        # If message is 'quit', exit from main loop
        if (message == 'quit'):
            stay_alive = False            
        # If message is 'start', begin the simulation (if it is anything else, do nothing)
        elif (message == 'start'):            
            if debug: print('[Kernel]-> waiting for initialized data\n')
            (charges, neutrals, ccp, parameters, options) = data_connection.recv()
            if debug: print('[Kernel]-> initialized data read\n')            
            # Start simulation loop
            run_simulation = True
            while run_simulation:                
                # Read the iteration time requested
                if debug: print('[Kernel]-> waiting for iteration time')
                time = time_connection.recv()
                parameters.dt_output = time                  
                if debug: print('[Kernel]-> iteration time is ' +str(parameters.dt_output))
                
                # Call the particle mover
                if debug: print('[Kernel]-> calling particle mover')
                particle_mover.move_particles(charges, neutrals, ccp, parameters, options)
                if (charges.n_active(0) <= 0):
                    if debug: print('[Kernel]-> no more electrons') 
                    run_simulation = False
                elif ( (parameters.sim_duration > 0) and (charges.time >= parameters.sim_duration) ):
                    if debug: print('[Kernel]-> time completed') 
                    run_simulation = False
                # Send to pipes signal of completed iteration and data
                if debug: print('[Kernel]-> sending message [completed]')
                signal_connection.send('completed')
                if debug: print('[Kernel]-> sending data')                        
                data_connection.send((charges, neutrals, ccp, parameters, options))
                if run_simulation:
                    # Read a signal 
                    if debug: print('[Kernel]-> waiting for message [continue/stop]')
                    message = signal_connection.recv()
                    if debug: print('[Kernel]-> message is \"' + message + '\"')
                    if (message == 'stop'):
                        if debug: print('[Kernel]-> simulation stopped')
                        run_simulation = False
                if debug: print('[Kernel]-> end loop (run_simulation = ' + str(run_simulation) +')\n')
            
    if debug: print('[Kernel]-> ***QUITTING***\n')
