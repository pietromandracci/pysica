# COPYRIGHT (c) 2020-2026 Pietro Mandracci

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

""" Simulation of a capacitively coupled plasma discharge: information functions """

from pysica.plasma.ccpla.ccpla_defaults import *
from pysica.managers import time_manager


def save_info_to_file(filename_info, options, time, gui=False):

    """ Save information about the process on a text file

        Parameters
        ----------
        filename_info:  name of the file to which information will be saved
        options:        command line options
        time:           a tuple containing in order: clock time,thread time, performance time
        gui:            if True, the GUI has been activated
    """
    
    info_file = open(filename_info, 'w')
    info_file.write('Process PID:              ' + str(PID) + EOL)
    info_file.write('CPU mode:                 ')
    if (options.cpu_multicore):
        info_file.write('multicore' + EOL)
        info_file.write('Number of threads:        ')
        if (options.cpu_threads == 0): info_file.write('auto' + EOL)
        else:                          info_file.write(str(options.cpu_threads) + EOL) 
    else:
        info_file.write('single core' + EOL)
    info_file.write('GUI activated:            ')
    if gui: info_file.write('YES' + EOL)
    else:   info_file.write('NO' + EOL)
        
    info_file.write('Elapsed clock time:       ' + time_manager.print_timestamp2time(time[0], printzeros=True)
                                                 + ' (' + str(time[0]) + ' ns)' + EOL)
    info_file.write('Elapsed thread time:      ' + time_manager.print_timestamp2time(time[1], printzeros=True)
                                                 + ' (' + str(time[1]) + ' ns)' + EOL)
    info_file.write('Elapsed performance time: ' + time_manager.print_timestamp2time(time[2], printzeros=True)
                                                 + ' (' + str(time[2]) + ' ns)' + EOL)
    info_file.close()
    

    
