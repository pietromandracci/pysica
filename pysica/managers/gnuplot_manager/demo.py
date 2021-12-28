#!/usr/bin/python

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

""" This is a small demo script, to show the main features of this package.

    If the package has been correctly installed, you should be able to run 
    this script from the system terminal by typing::

    $ python -m gnuplot_manager.demo

    or inside the python console by typing::

    >>> from gnuplot_manager.demo import main
    >>> main()

    The demo will run without need of input and last about 35 seconds.    
"""



from time import sleep
import numpy
from .functions import *

def wait():
    sleep(3)
    
def main():

    x = numpy.linspace(0, 100, 51)

    print('\n***Demo of the gnuplot_manager package***\n')
    
    print('In this demo, gnuplot output is silenced\n')
    p2d = new_plot(xpos=600, ypos=100, width=600, height=450,
                   plot_type='2D', title='2D PLOT WINDOW DEMO', redirect_output=None)        
    p3d = new_plot(xpos=680, ypos=200, width=600, height=450,
                   plot_type='3D', title='3D PLOT WINDOW DEMO', redirect_output=None)    
    
    print(' 1) -> Plotting a parabola from data\n')
    plot2d(p2d, x, x*x, 'parabola')
    wait()

    print(' 2) -> Plotting a parabola as a function\n')
    plot_function(p2d, 'x**2', replot=True)
    wait()
    
    print(' 3) -> Changing the scale limits\n')
    plot_set(p2d, xmin=1, xmax=5, ymin=1, ymax=25, replot=True)
    wait()

    print(' 4) -> Setting logarithmic x and y axes\n')
    plot_set(p2d, logx=True, logy=True, replot=True)
    wait()

    print(' 5) -> Resetting ranges and scales\n')
    plot_set(p2d, resetx=True, resety=True, logx=False, logy=False, replot=True)
    wait()    

    print(' 6) -> Adding a label to the plot\n')
    plot_label(p2d, x=20, y=20, label='THIS IS A LABEL', replot=True)
    wait()

    print(' 7) -> Removing the label\n')
    plot_label(p2d, label=None, erase=True, replot=True)
    wait()

    print(' 8) -> Plotting another curve with volatile data\n')
    plot2d(p2d, x, x*x*x/100, 'x^3/100 (volatile)', volatile=True, replot=True)
    wait()            

    print(' 9) -> Plotting a 3D curve from data\n')
    plot3d(p3d, x, x, 2*x*x, '3D curve')
    wait()

    print('10) -> Plotting a 3D surface as a function\n')
    plot_function(p3d, 'x**2+y**2', replot=True)
    wait()

    print('11) -> Setting logarithmic z axis\n')
    plot_set(p3d, logz=True, replot=True)
    wait()

    print('12) -> Listing the plots\n')
    plot_list()
    wait()

    print('13) -> Closing the 2d plot\n')
    plot_close(p2d)
    wait()

    print('14) -> Listing the plots again')
    plot_list()
    wait()    

    print('15) -> Closing all the plots, removing files and directories\n')
    plot_close_all()

    print('***End of demo***')
    
if __name__ == '__main__':
    main()
