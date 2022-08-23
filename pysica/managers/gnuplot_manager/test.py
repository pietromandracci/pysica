#!/usr/bin/python

""" This script is intended to test the functions module.

    If the package has been correctly installed, you should be able to run 
    this script from the system terminal by typing::

    $ python -m gnuplot_manager.test

    or inside the python console by typing::

    >>> from gnuplot_manager.test import main
    >>> main()
"""

import numpy
from .functions import *

def wait():
    return input('Press RETURN to continue...\n')

def main():

    x = numpy.linspace(0, 100, 101)
    a = numpy.linspace(0, 2.7, 271)

    print('\n**Test of gnuplot_manager functions**\n')
    wait()

    print('* -> Function new_plot')
    p2d = new_plot(xpos=480, ypos=100,
                   plot_type='2D', title='TEST 2D PLOT WINDOW', redirect_output=False)
    p3d = new_plot(xpos=500, ypos=120,
                   plot_type='3D', title='TEST 3D PLOT WINDOW', redirect_output=False)
    wait()
   
    print('* -> Function plot2d')
    plot2d(p2d, x, x*x, 'parabola')
    wait()

    print('* -> Function plot2d with replot')
    plot2d(p2d, x, x*x*x/100, 'x^3/100', replot=True)
    wait()

    print('* -> Function plot3d')
    plot3d(p3d, numpy.cos(4*a), numpy.sin(4*a), a, 'elicoidal (I)')
    wait()

    print('* -> Function plot3d with replot')
    plot3d(p3d, numpy.sin(4*a), numpy.cos(4*a), a, 'elicoidal (II)', replot=True)
    wait()
    
    print('* -> Function plot_curves (2D)')
    list1 = [ [x, x*x, 'parabola', 'with points'] , [x, x*x*x/100, 'x^3/100', 'with lines']  ]
    plot_curves(p2d, list1)
    wait()

    print('* -> Function plot_curves (3D)')
    list2 = [ [ numpy.cos(4*a), numpy.sin(4*a), a, 'elicoidal (I)', 'with points'] ,
              [ numpy.sin(4*a), numpy.cos(4*a), a, 'elicoidal (II)','with lines' ]  ]
    plot_curves(p3d, list2)
    wait()

    print('* -> Function plot_function (2d)')
    plot_function(p2d, 'x**2', None)
    wait()

    print('* -> Function plot_function (2d) with replot')
    plot_function(p2d, 'x**3/100', None, replot=True)
    wait()  

    print('* -> Function plot_function (3d)')
    plot_function(p3d, 'x**2+y**2', None)
    wait()

    print('* -> Function plot_function (3d) with replot')
    plot_function(p3d, '(x**3+y**3)/10', None, replot=True)
    wait()

    print('* -> Function plot_functions (2d)')
    list1 = [ ['x**2', None, 'with points'] , ['x**3/10', None, 'with lines']  ]
    plot_functions(p2d, list1)
    wait()

    print('* -> Function plot_functions (3d)')
    list2 = [ ['x**2+y**2', None, 'with points'] , ['(x**3+y**3)/10', None, 'with lines']  ]
    plot_functions(p3d, list2)
    wait()    
    
    print('* -> Function plot_set')
    plot_set(p3d, xmin=1, xmax=10, ymin=1, ymax=1000, logz=True, replot=True)
    wait()

    print('* -> Function plot_label')
    plot_label(p2d, x=20, y=20, label='THIS IS A LABEL', replot=True)
    wait()

    print('* -> Function plot_label (erase=True)')
    plot_label(p2d, label=None, erase=True, replot=True)
    wait()

    print('* -> Function plot_list')
    plot_list()
    wait()

    print('* -> Function plot_list with expanded option')
    plot_list(expanded=True)
    wait()

    print('* -> Function plot_close_all')
    plot_close_all()
    wait()

    print('**End of test**\n')
    
if __name__ == '__main__':
    main()
