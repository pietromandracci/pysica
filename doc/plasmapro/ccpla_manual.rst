
#############################################
CCPLA: Capacitively Coupled PLAsma simulation
#############################################

.. contents::

Introduction
============

This script allows the simulation of *capacitively coupled* plasma discharges at low pressure 
in a simple geometric configuration (square parallel electrodes).  It makes use of the PIC-MC (*particle in cell* plus *montecarlo*)
technique in the *1d3v* mode (one dimension for the electric field and 3 dimensions for the particle movement).


Installing
==========


Dependancies
------------

The ccpla script depends heavily on `numpy <https://numpy.org/>`_.
The GUI depends on `tkinter <https://docs.python.org/3/library/tkinter.html>`_ also,
and makes use of the the `gnuplot <http://www.gnuplot.info/>`_
progam to plot data 'on the fly' during simulation when using the GUI.

The script should work without *tkinter* and *gnuplot* installed, but without using the GIU (i.e. without the *-g* option).


How to install
--------------

The *ccpla* script is part of the *pysica* package, which is distribuited as a pypi wheel so,
if you have *pip* installed on your system, you can simply type at the console::
             
$ pip install pysica

In some linux distributions (including Debian-related ones, such as Ubuntu) you will have to install the package
inside a python virtual environment, since the operative system doesn't allow *pip* to install software
in the main file hierarchy.
You can find instructions on how to create
a virtual environment `here <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments>`_.

.. note:: The package has been developed and tested for use in linux. Some subpackages could probably be used under other systems also,
          but *they have not been tested on them* and there is no guarantee that they would work.

.. note:: The modules compiled from Fortran are linux libraries ('*.so*' files): if you want to use them in another operating system you need to
          recompile them using the *f2py* program and a Fortran compiler. The directories named *fortran* contain the Fortran source files,
          the compiled modules and the scripts used for the compilation (the name of which always start with 'f2py'), but the options
          used in the scripts to call *f2py* are specific for linux and the `gnu95 <https://gcc.gnu.org/fortran/>`_ fortran compiler.

Running
=======


Data files
----------

The *ccpla* script requires to read some data (configuration data and cross-section data) from external files
in order to perform the simulation: these files must be present in the directory from which the script is run,
otherwise it will give an error message.

Two plain ASCII files are needed for configuration:

*ccpla.conf* file
  contains the main physical and simulation parameters;

*ccpla.neutrals*
  contains the data about the gas mixture to simulate.

Additionaly, a set of *.csv* (comma separated value) files are needed, which contain cross-section data about the gases.

A set of files for the simulation of He, Ar and O2 plasmas (or their mixtures) can be downloded from the 
`data/ccpla <https://github.com/pietromandracci/pysica/tree/master/data/ccpla>`_ directory of the *pysica* *GitHub* repository.
You can copy them in the directory from which you want to run the script and then modify them at your will,
as explained in the comments inside each file.


How to run the script
---------------------

Once installed, you can run the *ccpla* script using the following command::

$ python3 -m pysica.plasmapro.ccpla [options]

where *[options]* states for a list of options you can give to the script.  A list of the available options can be obtained by::

$ python3 -m pysica.plasmapro.ccpla -h

The *-g* option runs the GUI instead of the text-based script::

$ python3 -m pysica.plasmapro.ccpla -g


Command line options 
---------------------

Here's a list of the available command line options and their meaning

*-h, --help*
    show a list of the available options
    
*-p, --print-only*
    prints the simulation parameters on the screen, but does not start the simulation,    
    this option will rise an error when calling the GUI (-g option)
    
*-s, --save-defaults*
    write the default parameters to a file named "ccpla.defaults" and then exit the program
    
*-b, --batch-mode*
    suppress all input from user, so that the script can be run in background, 
    this option has no effect when calling the GUI (-g option)
    
*-g, --gui-mode*
    start the GUI
    
*-W TEXT_WINDOW_WIDTH, --text-window-width=TEXT_WINDOW_WIDTH*
    set the width of the GUI window expressed in characters [120..200] (default=160), 
    this option has an effect only while calling the GUI (-g option)
    
*-H TEXT_WINDOW_HEIGHT, --text-window-height=TEXT_WINDOW_HEIGHT*
    set the height of GUI window expressed in characters GUI [20..80] (default=39), 
    this option has an effect only while calling the GUI (-g option)
    
*-F TEXT_WINDOW_FONT_SIZE, --text-window-font=TEXT_WINDOW_FONT_SIZE*
    set the font size in the GUI window [6..18] (default=12), 
    this option has an effect only while calling the GUI (-g option)
    
*-o, --redirect-output*
    redirect the screen output to a file named 'ccpla_output.log'
    
*-e, --redirect-errors*
    redirect the error messages to a file named 'ccpla_errors.log'
    
*-v VERBOSITY, --verbosity=VERBOSITY*
    set the verbosity level of the text output [0..3] (default=1), 
    this option has no effect when calling the GUI (-g option)
    
*-d DEBUG_LEV, --debug-level-python=DEBUG_LEV*
    Python debug level [0..2] (default=0)
    
*-D DEBUG_LEV_FOR, --debug-level-fortran=DEBUG_LEV_FOR*
    Fortran debug level [0..3] (default=0)
    
*-x, --graph-xsec*
    plot cross sections graphs before starting the program, 
    this option will rise an error when calling the GUI (-g option)


Graphical User Interface
========================

The GUI is run by using the *-g* options when callig the script

$ python3 -m pysica.plasmapro.ccpla -g

Main window
-----------

When the GUI starts, a main window is rised, together with a small window with licencing information,
which can be closed by pressing the "Dismiss" button.

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-main+splash.png
   :width:  809
   :height: 436

The main window has a menu on the top part and some buttons on the bottom.


The *File* menu
---------------

The *File* menu shows the following options:

*Reload configuration files*
    reload the content of the *ccpla.conf*
    
*Edit configuration files*
    open the *ccpla.conf* file in an external editor

*Quit*
    exit the program

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-main-menu-file.png


The *Parameters* menu
---------------------

The *Parameters* menu shows the following options:

*Show physical parameters*
    open a window with the physical parameters of the discharge
    
*Show simulation parameters*
    opne a window with the parameters used in the simulation

*Show output parameters*
    open a window with the parameters used for the data output 

*Show gas properties*
    open a window with the gas properties

*Show e-/neutral impact cross sections*
    open some gnuplot windows with the cross section plots for electron impact

*Show ion/neutral impact cross sections*
    open some gnuplot windows with the cross sections plots for ion impace

*Show e-/neutral impact parameters*
    open some gnuplot windows with other impact parameters (e.g. collision frequency) for electron collisions

*Show ion/neutral impact parameters*
    open some gnuplot windows with other impact parameters (e.g. collision frequency) for ion collisions 

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-main-menu-parameters.png
         

The *Runtime Plots* menu
------------------------

The *Runtime plots* menu allows to select which plots are shown during the simulation:

*Select all*
    select all plots for run-time plotting

*Unselect all*
    unselect all plots for run-time plotting

*Mean el energy and number vs time*
    energy and electron number as a function of simulation time (2 plots)

*Phase space plots*
    electron and ion energy vs angle and vs position along the z-axis (4 plots)
    
*Electric potential and charge*
    electric charge and electric potential as a function of the z-position (2plots)

*EEDF and IEDF*
    electron/ion energy distribution functions

*3D e- and ion positions*
    three dimensional plots of electron and ion positions
    
    



.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-main-menu-plots.png

The *Help* menu
---------------

The *Help* menu shows the following options.

*Online documentation (open in browser)*
    opens the online documentation (this file) inside a web browser

*About*
    shows a window with licencing information

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-main-menu-help.png


Buttons and sliders
-------------------

The following buttons are positioned at the botton of the main window, each one of them may be inactive (and appear in grey) in some
situations:

*RESET*
    pushing the button loads prepares the program for the simulation to start, it is inactive while the simulation is running

*START*
    starts the simulation, it is activated after *RESET* has been pressed and becomes inactive after the simulation has started

*Pause / Continue*
    pauses the simulation or continues it after it has been paused, the button label changes properly

*STOP / KILL*
    it is active only while the simulation is paused, if the label is *KILL* than the program is waiting for the kernel
    to finish the calculations for a simulatiom cycle (which is performed by the Fortran-compiled module) and can be
    interrupted by killing the kernel process only

In the bottom part of the main window there are two sliders also, by which it is possibile to change how often the otput
data are shown on the text window and on the runtime plots.    


.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/plasmapro/figure_gui-buttons.png
