
###############################################
*CCPLA*: Capacitively Coupled PLAsma simulation
###############################################

.. contents::

Introduction
============

The *ccpla* script allows the simulation of *capacitively coupled* plasma discharges at low pressure 
in a simple geometric configuration (square parallel electrodes).  It makes use of the PIC-MC (*particle in cell* plus *montecarlo*)
technique in the *1d3v* mode (one dimension for the electric field and 3 dimensions for the particle movement).


Installing
==========

Distribution method
-------------------

The *ccpla* script is part of the *pysica* package, which is distribuited as a *Python wheel*, and can be installed using the program *pip*. 
However, in some Linux distributions (e.g. Debian-related ones, including Ubuntu) the operative system doesn't allow *pip* to install software
in the main file hierarchy and  you will get an error message saying "externally managed environment" or something similar.
In this situation, you can install the package inside a *Python virtual environment* or in a local directory,
as described in the *pysica* documentation page, that you can find following the link given below.

Detailed instructions on how to install the *pysica* package can be found in the
`this page <https://github.com/pietromandracci/pysica/tree/master?tab=readme-ov-file#installing-and-importing-pysica>`_
of the *pysica* documentation on GitHub.

           
Dependencies
------------

The ccpla script depends heavily on `numpy <https://numpy.org/>`_.
The GUI depends on `tkinter <https://docs.python.org/3/library/tkinter.html>`_ also,
and makes use of the the `gnuplot <http://www.gnuplot.info/>`_ progam to plot data 'on the fly' during the simulation.

The script should work even if *tkinter* and *gnuplot* are not installed, but without using the GUI (i.e. without the *-g* option).


Fortran compiled modules
------------------------

The directory *pysica/plasma/ccpla/discharge/fortran* of the *pysica* package contains the Fortran source files,
the compiled modules for single core an multicore operation, and the scripts used for their compilation, 
which are named *f2py-fmodule* and *f2py-fmodule-parallel* respectively.
Note that the options used in those scripts to call *f2py* are specific for Linux and the
`GNU Fortran compiler <https://gcc.gnu.org/fortran/>`_.

.. note:: The script has been developed and tested to be used in a Linux environment. It may work in other systems also,
          but *it has not been tested on them* and there is no guarantee that it would work.
          The module compiled from Fortran is a linux library ('*.so*' file): if you want to use *ccpla* in another operating system,
          you need to recompile it using the *f2py* program and a Fortran compiler.


Running
=======


Data files
----------

The *ccpla* script needs to read some data (configuration data and cross-section data) from external files
in order to perform the simulation, and will give an error message if unable to find them.

Two plain ASCII files are needed for configuration, and must be present in the directory from which the script is run
(*working directory*):

*ccpla.conf*
  contains the main physical and simulation parameters;

*ccpla.neutrals*
  contains the data about the gas mixture to simulate.

Additionaly, a set of *.csv* (comma separated value) files are needed, which contain cross-section data about the types of collision
considered in the simulation.
These files must be present in a subdirectory named *ccpla.sigma* inside the working directory.

A set of examples files for the simulation of He, Ar and O2 plasmas (or their mixtures) can be downloded from the 
`data/ccpla <https://github.com/pietromandracci/pysica/tree/master/data/ccpla>`_ directory of the *pysica* *GitHub* repository.
You can create a subdirectory named *ccpla.sigma* the directory from which you want to run the script, and then copy these files inside it.
Of course, you can modify them at your will, as explained in the comments inside each file.

.. note:: the cross sections tabulated in these files have been taken from some publications or public databases,
          and are given as example only.  There is no guarantee that they are suited for a specific task.

How to run the script
---------------------

Once you have installed the *pysica* package, you can run the *ccpla* script invoking the Python interpreter.

- If you have installed *pysica* using the *pip install pysica* command in the global environment,
  you can run the *ccpla* script from any directory of your filesystem.
  
- If you have installed *pysica* using the *pip install pysica* command in a *Python virtual environment*,
  you can run the *ccpla* script from any directory of your filesystem, *after having activated the virtual environment*
  as described in 
  `this section <https://github.com/pietromandracci/pysica/tree/master?tab=readme-ov-file#How-to-install-in-a-virtual-Python-environment>`_
  of the *pysica* documentation.
  
- If you have installed *pysica* unzipping the package archive in a local directory,
  you will be able to run *ccpla* from the *pysica-x.y.z* directory created by the unzippig process, as described in
  `this other section <https://github.com/pietromandracci/pysica/tree/master?tab=readme-ov-file#How-to-install-in-a-local-directory>`_
  of the *pysica* documentation.
  In order to run *ccpla* from another directory, you have to create inside it a link to the directory named *pysica*,
  that is located inside the *pysica-x.y.z* directory.

Once you are in the proper directory, you can run the *ccpla* script from the console using the following command::

$ python3 -m pysica.plasma.ccpla.ccpla [options]

where *[options]* states for a list of options you can give to the script.  A list of the available options can be obtained by typing::

$ python3 -m pysica.plasma.ccpla.ccpla -h

The *-g* option runs the GUI instead of the text-based script::

$ python3 -m pysica.plasma.ccpla.ccpla -g



Command line options 
---------------------

Here's a list of the available command line options and their meaning

*-h, --help*
    show a list of the available options
    
*-p, --print-only*
    print the simulation parameters on the screen, but do not start the simulation,    
    this option rises an error when calling the GUI (-g option)
    
*-s, --save-defaults*
    write the default parameters to a file named "ccpla.defaults" and then exit the program
    
*-b, --batch-mode*
    suppress all input from user, so that the script can be run in background, 
    this option has no effect when calling the GUI (-g option)
    
*-g, --gui-mode*
    start the GUI

*-m, --multicore*
    run the Fortran module in parallel mode, using more than one CPU core
    
*-W TEXT_WINDOW_WIDTH, --text-window-width=TEXT_WINDOW_WIDTH*
    set the width of the GUI window expressed in characters, accepted values are in the range [120..200] (default=160), 
    this option has an effect only while calling the GUI (-g option)
    
*-H TEXT_WINDOW_HEIGHT, --text-window-height=TEXT_WINDOW_HEIGHT*
    set the height of GUI window expressed in characters, accepted values are in the range [20..80] (default=39), 
    this option has an effect only while calling the GUI (-g option)
    
*-F TEXT_WINDOW_FONT_SIZE, --text-window-font=TEXT_WINDOW_FONT_SIZE*
    set the font size in the GUI window, accepted values are in the range [6..18] (default=12), 
    this option has an effect only while calling the GUI (-g option)
    
*-o, --redirect-output*
    redirect the screen output to a file named 'ccpla_output.log'
    
*-e, --redirect-errors*
    redirect the error messages to a file named 'ccpla_errors.log'
    
*-v VERBOSITY, --verbosity=VERBOSITY*
    set the verbosity level of the text output [0..3] (default=1), 
    this option has no effect when calling the GUI (-g option)
    
*-x, --graph-xsec*
    plot cross sections graphs before starting the text-based script, 
    this option rises an error when calling the GUI (-g option)            
    
*-d DEBUG_LEV, --debug-level-python=DEBUG_LEV*
    Python debug level [0..2] (default=0), this is used for debugging purposes only

*-D DEBUG_LEV_FOR, --debug-level-fortran=DEBUG_LEV_FOR*
    Fortran debug level [0..3] (default=0), this is used for debugging purposes only


Graphical User Interface
========================

The GUI is run by using the *-g* options when callig the script::

  $ python3 -m pysica.plasma.ccpla.ccpla -g

When the GUI is started, it activates a main graphic window and a numerical kernel. The latter is an independent process,
which runs parallell to the GUI and is responsible for calling the Fortran-compiled module when a simulation cycle is requested.
The simulation data are transferred between the GUI and the kernel by means of pipes.

Main window
-----------

When the GUI starts, a main window is rised, together with a small window with licencing information,
which can be closed by pressing the "Dismiss" button.

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-main+splash.png
   :width:  809
   :height: 436

The main window has a menu bar at the top and some buttons and sliders at the bottom.

Menu bar
--------

The menu bar includes the following drop-down menus: *File*, *Parameters*, *Runtime Plots*, and *Help*.


*File*
,,,,,,

The *File* menu shows the following options:

*Reload configuration files*
    reload the content of the *ccpla.conf*. Note that some changes in the file may take effect only after the program is restarted,
    as explained in the comments inside the file itself.  The file *ccpla.neutrals* is *not* reloaded by this command,
    since any change to it becomes effective only after the program is restarted
    
*Edit configuration files*
    open the *ccpla.conf* file in an external editor. The file is reloaded after closing the editor.  The name of the editor to use is
    registered in the variable *EDITOR_NAME* in the file *pysica/plasma/ccpla/ccpla_defaults.py*, presently is *mousepad*.
    If the editor is not installed, an error window is opened

*Quit*
    exit the program.

.. note:: The *File* menu is not active while the simulation is running: in that case you have to press the *Pause* button,
   then the *STOP/KILL* button, in order to stop the simulation and have the menu active again.

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-main-menu-file.png


*Parameters*
,,,,,,,,,,,,

The *Parameters* menu shows the following options:

*Show physical parameters*
    open a window with the physical parameters of the discharge
    
*Show simulation parameters*
    open a window with the parameters used in the simulation

*Show output parameters*
    open a window with the parameters used for the data output

*Show output filenames*
    open a window with the names of the files where simulation data are saved [#a]_

*Show gas properties*
    open a window with the gas properties

*Show e-/neutral impact cross sections*
    open some gnuplot windows with the cross section plots for electron impact

*Show e-/ion recomb cross sections*
    open some gnuplot windows with the cross section plots for electron/ion recombination [#b]_

*Show ion/neutral impact cross sections*
    open some gnuplot windows with the cross section plots for ion impact

*Show e-/neutral impact parameters*
    open some gnuplot windows with other impact parameters (e.g. collision frequencies) for electron collisions

*Show ion/neutral impact parameters*
    open some gnuplot windows with other impact parameters (e.g. collision frequencies) for ion collisions

*Show e-/ion recomb cross parameters*
    open some gnuplot windows with other parameters (e.g. collision frequencies) for electron/ion recombination [#c]_
    

.. [#a] This option is activated after the *RESET* button has been pressed, and only if the simulation
        parameter *save_delay* in the file *ccpla.conf* is not zero.       

.. [#b] This option is activated only if the simulation parameter *isactive_recomb* in the file *ccpla.conf* is not zero.

.. [#c] This option is activated only if the simulation parameter *isactive_recomb* in the file *ccpla.conf* is not zero.
               
.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-main-menu-parameters.png
         

*Runtime Plots*
,,,,,,,,,,,,,,,

The *Runtime plots* menu allows to select which plots are shown during the simulation:

*Select all*
    select all plots for run-time plotting

*Unselect all*
    unselect all plots for run-time plotting

*Mean el energy and number vs time*
    number of electrons (real and computational) and mean electron energy vs simulation time (2 plots)

*Phase space plots*
    electron and ion energy vs angle and vs z-coordinate (4 plots)
    
*Electric potential and charge*
    electric charge and electric potential vs z-coordinate (2 plots)

*EEDF and IEDF*
    electron/ion energy distribution functions (2 plots)

*3D e- and ion positions*
    three dimensional plots of electron and ion positions (2 plots)

.. note:: if some active plots are deactivated while the simulation is running, they are not removed from the screen,
   but they are no longer refreshed until they are re-activated

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-main-menu-plots.png

*Help*
,,,,,,

The *Help* menu shows the following options.

*Online documentation (open in browser)*
    opens the online documentation (this file) inside a web browser. The name of the browser to use is
    registered in the variable *BROWSER_NAME* in the file *pysica/plasma/ccpla/ccpla_defaults.py*, presently *firefox*.
    If the browser is not installed, an error window is opened.

*About*
    shows a window with licencing information

.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-main-menu-help.png


Buttons and sliders
-------------------

The following buttons are positioned at the botton of the main window, each one of them may be inactive (and appear in grey) in some
situations:

*RESET*
    initializes the simulation data and plots, it is inactive while the simulation is running

*START*
    starts the simulation, it is activated after *RESET* has been pressed and becomes inactive after the simulation has started

*Pause / Continue*
    pauses the simulation or continues it after it has been paused, the button label changes properly

*STOP / KILL*
    stops the simulation. It is active only while the simulation is paused. If the button label changes froom *STOP* to *KILL*,
    the program is waiting for the kernel to finish the calculations for a simulation cycle
    (which is performed by the Fortran-compiled module) and can be interrupted by killing the kernel process only.
    A confirmation window is opened before killing the kernel.

In the bottom part of the main window there are two sliders also, by which it is possibile to change how often the output
data are shown on the text window and on the runtime plots. Beside them, there is a small text area in which the values
of these parameters are written, together with the timestep.  The latter can not be changed during the simulation,
but is determined by the *dt* entry in the *ccpla.conf* file.


.. image:: https://raw.githubusercontent.com/pietromandracci/pysica/master/doc/ccpla/figure_gui-buttons.png
