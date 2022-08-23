
#############################################
CCPLA: Capacitively Coupled PLAsma simulation
#############################################

.. contents::

Introduction
============

This script allows the simulation of *capacitively coupled* plasma discharges at low pressure 
in a simple geometric configuration (square parallel electrodes).  It makes use of the PIC-MC (*particle in cell* plus *montecarlo*)
technique in the *1d3v* mode (one dimension for the electric field and 3 dimensions for the particle movement).


Installing and running
======================


Dependancies
------------

The ccpla script depends heavily on `numpy <https://numpy.org/>`_.
The GUI depends on `tkinter <https://docs.python.org/3/library/tkinter.html>`_ also,
and makes use of the the `gnuplot <http://www.gnuplot.info/>`_
progam to plot data 'on the fly' during simulation when using the GUI.

The script should work without *tkinter* and *gnuplot* installed, but only in batch mode (no GUI and no plots).          


How to install
--------------

The *ccpla* script is part of the *pysica* package, which is distribuited as a pypi wheel so,
if you have *pip* installed on your system, you can simply type at the console::
             
$ pip install pysica

In Debian-related linux distributions you will have to install the package inside a python virtual environment, since the operative
system doesn't allow *pip* to install software in the main file hierarchy.  You can find instructions on how to create
a virtual environment `here <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments>`_.

.. note:: The package has been developed and tested for use in linux. Some subpackages could probably be used under other systems also,
          but *they have not been tested on them* and there is no guarantee that they would work.

.. note:: The modules compiled from Fortran are linux libraries ('*.so*' files): if you want to use them in another operating system you need to
          recompile them using the *f2py* program and a Fortran compiler. The directories named *fortran* contain the Fortran source files,
          the compiled modules and the scripts used for the compilation (the name of which always start with 'f2py'), but the options
          used in the scripts to call *f2py* are specific for linux and the `gnu95 <https://gcc.gnu.org/fortran/>`_ fortran compiler.


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
`data <https://github.com/pietromandracci/pysica/tree/master/data/ccpla>`_ directory of the *pysica* *GitHub* repository.
You can copy them in the directory from which you want to run the script and then modify them at your will,
as explained in the comments inside each file.


How to run the script
---------------------

Once installed, you can run the *ccpla* script using the following command::

$ python3 -m pysica.plasmapro.ccpla [options]

where *[options]* states for a list of options you can give to the script.  A list of the available options can be obtained by::

$ python3 -m pysica.plasmapro.ccpla -h

The *-g* options runs the GUI instead of the text-based batch script::

$ python3 -m pysica.plasmapro.ccpla -g




