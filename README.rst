
##################################################
*PYSICA*: PYthon tools for SImulation and CAlculus
##################################################

.. contents::

Introduction
============

This package contains a collection of tools developed for some specific simulation and calculus tasks
in some fields of physics, including nonthermal plasma discharges, as well as surface modification and analysis.


Package structure
=================

In the following, the main modules and subpackages are listed.
Additional documentation is available in the docstrings of each module and subpackage.


constants (module)
------------------

Contains some physical constants used in various modules and packages.


parameters (module)
-------------------

Contains some parameters used in various modules and packages.

    
analysis (package)
------------------

Contains some modules to manage distribution functions and data histograms.

*univariate (module)*
  tools for the statistical analysis of univariate samples;

*bivariate (module)*
  tools for the statistical analysis of bivariate samples;

*spectra (module)*
  tools for the analysis of different types of spectra, whith a special focus on:
    - optical data (e.g. transmission spectra) of thin films;
    - surface morphology data (e.g. surface roughness analysis).


  
functions (package)
-------------------

Contains some general purpose functions.

*fortran (package)*
  some general purpose functions, compiled from Fortran using f2py.
  They are collected in the *fmathematics* module.

*mathematics (module)*
  some general purpose mathematical funtions.

*statistics (module)*
  some generic statistics functions.  

*pdf (module)*
  some probabilty distribution functions (pdf).

*random_pdf (module)*
  functions useful to generate random numbers following specific pdfs.
  
*physics (module)*
  some general purpose funcions used in generic physics applications.
  
*optics (module)*
  some functions useful for optical applications.


managers (package)
------------------

Contains some modules and packages used to manage input/output of data from/to ascii files,
to print physical quantities managing the unit prefixes, and to plot data by means of the *gnuplot* program.

*io (package)*
  some modules used for generic input-output management.

*data_manager (module)*
  tools to manage data reading and writing from files;

*unit_manager (module)*
  tools to manage the output of numerical data with automatic managment of unit prefixies;

*gnuplot_manager (package)*
  a package to facilitate the use of gnuplot inside Python [#gnuplot_manager]_.

.. [#gnuplot_manager] *gnuplot_manager* is also available as a standalone package (without the rest of *pysica*) on
  `GitHub <https://github.com/pietromandracci/gnuplot_manager>`_  and
  `PyPi <https://pypi.org/project/gnuplot-manager>`_.


plasma (package)
-------------------

A package containing tools for the simulation of plasma discharges.

*ccpla (package)*
  a package containing scripts, modules, and subpackages used to simulate low pressure capacitively coupled discharges.
    

Installing and importing *pysica*
=================================


Dependencies
------------

This package depends heavily on `numpy <https://numpy.org/>`_ and `matplotlib <https://matplotlib.org/>`_,
while some specific modules and packages depend on `scipy <https://scipy.org/>`_ also.
Some packages make use of `tkinter <https://docs.python.org/3/library/tkinter.html>`_
and of the `gnuplot <http://www.gnuplot.info/>`_ progam, but they should work also without it,
although without some features. 

.. note:: The package has been developed and tested to be used in a Linux-based operative system.
          Some subpackages could probably be used under other systems also,
          but *they have not been tested on them* and there is no guarantee that they would work.

.. note:: The modules compiled from Fortran are linux libraries ('*.so*' files): if you want to use them in another operating system you need to
          recompile them using the *f2py* program and a Fortran compiler. The directories named *fortran* contain the Fortran source files,
          the compiled modules and the scripts used for the compilation (the name of which always start with 'f2py'), but the options
          used in the scripts to call *f2py* are specific for linux and the `gnu95 <https://gcc.gnu.org/fortran/>`_ Fortran compiler.


How to install in the global Python environement
------------------------------------------------

*pysica* is distribuited as a *Python wheel* so, if you have the program *pip* installed on your system, you can simply type at the console::

$ pip install pysica

in this way the Python interpreter will be able to use the *pysica* package regardless of the location from where it is invoked.

.. note::  In some Linux distributions (e.g. Debian-related ones) you will have to install the package inside a Python virtual environment,
           since the operative system doesn't allow *pip* to install software in the main file hierarchy.
           You can find instructions on how to create a virtual environment
           `in this page <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments>`_
           of the Python documentation.


How to install in a local directory
-----------------------------------

You can also install *pysica* in any directory of your system by dowloading the most recent zip archive from the *pysica* 
`GitHub page <https://github.com/pietromandracci/pysica/releases>`_ and unzipping it in a directory of your choice.

A new directory will be created, named *pysica-x.y.z*, where *x.y.z* identifies the version number.
In order to use *pysica* you will have to open a terminal, navigate to this directory, and call the Python interpreter from there.

          

How to import
-------------

Once you have installed *pysica*, you can run the Python interpreter from the console::

$ python3

and then import *pysica* using the *import* directive as usual:

>>> import pysica

Or you can import a single mudule or package that you need, such as:

>>> from pysica.managers import gnuplot_manager

or

>>> from pysica.analysis import spectra


Documentation
=============

Documentation about the modules and packages is available in the docstrings.  Additional documentation can be found in the
`doc <https://github.com/pietromandracci/pysica/tree/master/doc>`_ directory of the *GitHub* repository.


