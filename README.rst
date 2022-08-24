
##################################################
*PYSICA*: PYthon tools for SImulation and CAlculus
##################################################

.. contents::

Introduction
============

This package contains a collection of tools developed for some specific simulation and calculus tasks
in the fields of cold plasma processes, as well as thin-film and surface characterization.


Package structure
=================

In the following, the modules and subpackages are listed.
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


fortran (package)
-----------------

Contains some general purpose modules compiled from Fortran using f2py.


*fmathematics (module)*
  contains some general purpose mathematical functions used in other modules and
  compiled from fortran to improve speed.

  
functions (package)
-------------------

Contains some general purpose functions.


io (package)
------------

Contains some tools for input-output management.


*io_files (module)*
  contains some tools to operate on files;


*io_screen (module)*
  contains some tools to operate on the screen.


managers (package)
------------------

Contains some modules and packages used to manage input/output of data from/to ascii files,
to print physical quantities managing the unit prefixes, and to plot data by means of the *gnuplot* program.

*data_manager (module)*
  tools to manage data reading and writing from files;


*unit_manager (module)*
  tools to manage the output of numerical data with automatic managment of unit prefixies;


*gnuplot_manager (package)*
  a package to facilitate the use of gnuplot inside python [#gnuplot_manager]_.

.. [#gnuplot_manager] *gnuplot_manager* is also available as a standalone package (without the rest of *pysica*) on
  `GitHub <https://github.com/pietromandracci/gnuplot_manager>`_  and
  `PyPi <https://pypi.org/project/gnuplot-manager>`_.


plasmapro (package)
-------------------

A package containing tools for the simulation of plasma discharges.

*ccpla*
  a script to simulate a ccp cold plasma discharge by means of the PIC-MC (1d3v) technique;

*ccpla_analysis*
  a module containing tools to analyze the output data from *ccpla*.
    

Installing and importing *pysica*
=================================


Dependancies
------------

This package depends heavily on `numpy <https://numpy.org/>`_ and `matplotlib <https://matplotlib.org/>`_,
while some specific modules and packages depend on `scipy <https://scipy.org/>`_ also.
Some packages make use of `tkinter <https://docs.python.org/3/library/tkinter.html>`_
and of the `gnuplot <http://www.gnuplot.info/>`_ progam, but they should work also without it,
although without some features. 


How to install
--------------

*pysica* is distribuited as a pypi wheel so, if you have *pip* installed on your system, you can simply type at the console::

$ pip install pysica

In some linux distributions (e.g. Debian-related ones) you will have to install the package inside a python virtual environment,
since the operative system doesn't allow *pip* to install software in the main file hierarchy.  You can find instructions on how to create
a virtual environment `here <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments>`_.

.. note:: The package has been developed and tested for use in linux. Some subpackages could probably be used under other systems also,
          but *they have not been tested on them* and there is no guarantee that they would work.

.. note:: The modules compiled from Fortran are linux libraries ('*.so*' files): if you want to use them in another operating system you need to
          recompile them using the *f2py* program and a Fortran compiler. The directories named *fortran* contain the Fortran source files,
          the compiled modules and the scripts used for the compilation (the name of which always start with 'f2py'), but the options
          used in the scripts to call *f2py* are specific for linux and the `gnu95 <https://gcc.gnu.org/fortran/>`_ Fortran compiler.


How to import
-------------

Once installed, you can import *pysica* using the *import* directive as usual:

>>> import pysica

Or you can import a single mudule or package that you need, such as:

>>> from pysica.managers import gnuplot_manager

or

>>> from pysica.analysis import spectra


Documentation
=============

Documentation about the modules and packages is available in the docstrings.  Additional documentation can be found in the
`doc <https://github.com/pietromandracci/pysica/tree/master/doc>`_ directory of the *GitHub* repository.


