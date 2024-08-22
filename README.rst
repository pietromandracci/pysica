
##################################################
*PYSICA*: PYthon tools for SImulation and CAlculus
##################################################

.. contents::

Introduction
============

This package contains a collection of tools developed for some specific simulation and calculus tasks
in some fields of physics, including nonthermal plasma discharges, as well as surface modification and analysis.
The package was developed to fit my own needs, and it is released under the GNU GPL licence in the hope it can be helpful to someone else.
Of course, there are no warranties of any kind, as stated in the licence.

Installing and importing *pysica*
=================================

Distribution method
-------------------

Creating a specific package of *pysica* for each of the main Linux distributions (e.g. Debian, Red Hat, SUSE, and so on) would be too difficult
and time consuming, so I have decided to distribute it using the Python wheel system, that is distribution-independent.
Unfortunately, this installation method is concurrent with the one used by the specific Linux distribution (e.g. apt for Debian),
so many distributions allow to use it only inside a `Python virtual environment <https://docs.python.org/3/library/venv.html>`_,
to prevent the possibility of conflicts between the packages installed by different methods.
The procedure to create a virtual environment and install *pyisica* inside is described in the section
`How to install in a virtual Python environment`_ of this document. However, if you find uncomfortable to use a virtual environment,
you can also install the *pysica* package in a local directory of your choice,
as described in the section `How to install in a local directory`_.

.. note:: The package has been developed and tested in a Linux-based operative system.
          Some subpackages could probably be used under other systems also,
          but *they have not been tested on them* and there is no guarantee that they would work.
          The modules compiled from Fortran are linux libraries ('*.so*' files): if you want to use them in another operating system you need to
          recompile them using the *f2py* program and a Fortran compiler. The directories named *fortran* contain the Fortran source files,
          the compiled modules and the scripts used for the compilation (the name of which always start with 'f2py'), but the options
          used in the scripts to call *f2py* are specific for linux and the `gnu95 <https://gcc.gnu.org/fortran/>`_ Fortran compiler.          


Dependencies
------------

This package depends heavily on `numpy <https://numpy.org/>`_ ,
while some specific modules and packages depend on `scipy <https://scipy.org/>`_ and `matplotlib <https://matplotlib.org/>`_ also.
Some packages make use of `tkinter <https://docs.python.org/3/library/tkinter.html>`_
and of the `gnuplot <http://www.gnuplot.info/>`_ progam, but they should work also without them,
although without some features. 



        
How to install in the global Python environment
-----------------------------------------------

*pysica* is distribuited as a *Python wheel* so, if you have the program *pip* installed on your system
and the Linux distribution you are using allows you to do so, you can type at the teminal::

$ pip install pysica

in this way the Python interpreter will be able to use the *pysica* package regardless of the location from where it is invoked.
If *pip* is not installed on your system, you will have to install it following the specific method of the distribution you are using,
e.g. for Debian-related distibutions (including Ubuntu)::

  $ sudo apt-get install pip

.. note::  In several Linux distributions (including Debian-related ones, like Ubuntu) the operative system does not allow *pip*
           to install software in the main file hierarchy and
           you will get an error message saying "externally managed environment" or something similar.
           In this situation, you can install the package in a virtual Python environment, as described in the section
           `How to install in a virtual Python environment`_, or in a local directory, as described in the
           `How to install in a local directory`_ section.


How to install in a virtual Python environment
----------------------------------------------

Creating a *Python virtual environment* allows you to install *pysica* in distributions that do not allow *pip* to install Python packages
in the main filesystem. You can find detailed instructions on how to create and use a virtual environment
in `this page <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments>`_ of the Python documentation.
In the following, a step-by-step description is provided of the procedure I would recommend to create a virtual environment
for a specific project and install *pysica* inside it.

- Open the terminal and move to the directory you want to use for your project, let's assume the directory name is *myproject*

- Decide a name for the virtual environment you are going to create, since you can create several environments in parallel.
  Let's assume you have chosen the name *envpysica*.
  
- Create the environment using the following command at the Linux console::

    $ python3 -m venv envpysica

  this will create a new directory inside the *myproject* directory, named *envpysica*, inside which a copy of the Python system will be created.

- Now you have to activate the *envpysica* virtual environment by the following command::
    
    $ source ./envpysica/bin/activate

  or equivalently::

    $ . ./envpysica/bin/activate

  where the first dot is a short for the command *source*.
  
  After the activation, the shell prompt should change in this way::

    (envpysica) $

  showing that you are using the *envpysica* environment.
  While *envpysica* is activated, if you call the Python interpreter it will use the packages that are inside the *envpysica*
  directory instead of the ones stored in the global filesystem. Moreover, if you install a Python package using the *pip*
  program, it will be installed in that directory.
  
- You can now install the *pysica* package inside the *envpysica* environment, and its files will be stored in the *envpysica* directory
  instead of in the global Linux file structure. To do this, type the following command at the terminal::

    (envpysica) $ pip install pysica

  the installation process should start and a progress bar should be showed in the terminal.  
  After the installation is completed, you can run the Python interpreter 
  and import the *pysica* module as described in the section `How to import`_.

- When you have finished working with *pysica*, you can exit from the *envpysica* environment by typing at the console the command::

    (envpysica) $ deactivate

  and the shell prompt should return the normal one.

Now, each time you want to use *pysica*, you have to enter the *envpysica* environment by moving to the *myproject* directory and running
the command::
  
    $ . ./envpysica/bin/activate

and then run the Python interpreter.
    

How to install in a local directory
-----------------------------------

If your Linux distribution does not allow you to install the package in the global file structure,
and you do not want to use a virtual environment,
you can install *pysica* in any directory of your system by dowloading the most recent zip or tar.gz archive from the *pysica* 
`GitHub page <https://github.com/pietromandracci/pysica/releases>`_ and unzipping it in a directory of your choice.

A new directory will be created, named *pysica-x.y.z*, where *x.y.z* identifies the version number. 
In order to use *pysica*, you will have to and call the Python interpreter from this *pysica-x.y.z* directory,
so that it will be able to find the package files. 

You can also use the package calling the Python intepreter from another directory,
but you have to create in that directory a symbolic link to the directory named *pysica*,
which is inside the *pysica-x.y.z* directory created during the zip archive extraction.
          

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

Documentation about the modules and packages is available in the docstrings, which can be accessed inside the Python interpreter
using the *help* function, after you have imported them. As an example, to read the docstrings of the subpackage named *analysis* you can type::

  >>> import pysica.analysis
  >>> help(pyisca.analysis)

or::

  >>> from pysica import analysis
  >>> help(analysis)

Note that, due to the Python importing mechanism, in order for the help function to work you must import the specific subpackage 
on which you want help: e.g. if you import the main package *pysica*, the help fucntion will not work on the subpackage *pysica.analysis*.

For some packages, additional documentation can be found in the
`doc <https://github.com/pietromandracci/pysica/tree/master/doc>`_ directory of the *GitHub* repository.  In this case,
a direct link to the documentation is given in the corresponding paragraph of the section `Package structure`_.



Package structure
=================

In the following, the main modules and subpackages are listed.


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
  some general purpose functions, compiled from Fortran using f2py,
  they are collected in the *fmathematics* module;

*mathematics (module)*
  some general purpose mathematical funtions;

*statistics (module)*
  some generic statistics functions;

*pdf (module)*
  some probabilty distribution functions (pdf);

*random_pdf (module)*
  functions useful to generate random numbers following specific pdfs;
  
*physics (module)*
  some general purpose funcions used in generic physics applications;
  
*optics (module)*
  some functions useful for optical applications.


managers (package)
------------------

Contains some modules and packages used to manage input/output of data from/to ascii files,
to print physical quantities managing the unit prefixes, and to plot data by means of the *gnuplot* program.

*io (package)*
  some modules used for generic input-output management;

*data_manager (module)*
  tools to manage data reading and writing from files;

*unit_manager (module)*
  tools to manage the output of numerical data with automatic managment of unit prefixies;

*gnuplot_manager (package)*
  a package to facilitate the use of gnuplot inside Python [#gnuplot_manager]_.

  
.. [#gnuplot_manager] This package is also available as a standalone package (without the rest of *pysica*) on
   `its own GitHub page <https://github.com/pietromandracci/gnuplot_manager>`_,
   where you can find extensive documentation about it.
   


plasma (package)
-------------------

A package containing tools for the simulation of plasma discharges.

*ccpla (package)*
  a package containing scripts, modules, and subpackages used to simulate low pressure capacitively coupled discharges.
  Documentation about this package is given `here <https://github.com/pietromandracci/pysica/tree/master/doc/ccpla/ccpla_manual.rst>`_. 

  

