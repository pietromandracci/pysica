=============
Version 0.5.0
=============

* In the ccpla subpackage some changes were made in the format the data files are saved, moreover now the cross section files
  are loaded from a the ccpla.sigma subdiretory, allowing for a more clean working directory

=============
Version 0.4.1
=============

* Cleaning up of unused files and minor updates

=============
Version 0.4.0
=============

* The subpackages structure was rearranged in more rational way (or at least I hope so)

=============
Version 0.3.0
=============

* The fortran module of the plasmapro.discharge subpackage was modified to allow multithread execution
  of the most cpu-intensive part:
  it can be activated running the ccpla script with the -m comman line argument

* A bug was fixed in the fortran module of the plasmapro.discharge subpackage (missing array initialization)

* The ccpla.py and ccpla_gui.py scripts have been modified to allow that distribution data are saved only a reduced
  number of times to reduce space consumption on disk

=============
Version 0.2.2
=============

* Minor changes in the documentation

=============
Version 0.2.1
=============

* Update of the file headers with copyright info

=============
Version 0.2.0
=============

* Major changes in the script ccpla
  - ion-neutral scattering is now calculated in the center of mass frame of reference
  - more data are saved to file during the simulation, including potential distribution and electric current

* Changes in the script ccpla_analysis
  - it is now possible to plot other data, such as the  potential distribution
 
=============
Version 0.1.0
=============

* Major changes in the script ccpla
  - electron impact excitation processes have been added to the simulation
  - some changes in the GUI have made to correct some erros that did happen when reloading configuration files from the GUI  

=============
Version 0.0.1
=============

First release
