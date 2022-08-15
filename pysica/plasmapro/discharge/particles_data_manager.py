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

""" Functions for reading the particle properties from an ASCII file.

    This module contains functions for managing the import of particle properties from
    properly formatted ASCII files.

    Documentation is also available in the docstrings.
"""


# +--------------------------------+
# | Import required Python modules |
# +--------------------------------+

from pysica.parameters import *
from pysica.constants import *
from pysica.plasmapro.ccpla_defaults import *

# +--------------------------------------------------+
# | Read target particles properties from ascii file |
# +--------------------------------------------------+

def read_file_ntypes(filename, particles):
        """Gets the number of neutral types from the neutrals ascii file.

                Parameters
                ----------

                filename:  name of file containing neutral gases characteristics
                particles: an instance of the class TargetParticles defined in the module target_particles.py

                Returns
                -------

                (status, message)
                        status:  0 = no error
                                 1 = could not open file
                                 2 = no gas type given
                        message: a string containing an error message or 'Ok'
        """

        status, message  = 0, 'OK'              # This variable is used for error checking      

        try:
                f_input = open(filename,'r')
        except IOError:
                status  = 1
                message = 'could not open file'
                return (status, message)

        line            = 0     # Line of the ASCII file
        particles.types = 0
        for string in f_input:
                line = line + 1
 
                # Cut newline charater from the end of line
                string = string.strip(EOL)

                # Cut line portion at right of '#' char
                string = string.split('#')[0]

                # Cut spaces from beginning and end of line
                string = string.strip()

                # If the string is not empty, then the line defines a neutral type
                if (len(string)!=0):
                        particles.types= particles.types +1 

        f_input.close()

        if ( particles.types == 0 ):
                status  = 2
                message = 'no gas type was defined'

        return (status, message)



def read_file_properties(filename, sep, particles, debug=False):
        """Reads the neutrals characteristics from a ascii file and calculates some properties.

                        Parameters
                        ----------

                        filename:         name of file containing neutral gases characteristics
                        sep:              character used in the ASCII file to separate values (usually the tab char)
                        particles:        an instance of the class TargetParticles defined in the module target_particles.py

                        Initialized data attributes
                        ---------------------------

                        Returns
                        -------

                        (status, message)
                                status:  0 = no error
                                         1 = could not open file
                                         2 = missing values
                                         3 = numerical value expected
                                         4 = positive value expected
                                         5 = unknown molecule type
                                         6 = number of types exceedes max allowed
                                message: a string containing an error message or 'Ok'
                """

        status, message  = 0, 'OK'              # This variable is used for error checking      


        try:
                f_input = open(filename,'r')
        except IOError:
                status  = 1
                message = 'could not open file'
                return (status, message)

        line = 0        # Line of the ASCII file
        row  = 0        # Row of the file (excluding empty lines and comments)

        if debug: print('\nReading gas data from file :\"' + filename + '\"\n')

        for string in f_input:

                line = line + 1
                if debug: print("Line: " +str(line)+ " = \'" +string+ "\'")
 
                # Cut newline character from the end of line
                string = string.strip(EOL)

                # Cut line portion at right of '#' char
                string = string.split('#')[0]

                # Cut spaces from beginning and end of line
                string = string.strip()

                # Ignore the line if it is empty
                if (len(string)!=0):

                        # Get a list of the strings in the line
                        l = string.split(sep)
                        if (len(l) < 8):  # 8 is the number of strings expected for 1 excitation process and no dissociation
                                status  = 2
                                message = "line " + str(line) + ", missing values"
                                f_input.close()
                                return (status, message)
                        if debug: print("line, row, l = ", line, row, l)

                        # Get neutral name
                        particles.names[row] = l[0].strip()

                        # Get neutral molecule type
                        particles.molecule_type[row] = l[1].strip()
                        if (particles.molecule_type[row] not in MOLECULE_TYPES):        
                                status  = 5
                                message = ( "line "
                                            + str(line)
                                            + ", unknown molecule type "
                                            + "\""
                                            + str(particles.molecule_type[row])
                                            + "\" for "
                                            + str(particles.names[row]) )
                                f_input.close()
                                return (status, message)

                        # Get gas flow in arbitrary units
                        try:
                                value = float(l[2].strip())
                        except ValueError:
                                status  = 3
                                message = ( "line "
                                            + str(line)
                                            + ", numerical value expected for "
                                            + particles.names[row]
                                            + " gas flow" )
                                f_input.close()
                                return (status, message)
                        if (value <= 0):
                                status  = 4
                                message = ( "line "
                                            + str(line)
                                            + ", positive value expected for "
                                            + str(particles.names[row])
                                            + " gas flow" )
                                f_input.close()
                                return (status, message)                                
                        particles.gas_flow[row] = value

                        # Get neutral mass
                        try:
                                value = float(l[3].strip())
                        except ValueError:
                                status  = 3
                                message = "line "+str(line)+", numerical value expected for "+\
                                          particles.names[row]+" molecular mass"
                                f_input.close()
                                return (status, message)
                        if (value <= 0):
                                status  = 4
                                message = "line "+str(line)+", positive value expected for "+\
                                          particles.names[row]+" molecular mass"
                                f_input.close()
                                return (status, message)                        
                        particles.mass[row] = value

                        # Get secondary emission coefficient
                        try:
                                value = float(l[4].strip())
                        except ValueError:
                                status  = 3
                                message = "line "+str(line)+", numerical value expected for "+\
                                          particles.names[row]+" secondary emission coefficient"
                                f_input.close()
                                return (status, message)
                        if (value <= 0):
                                status  = 4
                                message = "line "+str(line)+", positive value expected for "+\
                                          particles.names[row]+" secondary emission coefficient"
                                f_input.close()
                                return (status, message)                                
                        particles.secondary_emission[row] = value

                        # Get first ionization energy
                        try:
                                value = float(l[5].strip())
                        except ValueError:
                                status  = 3
                                message = "line "+str(line)+", numerical value expected for "+\
                                          particles.names[row]+" ionization energy"
                                f_input.close()
                                return (status, message)
                        if (value <= 0):
                                status  = 4
                                message = "line "+str(line)+", positive value expected for "+\
                                          particles.names[row]+" ionization energy"
                                f_input.close()
                                return (status, message)                                
                        particles.ionization_energy[row] = value

                        # Get excitation energies
                        try:
                                value = int(l[6].strip())
                        except IndexError:
                                status  = 2
                                message = "line "+str(line)+", missing values"
                                f_input.close()
                                return (status, message)
                        except ValueError:
                                status  = 3
                                message = "line "+str(line) + \
                                          ", numerical value expected for number of excitation types "+\
                                          "of "+particles.names[row]
                                f_input.close()
                                return (status, message)
                        if (value < 1): 
                                status  = 4
                                message = "line "+str(line)+", positive value expected for number of excitation types "+\
                                          "of "+particles.names[row]
                                f_input.close()
                                return (status, message)
                        if (value > MAX_EXCITATION_TYPES): 
                                status  = 6
                                message = "line "+str(line)+", number of excitation types exceeds maximum allowed"+\
                                          " (" + str(MAX_EXCITATION_TYPES) + ')'
                                f_input.close()
                                return (status, message)                        
                        if (len(l) < value + 7): 
                                status  = 2
                                message = "line "+str(line)+", missing excitation energy values for "+\
                                          particles.names[row]
                                f_input.close()
                                return (status, message)
                        particles.excitation_types[row] = abs(value)
                        for i in range(value):
                                try:
                                        value2 = float(l[7+i].strip())
                                except ValueError:
                                        status  = 3
                                        message = "line " + str(line) + ", numerical value expected"+\
                                                  " for excitation energy of "+\
                                                  particles.names[row]
                                        f_input.close()
                                        return (status, message)
                                if (value2 <= 0):
                                        status  = 4
                                        message = "line " + str(line) + ", positive value expected"+\
                                                  " for excitation energy of "+\
                                                  particles.names[row]
                                        f_input.close()
                                        return (status, message)                                                
                                particles.excitation_energy[row,i] = abs(value2)
                        i_diss = 7 + value  # First index for dissociation data (if any)
                        if debug: print('value=', value, 'i_diss=', i_diss)

                        # If the molecule is not monoatomic, get dissociation energies
                        if (particles.molecule_type[row] != 'a'):
                                try:
                                        value = int(l[i_diss].strip())
                                except IndexError:
                                        status  = 2
                                        message = "line "+str(line)+", missing values"
                                        f_input.close()
                                        return (status, message)
                                except ValueError:
                                        status  = 3
                                        message = "line "+str(line) + \
                                                  ", numerical value expected for number of dissociation types "+\
                                                  "of "+particles.names[row]
                                        f_input.close()
                                        return (status, message)
                                if (value < 1): 
                                        status  = 4
                                        message = "line "+str(line)+", positive value expected for number of dissociation types "+\
                                                  "of "+particles.names[row]
                                        f_input.close()
                                        return (status, message)
                                if (value > MAX_DISSOCIATION_TYPES): 
                                        status  = 6
                                        message = "line "+str(line)+", number of dissociation types exceeds maximum allowed"+\
                                                  " (" + str(MAX_DISSOCIATION_TYPES) + ')'
                                        f_input.close()
                                        return (status, message)                                
                                if (len(l) < value + i_diss + 1): 
                                        status  = 2
                                        message = "line "+str(line)+", missing dissociation energy values for "+\
                                                  particles.names[row]
                                        f_input.close()
                                        return (status, message)
                                particles.dissociation_types[row] = abs(value)
                                for i in range(value):
                                        try:
                                                value2 = float(l[i_diss+i+1].strip())
                                        except ValueError:
                                                status  = 3
                                                message = "line " + str(line) + ", numerical value expected"+\
                                                          " for dissociation energy of "+\
                                                          particles.names[row]
                                                f_input.close()
                                                return (status, message)
                                        if (value2 <= 0):
                                                status  = 4
                                                message = "line " + str(line) + ", positive value expected"+\
                                                          " for dissociation energy of "+\
                                                          particles.names[row]
                                                f_input.close()
                                                return (status, message)                                                
                                        particles.dissociation_energy[row,i] = abs(value2)
                        row = row + 1

                if (row >= particles.types): break

        f_input.close()

        return (status, message)
