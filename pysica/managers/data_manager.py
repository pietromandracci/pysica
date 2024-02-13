# COPYRIGHT (c) 2020-2024 Pietro Mandracci

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

""" Tools to import/export numerical data from/to data files in ASCII format

    In this module some classes and methods are provided to deal with 
    ASCII data files of different types.

    * Configuration files, containing the values to assign 
      to some variables, in the form:
          # Comment
          variable1 = <value1> #comment
          <empty lines>
          variable2 = <value2> #comment
      empty lines and comments starting with '#' are allowed and ignored.

      Data are registered in a dictionary: 
          { 'variable1': [<value1>], 'variable2': [<value2>], ... }

    * Data sequence, containing a sequence of ordered values, 
      in the following format:
          value1  <separator> value2  <separator> value3 ... <separator> valueN

      Almost any type of ASCII character sequence can be used as separator, 
      but usually EOL, "," or TAB are more convenient.

      Space is *not* usable as separation character, since spaces are 
      removed before text parsing

    * Data grids, i.e. two dimensional matrix of data, expressed in ASCII format
          #Comment
          #Comment
          <Value11> <separator> <Value12> <separator> <Value13>   #Comment
          <Value21> <separator> <Value22> <separator> <Value23>
      empty lines and comments starting with '#' are allowed and ignored.

      Data are stored in a 2-d numpy array, just like the following
          [ [<Value11>, <Value12>, <Value13>], 
            [<Value21>, <Value22>, <Value23>] ]
      or the following 
          [ [<Value11>, <Value21>], 
            [<Value12>, <Value22>], 
            [<Value13>, <Value23>] ]

      You may also specify a number of lines to skip at the beginning 
      of the ASCII file: this is often useful in order to skip the header 
      of some data file formats.

    Documentation is also available in the docstrings.
"""

import numpy

from ..parameters import *

#+--------------------+
#| Configuration Data |
#+--------------------+

class ConfigurationData:
    """ Configuration data.
 
        This class defines a dictionary used to store variable names 
        and associated values: typical use is reading the values of 
        some variables from a configuration file.

        Format of the configuration file is like the following  
            #Comment
            variable1 = 123.45E-6
            variable2 = 7654.3E-2  #Comment
        empty lines and comments starting with '#' are allowed and ignored.
                
        The data are stored in a dictionary with the following form: 
            { <variable>: [list] }, where 
            * <variable> -> string describing a variable of some kind 
            * [list]     -> list contaning the variable value, 
                            and optionally other information
                            (see docstring of method read_file())
    """

    def __init__(self):
        """ Creates a dictionary for configuration data.

            Creates an empty dictionary, which will contain the couples 
            { <variable-name>: [list] }.
            To fill the dictionary with values read from a ascii file, 
            use the read_file() method.

            Initialized data attributes
            ---------------------------

            self.d:         an empty dictionary
        """

        # Dictionary in which variable names and corresponding
        # default values must be stored
        self.d = {}     

        
    def read_file(self, filename, check_values=False):
        """ Read configuration data from a file.

            This method reads variable names and values from 
            a ASCII file and stores them in the dictionary self.d 
            (which was created at the instance creation).

            It is also possible to ask the method to check 
            variable names and values and return an
            error if they do not meet some constraints, 
            by setting to True the parameter check_values.
            In this case, *before* calling the method, 
            you must fill the dictionary self.d 
            with couples of the type 
                { <variable1>: [list1], <variable2>: [list2], ... }
            where each list must have the format: 
                [ <value>, <string>, <code>, <value1>, <value2> ]
            * <value>:  numeric value associated to the variable, 
                        which will be overwritten if the ASCII file 
                        contains a different value for that variable
                        typically, this is the default value
            * <string>: a string describing the variable meaning
            * <code>:   a string describing the allowed range 
                        for <value>
                        'checkmin':            <value> >= <value1>
                        'checkmax':            <value> <= <value1>
                        'checkminstrict':      <value> >  <value1>
                        'checkmaxstrict':      <value> <  <value1>
                        'checkrange':          <value1> <= <value> <= <value2>
                        'checkrangestrict':    <value1> <  <value> <  <value2>
                        'checkrangestrictmin': <value1> <  <value> <= <value2>
                        'checkrangestrictmax': <value1> <= <value> <  <value2>
                        'boolean':             <value> = 0 or 1
            * <value1>: minimum, maximum, or lower boundary of allowed range 
            * <value2>: upper boundary of allowed range

            Example:
                particle = ConfigurationData()
                particle.dic = {        
                    "position":     [ 1.0, 'particle position' ] , 
                    "velocity":     [ 2.0, 'particle velocity', 
                                      'checkmin', 0.0 ], 

                    "energy":       [ 3.0, 'particle acceleration, ']
                    "mass":         [ 1.0, 'particle mass', 'checkrange', 
                                      1.0, 40.0 ],
                    "active":       [ 1, 'active particle', boolean ]
                                }
                particle.read_file("particledata.dat", check_values=True)


            Format of the ASCII file is like the following                     
                #Comment
                variable1 = 123.45E-6
                variable2 = 7654.3E-2  #Comment

            Parameters
            ----------

            filename:     a string containing the name of input data file
            check_values: if set to True, checks that variable names in file 
                          are already present in dictionary
                          and that their values are inside specified ranges

            Initialized data attributes
            ---------------------------

            self.d: dictionary containing variables and associated values

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = could not open file
                     2 = syntax error
                     3 = invalid variable assignment
                     4 = unknown variable
                     5 = bad range specification
                     6 = variable out of allowed range
                     7 = bad keyword in range specification
            message: a string containing an error message or 'Ok'
        """

        # These variables are used for error checking            
        status   = 0
        message  = OK

        try:
            f_input = open(filename, 'r')
        except IOError:
            status  = 1
            message = 'could not open file'
            return (status, message)

        # Read command lines from input file
        line = 0
        for string in f_input:
            line = line + 1
            # Cut newline charater from the end of command line
            string = string.strip(EOL)
            # Cut line portion at right of '#' char
            string = string.split('#')[0]
            # Cut spaces from beginning and end of command line
            string = string.strip()
            # Ignore line if it is empty
            if (len(string) != 0):    
                # Separate variable name from the value we want to assign to it
                l = string.split("=")   
                # If line doesn't contain "=" character, we have a syntax error
                if len(l) < 2:
                    status = 2
                    message = "line "+str(line)+": invalid syntax" 
                    break
                left_string  = l[0].strip()     # Variable name
                right_string = l[1].strip()     # Variable value

                # Try to transform right part of assignemnt to floating-point
                try:
                    value = float(right_string)
                except ValueError:
                    status = 3
                    message = ("line "+str(line)
                               + ": assignment to variable \""
                               + left_string + "\" is not a number")
                    break

                # Try to recognise variable name among dictionary
                if (check_values):      
                    check = False
                    for (v, l) in self.d.items():
                        if (left_string == v): 
                            self.d[left_string][0] = value
                            check = True
                            break                                   
                    if not check:
                        status = 4
                        message = ("line " + str(line)
                                   + ": unknown variable \"" + left_string
                                   + "\"")
                        break

                # If check is not required, just add the (variable: [value]) couple to dictionary
                else:
                    self.d[left_string] = [value]

                # Check if value is in the valid range
                # (if check was requested and range was specified)
                if (check_values & (len(self.d[left_string]) > 2)):
                    l = self.d[left_string]
                    e5 = ('variable \'' + left_string
                          + '\': bad specification of range boundaries')
                    e6 = l[1] + ' is out of allowed range '
                    e7 = ('variable \'' + left_string
                          + '\': unknown range specification keyword')
                    if (l[2] == 'checkmin'):
                        if (len(l) < 4):
                            status = 5
                            message = e5
                            break                                                   
                        elif (l[0] < l[3]): #value < min
                            status = 6
                            message = e6 + '( '+left_string+' >= '+str(l[3])+' )'
                            break
                    elif (l[2] == 'checkminstrict'):
                        if (len(l) < 4):
                            status = 5
                            message = e5
                            break                                                   
                        elif (l[0] <= l[3]): #value <= min
                            status = 6
                            message = e6 + '( '+left_string+' > '+str(l[3])+' )'
                            break
                    elif (l[2] == 'checkmax'):
                        if (len(l) < 4):
                            status = 5
                            message = e5
                            break                                                   
                        elif (l[0] > l[3]):                               #value > max
                            status = 6
                            message = e6 + '( '+left_string+' <= '+str(l[3])+' )'
                            break
                    elif (l[2] == 'checkmaxstrict'):
                        if (len(l)<4):
                            status = 5
                            message = e5
                            break                                                   
                        elif (l[0] >= l[3]):                      #value >= max
                            status = 6
                            message = e6 + '( '+left_string+' < '+str(l[3])+' )'
                            break
                    elif (l[2] == 'checkrange'):
                        if (len(l) < 5):
                            status = 5
                            message = e5
                            break                                                   
                        elif ((l[0] < l[3]) | (l[0] > l[4])):     #value<min or value>max
                            status = 6
                            message = (e6 + '( '+str(l[3]) + ' <= '
                                       + left_string + ' <= ' + str(l[4]) + ' )')
                            break                                   
                    elif (l[2]=='checkrangestrict'):
                        if (len(l) < 5):
                            status = 5
                            message = e5
                            break                                                   
                        elif ((l[0] <= l[3]) | (l[0] >= l[4])):   #value<=min or value>=max
                            status = 6
                            message = (e6 + '( ' + str(l[3]) + ' < '
                                       + left_string + ' < ' + str(l[4])+' )')
                            break
                    elif (l[2]=='checkrangestrictmin'):
                        if (len(l)<5):
                            status = 5
                            message = e5
                            break                                                   
                        elif ( (l[0]<=l[3]) | (l[0]>l[4]) ):    #value<=min or value>max
                            status = 6
                            message = (e6 + '( ' + str(l[3]) + ' < '
                                        + left_string + ' <= ' + str(l[4]) + ' )')
                            break
                    elif (l[2] == 'checkrangestrictmax'):             
                        if (len(l) < 5):
                            status = 5
                            message = e5
                            break                                                   
                        elif ((l[0] < l[3]) | (l[0] >= l[4])):    #value<min or value>=max
                            status = 6
                            message = (e6 + '( ' + str(l[3]) + ' <= '
                                       + left_string + ' < ' + str(l[4]) + ' )')
                            break
                    elif (l[2] == 'boolean'):
                        if (l[0] not in range(2)):
                            status = 6
                            message = e6 + '( '+left_string+' must be 0 or 1 )'
                            break   
                    else:
                        status = 7
                        message = e7 + ' (\''+l[2]+'\')'
                        break
        f_input.close()

        return (status, message)


    def write_file(self, filename):
        """ Write configuration data to a file.

            This method creates an ASCII file with the name given and writes 
            to it variable names and values from the dictionary self.d
            File format is like the following
            variable1 = 123.45E-6

            variable2 = 7654.3E-2

            Parameters
            ----------

            filename:       a string containing the name of input data file

            Used data attributes
            --------------------

            self.d:         dictionary containing variables and associated values

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = could not open file
                     2 = error during file writing (this *should* never happen...)
            message: a string containing an error message or 'Ok'
        """

        status  = 0
        message = OK    
        try:
            f_output = open(filename,'w')
        except IOError:
            status  = 1
            message = 'could not open file'
            return (status, message)

        # Write names and values to file
        for l in range(len(self.d)):
            item_list = list(self.d.items())
            variable = item_list[l][0]
            value    = item_list[l][1][0]
            comment  = '#' + item_list[l][1][1] + EOL
            if   (value == True):  value = 1 
            elif (value == False): value = 0
            string = variable + " = " + str(value) + EOL + EOL
            try:
                f_output.write(comment)
                f_output.write(string)
            except IOError:
                status = 2
                message = "line "+str(line)+": error while writing line"
                break

        f_output.close()

        return (status, message)


#+---------------+
#| Data sequence |
#+---------------+

class DataSequence:
    """ A monodimensonal sequence of values, to be read from a ASCII file """

    def read_file(self, filename, sep=EOL):
        """ Read a monodimensonal sequence of values from a ASCII file.

            filename:               a string containing the name of the file
            sep:                    character used as separator between values (usually the end-of-line character)

            Initialized data attributes
            ---------------------------

            self.data_array:        array of values
            self.n_data:            number of data in the array

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = could not open file
                     2 = error reading file
                     3 = no data
            message: a string containing an error message or 'Ok'
        """

        status  = 0                             # This variable is used for error checking
        message = OK

        try:
            f_input = open(filename,'r')
        except IOError:
            status  = 1
            message = 'could not open file: "' + filename + '"'
            return (status, message)

        # Read data from input file
        try:
            self.data_array = numpy.fromfile(file=f_input, dtype=float, count=-1, sep=sep)
        except (IOError or ValueError):
            status  = 2
            message = 'error reading file: "' + filename + '"'
            return (status, message)                        
        self.n_data = len(self.data_array)
        if (self.n_data==0):
            status  = 3
            message = 'no valid data retrieved'
            return (status, message)

        return (status, message)


    def write_file(self, filename, sep=EOL):
        """ Write a monodimensional numpy array to a ASCII file.

                Parameters
                ----------

                filename:               a string containing the name of the file

                sep:                    character used as separator between values (usually the end-of-line character)

                Used data attributes
                --------------------

                self.data_array:        array of values

                Returns
                -------

                (status, message)
                status:  0 = no error
                         1 = could not open file
                         2 = error during file writing (this *should* never happen...)
                message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        # Open file for writing 
        try:
            f_output = open(filename,'w')
        except IOError:
            status  = 1
            message = 'could not open file: "' + filename + '"'
            return (status, message)

        self.data_array.tofile(file=f_output, sep=sep, format="%s")

        f_output.close()

        return (status, message)

       

#+-----------+
#| Data grid |
#+-----------+

class DataGrid:
    """ A 2-dimensional array of data.

        This class defines a 2-dimensional numpy into which a numerical grid (two dimensional matrix) 
        can be read from a ascii file.
    """


    def read_file(self, filename, n_rows=0, n_columns=0,
                  sep=SEP, pad_value=0, transpose=False, skip=0, debug=False):
        """ Read a 2-d array of values from file.

            Read a table of values from a ASCII file into a numpy bidimensional array
            Empty lines, as well as lines having '#' as first non-blank character are ignored

            filename:   a string containing the name of the ASCII file
            n_rows:     number of rows (0 = autodetect)
                        if it exceeds the number of rows in the data, 
                        the array is padded with the value defined by the pad_value parameter
            n_columns:  number of columns (0 = autodetect)
                        if it exceeds the number of columns in the data, 
                        the array is padded with with the value defined by the pad_value parameter
            sep:        character used as separator between columns of the table 
                        (usually '\t' or ',')
            pad_value:  value used to fill the empty cells of the grid
            transpose:  if set to 
                            False: the array will be of the form [row, column] 
                            True:  the array will be of the form [column, row]
            skip:       number of lines to skip at file beginning
                        useful to skip the header of some kind of data file formats
            debug:      if set to True, debug information is printed

            Initialized data attributes
            ---------------------------

            self.data_array:        2-d numpy array

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = could not open file
                     2 = syntax error
                     3 = non-numerical value
                     4 = invalid number of rows/columns
                     5 = invalid number of rows to skip
            message: a string containing an error message or 'Ok'
        """

        status, message  = 0, OK                # This variable is used for error checking      

        # Check that the given numbers of rows ad columns are valid
        if (n_rows < 0) | (n_columns < 0):
            status  = 4
            message = 'invalid number of rows/columns'
            return (status, message)

        # Check that the number of rows to skip (if given) is positive
        if (skip < 0):
            status  = 5
            message = 'invalid number of rows to skip'
            return (status, message)

        # Open ASCII file
        try:
            f_input = open(filename,'r')
        except IOError:
            status  = 1
            message = 'could not open file'
            return (status, message)

        # If numbers of rows and/or columns are unknown, get them
        if (n_rows == 0) | (n_columns == 0):    # Get the number of rows and/or columns
            (rows, cols) = _get_rows_columns(f_input, sep, skip, debug)
            if debug: print("rows, cols = ", rows, cols)
            if (n_rows    == 0):    n_rows    = rows 
            if (n_columns == 0):    n_columns = cols                        

        # Read the data
        (status, message) = self._read_file(f_input, n_rows, n_columns,
                                            sep, pad_value, transpose, skip, debug)
        f_input.close()

        return (status, message)


    def _read_file(self, f_input, n_rows, n_columns, sep, pad_value, transpose,
                   skip, debug=False):
        """ Read a 2-d array of data from file (no check).

            Read a table of values from a ASCII file into a numpy bidimensional array
            WARNING: this method is not intended to be called directly: use read_file instead

            Parameters
            ----------

            filename:               a string containing the name of the file
            sep:                    character used as separator between columns of the table 
                                            (usually '\t' or ',')
            n_rows:                 number of rows
            n_columns:              number of columns
            transpose:              is set to 
                                            False: the array will be of the form [row, column] 
                                            True:  the array will be of the form [column, row]

            skip:                   number of lines to skip at file beginning
            debug:                  if True, debug information is printed

            Initialized data attributes
            ---------------------------

            self.data_array:        2-d numpy array

            Returns
            -------

            (status, message)
            status:  0 = no error
                     2 = syntax error
                     3 = non-numerical value
            message: a string containing an error message or 'Ok'

        """

        status, message  = 0, OK

        # Initialize numpy array
        if transpose:
            self.data_array = pad_value * numpy.ones((n_columns, n_rows))
        else:           
            self.data_array = pad_value * numpy.ones((n_rows, n_columns))

        line = 0        # Line of the ASCII file
        row  = 0        # Row of the table
        for string in f_input:
            line = line + 1
            if debug: print("Line: " +str(line)+ " = \'" +string+ "\'")
            if (line >= skip+1):
                # Cut newline charater from the end of line
                string = string.strip(EOL)
                # Cut line portion at right of '#' char
                string = string.split('#')[0]
                # Cut spaces from beginning and end of line
                string = string.strip()

                # Ignore line if it is empty
                if (len(string) != 0):
                    # Get a list of the strings in the line
                    l = string.split(sep)
                    for column in range(n_columns):
                        if (column <= len(l)-1):
                            # Try to transform string value to floating-point
                            try:
                                value = float(l[column].strip())
                            except ValueError:
                                status  = 3
                                message = "line "+str(line)+": non-numerical value"
                                break
                            if transpose:
                                self.data_array[column, row] = value
                            else:
                                self.data_array[row, column] = value
                    row = row + 1
                if (row > n_rows-1): break

        return (status, message)


    def write_file(self, filename, sep=SEP, transpose=False):
        """ Write a bidimensional numpy array to a ASCII file.

            Parameters
            ----------

            filename:               a string containing the name of the file
            sep:                    character used as separator between the two columns of the table 
                                    (usually <TAB> or ',')
            transpose:              is set to 
                                    * False: the array is saved in the form [row, column] 
                                    * True:  the array is saved in the form [column, row]
                                    where column is in range (0,1)

            Used data attributes
            --------------------

            self.data_array:        2-dimensional array of values

            Returns
            -------

            (status, message)
            status:  0 = no error
                     1 = could not open file
                     2 = table inconsistency
                     3 = error during file writing (this *should* never happen...)
            message: a string containing an error message or 'Ok'
        """

        status, message = 0, OK 

        # Open file for writing 
        try:
            f_output = open(filename,'w')
        except IOError:
            status  = 1
            message = 'could not open file'
            return (status, message)

        if transpose:   (n_columns, n_rows) = self.data_array.shape
        else:           (n_rows, n_columns) = self.data_array.shape

        # Write table to file
        for row in range(n_rows):
            string = ""
            for col in range(n_columns):
                if transpose: string = string + str(self.data_array[col, row])
                else:         string = string + str(self.data_array[row, col])
                if (col < n_columns-1): string = string + sep
            string = string + EOL
            try:
                f_output.write(string)
            except IOError:
                status = 2
                message = "line "+str(line)+": error while writing"
                break

        f_output.close()

        return (status, message)

#+-----------+
#| Utilities |
#+-----------+

def _get_rows_columns(f_input, sep=SEP, skip=0, debug=False):
    """ Gets number of rows and colums from a ascii data-file.

        Gets the number of rows and colums of matrix saved in a ASCII data file
        Empty lines, as well as lines starting with the `#' character are not considered

        Parameters
        ----------

        f_input:        data file identifier
        sep:            character used as separator between data inside rows
                        (usually '\t' or ',')
        skip:           number of ASCII lines to skip at the file beginning
                        this may be useful to skip the header of some kind of data file formats
        

        Returns
        -------

        (rows, columns)
    """

    line = 0
    rows = 0
    cols = 0
    for string in f_input:
        line = line + 1
        if (line >= skip+1):
            # Cut newline charater from the end of line
            string = string.strip(EOL)
            # Cut line portion at right of '#' char
            string = string.split('#')[0]
            # Cut spaces from beginning and end of line
            string = string.strip()
            if len(string)!=0:  
                rows = rows+1
                l = len(string.split(sep))
                if (l > cols): cols = l
    f_input.seek(0)
    
    return (rows, cols)


