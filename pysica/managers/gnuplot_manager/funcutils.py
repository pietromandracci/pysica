# COPYRIGHT (c) 2020-2022 Pietro Mandracci

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

""" Some utility functions. """

from .global_variables import *


def correct_filename(filename):
    """ Remove invalid characters from a filename.

        If the string contains characters in INVALID_CHARS, they are
        substituted with the char in SUBSTITUTE_CHAR.

        Parameters
        ----------

        filename: a string, representinge a filename

        Returns
        -------

        The modified string
    """

    for c in INVALID_CHARS:
        filename = filename.replace(c, SUBSTITUTE_CHAR)

    return filename


def data2string_1d(data, eol=EOL):
    """ Create a string with the given data, to be sent to gnuplot.

        The string can be used to pass data to gnuplot using the
        special filename "-". The form of the string is
        "data[0]\ndata[1]\n ... data[len(data)]\ne\n"

        Parameters
        ----------

        data:  data representing the points to plot
        eol:   end-of-line character

        Returns
        -------
 
        A string with the data in the form requested by gnuplot
    """
    #s = ''
    #for value in data:
    #    s += str(value) + eol
    #s += 'e' + eol
    #return s
    
    l = []
    for value in data:
        l.append(str(value))
    l.append('e')
    return eol.join(l) + eol


def data2string_2d(x_data, y_data, eol=EOL):
    """ Create a string with the given data, to be sent to gnuplot.

        The string can be used to pass data to gnuplot using the
        special filename "-". The form of the string is
        "x_data[0] y_data[0]\nx_data[1] y_data[1]\n ... 
         ... x_data[len(data)] y_data[len(data)]\ne\n"

        Parameters
        ----------

        x_data:   data representing x coordinates of the points to plot
        y_data:   data representing y coordinates of the points to plot
        eol:      end-of-line character

        Returns
        -------
 
        A string with the data in the form requested by gnuplot
    """
    #s = ''
    #for i in range(len(x_data)):
    #    s += str(x_data[i]) + ' ' + str(y_data[i]) + eol       
    #s += 'e' + eol    
    #return s

    l = []
    for i in range(len(x_data)):
        l.append( ' '.join( (str(x_data[i]), str(y_data[i])) ) )
    l.append('e')
    return eol.join(l) + eol
               

def data2string_3d(x_data, y_data, z_data, eol=EOL):
    """ Create a string with the given data, to be sent to gnuplot.

        The string can be used to pass data to gnuplot using the
        special filename "-". The form of the string is
        "x_data[0] y_data[0] z_data[0]\nx_data[1] y_data[1]  z_data[1]\n ... 
         ... x_data[len(data)] y_data[len(data)] z_data[len(data)]\ne\n"

        Parameters
        ----------

        x_data:   data representing x coordinates of the points to plot
        y_data:   data representing y coordinates of the points to plot
         z_data:   data representing z coordinates of the points to plot
        eol:      end-of-line character

        Returns
        -------
 
        A string with the data in the form requested by gnuplot
    """
    #s = ''
    #for i in range(len(x_data)):
    #    s += str(x_data[i]) + ' ' + str(y_data[i]) + ' ' + str(z_data[i]) + eol       
    #s += 'e' + eol
    #return s
    l = []
    for i in range(len(x_data)):
        l.append( ' '.join( (str(x_data[i]), str(y_data[i]), str(z_data[i])) ) )
    l.append('e')
    return eol.join(l) + eol
    

def data_file_1d(x_data, filename, eol=EOL):
    """ Create a file and write 1d data to it in csv format.

        The file will contain 1d data in ascii format, in a single column:

        x1
        x2
        ..

        Parameters
        ----------

        x_data:   data representing the points to plot
        filename: the name of the data file to be created
        eol:      end-of-line character
    """

    data_file = open(filename, 'w') 
    for x in x_data:
        data_file.write(str(x) + eol)
    data_file.close()


def data_file_2d(x_data, y_data, filename, sep=SEP, eol=EOL):
    """ Create a file and write 2d data to it in csv format.

        The file will contain data in ascii format, in the form:

        x1  y1
        x2  y2
        ..  ..

        Parameters
        ----------

        x_data:   data representing x coordinates of the points to plot
        y_data:   data representing y coordinates of the points to plot
        filename: the name of the data file to be created
        sep:      separator character between colums (vaues in a line)
        eol:      end-of-line character
    """

    data_file = open(filename, 'w') 
    for i in range(len(x_data)):
        data_file.write( str(x_data[i])
                         + sep
                         + str(y_data[i])
                         + eol )
    data_file.close()


def data_file_3d(x_data, y_data, z_data, filename, sep=SEP, eol=EOL):
    """ Create a file and write 3d data to it in csv format.

        The file will contain data in ascii format, in the form:

        x1  y1  z1
        x2  y2  z2
        ..  ..  ..

        Parameters
        ----------

        x_data:   data representing x coordinates of the points to plot
        y_data:   data representing y coordinates of the points to plot
        z_data:   data representing z coordinates of the points to plot
        filename: the name of the data file to be created
        sep:      separator character between colums (vaues in a line)
        eol:      end-of-line character
    """

    data_file = open(filename, 'w') 
    for i in range(len(x_data)):
        data_file.write(  str(x_data[i])
                          + sep
                          + str(y_data[i])
                          + sep
                          + str(z_data[i])                                
                          + eol )
    data_file.close()
