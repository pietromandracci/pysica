# COPYRIGHT (c) 2020-2026 Pietro Mandracci

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

""" PYthon tools for SImulation and CAlculus: functions for dealing with time units.

    This module contains some functions which allow to manage time units.

    Documentation is also available in the docstrings.
"""

#+-------------------------+
#| Import required modules |
#+-------------------------+

import math

#+---------------------+
#| Numerical constants |
#+---------------------+

NS_IN_US   = int(1E3)
NS_IN_MS   = int(1E6)
NS_IN_S    = int(1E9)
NS_IN_MIN  = 60 * NS_IN_S
NS_IN_HOUR = 60 * NS_IN_MIN


#+-----------------------------------------------------------------+
#| Functions used to convert timestamps to readable time durations |
#+-----------------------------------------------------------------+
    

def timestamp2time(timestamp):
    """ Get a timestamp expressed in nanoseconds and converts it in hours + minutes + seconds + milliseconds + microseconds + nanoseconds

        Parameters
        ----------

        timestamp:  integer number expressing a time duration in nanoseconds
                    if it is not integer, it will be converted to integer

        Returns
        -------

        a tuple (hours, minutes, seconds, millisecons,  microsecond, nanoseconds)
    """

    t = int(timestamp)

    # Initialize time quantities
    hours        = 0
    minutes      = 0
    seconds      = 0
    milliseconds = 0
    microseconds = 0
    nanoseconds  = 0
    
    # Calculate the whole numbers of h, min, s, ms, us, ns
    if ( t >= NS_IN_HOUR):
       hours        = int(t /            NS_IN_HOUR)
       t            = t - hours        * NS_IN_HOUR
    if ( t >= NS_IN_MIN):
       minutes      = int(t /            NS_IN_MIN)
       t            = t - minutes      * NS_IN_MIN
    if ( t >= NS_IN_S):
       seconds      = int(t /            NS_IN_S)
       t            = t - seconds      * NS_IN_S
    if ( t >= NS_IN_MS):
       milliseconds = int(t /            NS_IN_MS)
       t            = t - milliseconds * NS_IN_MS       
    if ( t >= NS_IN_US):
       microseconds = int(t /            NS_IN_US)
       t            = t - microseconds * NS_IN_US
    nanoseconds = t

    return (hours, minutes, seconds, milliseconds, microseconds, nanoseconds)


def print_timestamp2time(timestamp, printzeros=False):
    """ Get a timestamp expressed in nanoseconds and prints it in hours + minutes + seconds + milliseconds + microseconds + nanoseconds

        Parameters
        ----------

        timestamp:  integer number expressing a time duration in nanoseconds
                    if it is not integer, it will be converted to integer
    
        printzeros: print also the zero quantities
                    e.g. '0 h 10 min 3 s 0 ms 100 us 25 ns' is returned instead of '10 min 3 s 100 us 25 ns'

        Returns
        -------

        a string espressing the time in hours, minutes, seconds, milliseconds, microseconds, nanoseconds
    """
    
    (hours, minutes, seconds, milliseconds, microseconds, nanoseconds) = timestamp2time(timestamp)
    
    s = ''
    if ( (hours        > 0) or printzeros ): s = s + str(hours).rjust(3)        + ' h '
    if ( (minutes      > 0) or printzeros ): s = s + str(minutes).rjust(3)      + ' min '
    if ( (seconds      > 0) or printzeros ): s = s + str(seconds).rjust(3)      + ' s '
    if ( (milliseconds > 0) or printzeros ): s = s + str(milliseconds).rjust(3) + ' ms '
    if ( (microseconds > 0) or printzeros ): s = s + str(microseconds).rjust(3) + ' us '
    if ( (nanoseconds  > 0) or printzeros ): s = s + str(nanoseconds).rjust(3)  + ' ns '

    return s    
