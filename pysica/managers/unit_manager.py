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

""" PYthon tools for SImulation and CAlculus: functions for the expression of quantities in SI units using prefixes.

    This module contains some functions which allow to manage the use prefixes when printing a quantity 
    in SI units, just like in the following situation:
    * you want to print a quantity expressed in a certain SI unit (e.g. an energy expressed in joules)
    * the quantity may vary among several orders of magnitude (e.g. from 1E-15 J to 1E6 J)
    * you want to express the quantity using SI prefixies (e.g. 1 fJ and 1 MJ instead of 1E-15 J and 1E6 J)
    * you may also want to limit the number of significant digits to be shown (e.g. 1.15 nJ instead of 
      1.149456E-9 J)

    The most straightforward use of this module is to use the function print_unit, 
    giving to it a numerical value and a string that contains the symbol of the SI unit 
    in which the value is expressed.  
    The function will return a string with the value expressed in the unit, using prefixies, and limiting the 
    number of significant digits if requested:
        print_unit(1.123E-3, 'J')    -> '1.123 mJ'
        print_unit(1123456, 'J')     -> '1.123456 MJ'
        print_unit(1123456, 'J', 2)  -> '1.1 MJ'

    As an option, you may also request a check that the string you provide is a SI unit.

    Documentation is also available in the docstrings.
"""

#+-------------------------+
#| Import required modules |
#+-------------------------+

import math


#+---------------------+
#| Numerical constants |
#+---------------------+

units = (
    # Fundamental SI units
    'm',
    'kg',
    's',
    'A',
    'K',
    'mol',
    'cd',
    'rad',
    # (Some) derived SI units
    'sr',
    'Hz',
    'N',
    'Pa',
    'J',
    'W',
    'C',
    'V',
    'F',
    'ohm',
    'S',
    'Wb',
    'T',
    'H',
    'lm',
    'lx',
    'Bq',
    'Gy',
    'Sv',
    'kat',
    # Some selected non-SI units, accepted for use with the SI
    'eV',           
    'l', 'L',       # liter
    't',            # ton
    'u',            # atomic mass unit
    'ua',           # astronomical unit
    # Some selected non-SI units, not accepted for use with SI (but included in this module)
    'torr',
    'Torr',
    'bar',
    'atm'
)

EMAX = 24

mul = {
    1: " ",
    3: "k",
    6: "M",
    9: "G",
    12: "T",
    15: "P",
    18: "E",
    21: "Z",
    24: "Y"
}

div = {
    1: " ",
    3: "m",
    6: "u",
    9: "n",
    12: "p",
    15: "f",
    18: "a",
    21: "z",
    24: "y"
}



#+---------------------------------------------------------+
#| Functions used to print values, uncertainties and units |
#+---------------------------------------------------------+

def print_exp(value, digits=0, align=False, debug=False):
    """ Returns a string expressing a number in exponential format 
        with a maximum number of significant digits in the mantissa

        This function returns a string representing a number in exponential format, with a maximum number 
        of significant digits in the mantissa
        If the given number is int, it will be converted to float
        Examples:
            print_exp(123456,    1) -> '1.0e5 '
            print_exp(-175836,   3) -> '-1.76e5'
            print_exp(0.123456,  4) -> '1.235e-1'
            print_exp(-0.912356, 1) -> '-9.0e-1'
            print_exp(312346149, 5) -> '3.1235e7'

        Parameters
        ----------

        value:  numerical value to be converted
        digits: maximum number of significant digits allowed for the mantissa
                0 or a negative value means keep all the digits
        align:  is set to True, add a leading space to integer number
                to mantain alignment to upper rows in columnar output

        Returns
        -------

        a string expressing the number in exponential format with the significant digits in the mantissa
    """

    value = float(value)

    # If the given number is nan, or inf, return the value as string
    if ( math.isnan(value) or math.isinf(value) ): return str(value)

    # If the number is negative apply the procedure to the positive part and then change sign
    if (value < 0):  return '-'+print_exp(-value, digits, debug)

    if (value == 0):  (mantissa, esp) = (0.0, 0)
    else:             (mantissa, esp) = _sep_mantissa_exp(value, digits=digits, debug=debug)

    ms=str(mantissa)
    es=str(esp)
    if (digits > 0): ms = ms.ljust(digits+1,'0')

    # If the mantissa string ends with '.0' cut the last two characters        
    if ( (ms[-1] == '0') and (ms[-2] == '.') ): ms = ms[0:-2]        

    # If the number is not decimal, add a leading space to preserve alignment
    if ( align and  ('.' not in ms)) : ms = ' ' + ms

    return ms + 'e' + es


def print_unit(value, unit, digits=0, kg_check=True, unit_check=False, align=False, debug=False):
    """ Returns a string expressing a physical quantity in SI (or compatible) units using prefixies.

        This function takes as arguments a numerical value and a string representig a measure unit
        It will convert the value using SI prefixs
        If the given number is of type int, it will be converted to float
        You can also specify how many significant digits must be provided in the result
        Examples:
            print_unit(1.23456E5,  'Hz')   -> '123.456 kHz'
            print_unit(1.23456E-2, 's', 3) -> '12.3 ms'
            print_unit(1.23456E-7, 'm', 3) -> '123.0 nm'

        WARNING: the function will (try to) treat nicely the case of 'kg', unless you disable this feature.
        However, you must give the unit correctly (e.g. you must provide 'kg' and *not* 'Kg')

        Parameters
        ----------

        value:          value of the physical quantity
        unit:           measure unit in which the physical quantity is expressed
        digits:         maximum number of digits to be shown (zero or negative number = no limit)
        kg_check:       if set to True, treat nicely the case of 'kg'
        unit_check:     if set to True, checks that the string given is included in the unit list
                        and returns an error message if it is not
        debug:          if set to True, print some additional information for debugging purposes

        Returns
        -------

        a string containing the same measure expressed with SI prefixs
    """

    if debug: print("*** function print_unit ***")    
    # If required, check that the given unit is included in the list
    if unit_check:
        if unit not in units: return('ERROR: \"'+unit+'\" is not a valid unit')

    # Try to treat nicely the ugly case of kg
    if kg_check & (unit == 'kg'):
        value = value * 1.0E3
        unit = 'g'

    # Do the conversion
    (v,u) = change_unit(value, unit, digits, debug)

    if debug:
        print("value  = ", v)
        print("unit   = ", u)

    # If nan or inf is returned
    if ( math.isnan(v) or math.isinf(v)):
        vs = str(v)
        us = str(u)
        if (digits > 0): vs = vs.rjust(digits+1)
        return vs + '  ' + us

    if (digits > 0):
        if (v > 0): vs = str(v).ljust(digits+1,'0')
        else:       vs = str(v).ljust(digits+2,'0')
    else:
        vs = str(v)
        
    us = str(u)

    if debug:
        print("vs     = ", vs)
        print("us     = ", us)

    # If the string value ends with '.0', and the final zero is not required, cut the last two characters
    if ( (vs[-1] == '0') and (vs[-2] == '.') and (len(vs) > digits+1) ): vs = vs[0:-2]

    # If the number is not decimal, add a leading space to preserve alignment
    if ( align and ('.' not in vs) ) : vs = ' ' + vs
        
    return vs + ' ' + us    


def print_uncertainty(value, uncertainty, digits=2, debug=False):
    """ Return a string expressing (value+/-uncertainty), properly fixing the number of significant digits.

        This function takes as arguments two numerical values representig a value and the associated uncertainty
        It will return a string representing the couple value +/- uncertainty, 
        fixing the numer of significant digits of the uncertainty as requested.

        Parameters
        ----------

        value:  numerical value to be converted
        digits: maximum number of significant digits allowed for the uncertainty
                0 or a negative value means keep all the digits
        align:  is set to True, add a leading space to integer number
                to mantain alignment to upper rows in columnar output

        Returns
        -------

        a string expressing the number in exponential format with the significant digits in the mantissa
    """

    (value, uncertainty, esp) = fix_uncertainty(value, uncertainty, digits, debug)
    vs = str(value)
    us = str(uncertainty)
    es = str(esp)
    # If the mantissa string ends with '.0' cut the last two characters        
    if ( (vs[-1] == '0') and (vs[-2] == '.') ): vs = vs[0:-2]
    if ( (us[-1] == '0') and (us[-2] == '.') ): us = us[0:-2] 
    # If the number is not decimal, add a leading space to preserve alignment
    if ('.' not in vs): vs = ' ' + vs
    if ('.' not in us): us = ' ' + us        
    string = '(' + vs + '+/-'+ us + ')'
    if (esp != 0): string = string +  'e' + es

    return string


#+----------------------------------------------------------+
#| Functions used to modify values, uncertainties and units |
#+----------------------------------------------------------+

def change_unit(value, unit, digits=0, debug=False):
    """ Express a physical quantity in SI (or compatible) units using prefixies.
                
        This function takes as arguments a numerical value and a string representing a measure unit
        It will convert the value using SI prefixes
        If the given number is int, it will be converted to float
        You can also specify how many significant digits must be provided in the result
        Examples:
            change_unit(1.23456E5,  'Hz')   -> (123.456, 'kHz')
            change_unit(1.23456E-2, 's', 3) -> (12.3, 'ms')
            change_unit(1.23456E-7, 'm', 3) -> (123.0, 'nm')

        WARNING: no check is made that the given string is a valid unit; 
        as a consequence the "ugly" case of "kg" is *not* treated nicely
        Examples:
            change_unit(1.23456E5,  'xy')   -> (0.123456, 'Mxy')

        Parameters
        ----------

        value:          value of the physical quantity
        unit:           measure unit in which the physical quantity is expressed
        digits:         maximum number of digits to be shown (zero or negative number = no limit)
        debug:          if set to True, print some additional information for debugging purposes

        Returns
        -------
        (new-value, new-unit)
        new-value:      the value of the same quantity expressed in the new unit
        new-unit:       the new unit, which is the old one multiplied by 10**(3*n)
                        where 3*n can be in range (-24, -21, -18, ..., 18, 21, 24)
    """

    if debug: print("*** function change_unit ***")
    
    value = float(value)

    # If the value is 0, inf, or nan, no conversion is needed
    if (math.isnan(value) or math.isinf(value) or (value==0.0)): return(value, unit)

    # If the value is negative, apply conversion to -value and than return -m
    if (value < 0):
        (m, prefix) = change_unit(-value, unit, digits=digits, debug=debug)
        return (-m, prefix)

    # Separate mantissa and exponent, taking care that mantissa is > 1
    (mantissa, esp) = _sep_mantissa_exp(value, digits=digits)

    if debug:
        print("mantissa = ", mantissa)
        print("exponent = ", esp)

    # Exponent is greater then the maximum value for which a prefix exist
    if (esp > EMAX):
        e = EMAX
        k = esp-e
        m = mantissa * 10**k
    # Exponent is lower then the minimum value for which a prefix exist
    elif (esp < -EMAX):
        e = -EMAX       
        k = esp-e
        m = mantissa * 10**k
    # Exponent is between minimum and maximum values: generate new mantissa and exponent
    else:
        k = esp % 3
        if (k==0):
            e = esp
            m = mantissa
        else:
            e = int(math.floor(esp / 3)) * 3
            m = mantissa * 10**k

    # If requested, reduce the number of decimal digits in the mantissa
    if (digits > 0):
        m = fix_digits(m, digits)
              
    if    (e > 0): prefix = mul[e]  + unit
    elif  (e < 0): prefix = div[-e] + unit
    else:          prefix = mul[1]  + unit

    # This is for debugging purposes
    if debug:
        print("m        = ", m)
        print("k        = ", k)
        print("e        = ", e)
        print("prefix   = ", prefix)

    if debug: print("*** end of function change_unit ***")

    return (m, prefix)
        

def fix_digits(value, digits, debug=False):
    """ Approximates a value to a given number of significant digits.
                
        This function approximate a real number to the number of significant digits specified
        If the given number is int, it will be converted to float
        Examples:
            fix_digits(1.23456,   1) ->  1.0
            fix_digits(-1.75836,  3) -> -1.76
            fix_digits(0.123456,  4) ->  0.1235
            fix_digits(-0.912356, 1) -> -1.0
            fix_digits(312346149, 5) ->  312350000.0

        Parameters
        ----------

        value:          numerical value to be converted
        digits:         number of signifincat digits required 
                        0 or a negative value means "do nothing at all"
        debug:          if set to True, print some additional information for debugging purposes

        Returns
        -------
                
        the nearest number with the requested significant digits
    """

    value = float(value)

    # If the number of digits requested is zero or less then zero, return the number without conversion
    if (digits <= 0):  return value
    # If the given number is zero, nan, or inf, do nothing at all
    if ( math.isnan(value) or math.isinf(value) or (value==0.0) ): return value
    # If the number is negative apply the procedure to the positive part and then change sign
    if (value < 0):  return -fix_digits(-value, digits, debug)
    # Apply the conversion
    (mantissa, esp) = _sep_mantissa_exp(value, digits=digits, debug=debug)

#    return mantissa*10**esp ###################################################################################
    return round(mantissa*10**esp, digits-esp)


def fix_uncertainty(value, uncertainty, digits=2, debug=False):
    """ Approximates a couple value to a given number of significant digits.
                
        Given a couple of real numbers, which represent the value and uncertainty of a physical quantity,
        this function approximate the uncertainty to the number of significant digits specified
        and the value to the  number of significant digits covered by the uncertainty.
        If the given number is int, it will be converted to float
        Examples:
            fix_uncertainty(1.23456,  0.034567, 1) ->   1.23,  0.03
            fix_uncertainty(-1.23456, 0.034567, 1) ->  -1.23,  0.03
            fix_uncertainty(-1.75836, 0.2345,   3) ->  -1.758, 0.234

        Parameters
        ----------

        value:          value of the physical quantity
        uncertainty:    uncertainty of the physical quantity, only the absolute value will be taken
        digits:         number of significant digits required for +the unceratinty 
                        0 or a negative value means "do nothing at all"
        debug:          if set to True, print some additional information for debugging purposes

        Returns
        -------
        the nearest number with the requested significant digits
    """

    value       = float(value)
    uncertainty = abs(float(uncertainty))
    digits      = int(digits)

    # If the number of digits requested is zero or less then zero, return the number without conversion
    if (digits <= 0):  return (value, uncertainty)
    # If the value or uncertaninty is zero, nan, or inf, do nothing at all
    if (math.isnan(value) or math.isinf(value) or math.isnan(uncertainty) or math.isinf(uncertainty)):
        return value, uncertainty
    # If the value is negative, apply the procedure to the positive part and then change sign
    if (value < 0):
        (v_mantissa, u_mantissa, v_esp) = fix_uncertainty(-value, uncertainty, digits, debug)
        return (-v_mantissa, u_mantissa, v_esp)
    # Fix the n
    (v_mantissa, v_esp) = _sep_mantissa_exp(value,       digits=None,   debug=debug)
    (u_mantissa, u_esp) = _sep_mantissa_exp(uncertainty, digits=digits, debug=debug)
    v_digits = v_esp - u_esp + digits
    if (v_mantissa >= 1): v_mantissa = round(v_mantissa, v_digits-1)
    else:                 v_mantissa = round(v_mantissa, v_digits)
    if (v_esp > u_esp):
        u_mantissa = u_mantissa * 10**(u_esp-v_esp)
        return (v_mantissa, u_mantissa, v_esp)
    else:
        v_mantissa = v_mantissa * 10**(v_esp-u_esp)
        return (v_mantissa, u_mantissa, u_esp)

        
def _sep_mantissa_exp(value, digits=None, debug=False):
    """ Separate mantissa and exponent of a floating point number, 
        limiting the number of significant digits in the mantissa.

        The function returns mantissa and exponent of a given number
        If required, also limits the number of the significant digits in the mantissa.

        Parameters
        ----------

        value:          numerical value to be converted
        digits:         maximum number of significant digits allowed for the mantissa
                        0 or a negative value means keep all the digits
        debug:          if set to True, print some additional information for debugging purposes

        Returns
        -------

        (mantissa, esp)
        mantissa: the mantissa (float)
        esp:      the exponent (int)
    """

    value = float(value)
        
    # If the given number is zero, nan, or inf, do nothing at all
    if ( math.isnan(value) or math.isinf(value) or (value==0.0) ): return value

    # Separate mantissa and exponent, taking care that mantissa is > 1
    if (value<0):
        esp      = int(math.floor(math.log10(-value)))
    else: 
        esp      = int(math.floor(math.log10(value)))
    mantissa = value / 10**esp

    # This is for debugging purposes
    if debug: print(mantissa, esp)

    # Reduce the number of decimal digits in the mantissa
    if (not( (digits is None) or math.isnan(digits) or math.isinf(digits) ) and (digits > 0)):
        digits = int(digits)
        if (abs(mantissa) >= 1):
            mantissa = round(mantissa, digits-1)
        else:
            mantissa = round(mantissa, digits)
        
    return (mantissa,esp)
