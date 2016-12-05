"""Module containing the data for various IGM metal lines.
Species covered are:
H, He,C, N, O, Ne, Mg, Si, Fe
Reads the data from VPFIT's atom.dat file
Transitions are various.

"""

import re
import os.path

_lines_cache_ = None

class LineData(object):
    """Class to aggregate a number of lines from VPFIT tables. Reads from atom.dat,
    to get the results call get_line(species, ion), where ion is the transition number.
    At the moment only lowest level transition in each sequence is read."""
    def __init__(self, vpdat = None):
        self.species = ('H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe')
        #Put in the masses by hand for simplicity
        self.masses = {'H': 1.00794,'He': 4.002602,'C': 12.011,'N': 14.00674,'O': 15.9994,'Ne': 20.18,'Mg': 24.3050,'Si': 28.0855,'Fe': 55.847 }
        #9 empty lists, one list per species
        self.lines = {}
        if vpdat is None:
            global _lines_cache_
            vpdata = os.path.join(os.path.dirname(os.path.realpath(__file__)), "atom.dat")
            if _lines_cache_ is None:
                _lines_cache_ = read_vpfit(vpdata,self.species)
            self.lines = _lines_cache_
        else:
            self.lines = read_vpfit(vpdata,self.species)

    def __getitem__(self,specion):
        """Get data for a particular line.
        specion: a tuple of (specie, ion)
        specie: number of species
        ion: transition number (from 1)"""
        lines = self.lines[specion]
        return lines

    def __len__(self):
        return len(self.lines)

    def get_mass(self, specie):
        """Get the mass in amu for a species"""
        return self.masses[specie]

class Line(object):
    """Class to store the parameters of a single line.
    Data is:
        lambda_X - line width in angstroms (10^-10 m)
        fosc_X - oscillator strength (dimensionless)
        gamma_X - line's gamma in Hz"""
    def __init__(self, lambda_X, fosc_X, gamma_X):
        self.lambda_X = lambda_X
        self.fosc_X = fosc_X
        self.gamma_X = gamma_X

def read_vpfit(vpdat, species):
    """Read VPFIT's atom.dat file for the 9 species"""
    vp = open(vpdat)
    lines = {}
    inline = "--"
    while inline != "":
        inline = vp.readline()
        try:
            (specie, ion) = find_species(inline)
        except ValueError:
            continue
        #is this a line we care about?
        if species.count(specie) == 0 or ion < 0:
            continue
        (lambda_X, fosc_X, gamma_X) = parse_line_contents(inline)
        line = Line(lambda_X, fosc_X, gamma_X)
        #Only add the first transition for lines
        if (specie,ion) not in lines:
            lines[(specie,ion)]={int(lambda_X):line}
        else:
            lines[(specie,ion)][int(lambda_X)] = line
#                 print "Read line: ",specie,ion
    return lines

def find_species(line):
    """Parse a line to find which species it is: species is given by first two characters, unless the second character is a capital.
    There may be whitespace.
    Ionisation number is then a capital letter followed by three characters."""
    #Match a capital and possibly a lower case, followed by a capital followed by any three characters.
    mat = re.match(r"([A-Z]\s*[a-z]?)([A-Z].{3})",line)
    if mat is None:
        raise ValueError
    species = re.sub(r"\s","",mat.groups()[0])
    ion = mat.groups()[1]
    ion = re.sub(r"\s","",mat.groups()[1])
    ion_r = roman_to_int(ion)
    return (species, ion_r)

def parse_line_contents(line):
    """Extract the sigma, gamma, lambda and m from the line:
       just read until we get something looking like a float and then store the first three of them"""
    #Split by spaces
    spl = line.split()
    res = []
    for ss in spl:
        try:
            float(ss)
            res.append(float(ss))
        except ValueError:
            pass
        if len(res) == 3:
            break
    return res[0], res[1], res[2]



def roman_to_int(roman):
    """
    Convert a roman numeral to an integer.

    >>> r = range(1, 4000)
    >>> nums = [int_to_roman(i) for i in r]
    >>> ints = [roman_to_int(n) for n in nums]
    >>> print r == ints
    1

    >>> roman_to_int('VVVIV')
    Traceback (most recent call last):
     ...
    ValueError: roman is not a valid roman numeral: VVVIV
    >>> roman_to_int(1)
    Traceback (most recent call last):
     ...
    TypeError: expected string, got <type 'int'>
    >>> roman_to_int('a')
    Traceback (most recent call last):
     ...
    ValueError: roman is not a valid roman numeral: A
    >>> roman_to_int('IL')
    Traceback (most recent call last):
     ...
    ValueError: roman is not a valid roman numeral: IL

    Thanks Paul Winkler, on the internet.
    """
    if not isinstance(roman, str):
        raise TypeError("expected string, got %s" % type(roman))
    #roman = roman.upper()
    nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
    ints = [1000, 500, 100, 50,  10,  5,   1]
    places = []
    for c in roman:
        if not c in nums:
            raise ValueError("roman is not a valid roman numeral: %s" % roman)
    for i in range(len(roman)):
        c = roman[i]
        value = ints[nums.index(c)]
        # If the next place holds a larger number, this value is negative.
        try:
            nextvalue = ints[nums.index(roman[i +1])]
            if nextvalue > value:
                value *= -1
        except IndexError:
            # there is no next place.
            pass
        places.append(value)
    tot = 0
    for n in places:
        tot += n
    return tot
