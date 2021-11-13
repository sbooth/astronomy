#!/usr/bin/env python3
#
#    MIT License
#
#    Copyright (c) 2019-2021 Don Cross <cosinekitty@gmail.com>
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#
"""Astronomy Engine by Don Cross

See the GitHub project page for full documentation, examples,
and other information:

https://github.com/cosinekitty/astronomy

"""

import math
import datetime
import enum
import re

KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit.
C_AUDAY   = 173.1446326846693       #<const> The speed of light expressed in astronomical units per day.

# Jupiter radius data are nominal values obtained from:
# https://www.iau.org/static/resolutions/IAU2015_English.pdf
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

JUPITER_EQUATORIAL_RADIUS_KM = 71492.0  #<const> The equatorial radius of Jupiter, expressed in kilometers.
JUPITER_POLAR_RADIUS_KM = 66854.0       #<const> The polar radius of Jupiter, expressed in kilometers.
JUPITER_MEAN_RADIUS_KM = 69911.0        #<const> The volumetric mean radius of Jupiter, expressed in kilometers.

# The radii of Jupiter's 4 largest moons were obtained from:
# https://ssd.jpl.nasa.gov/?sat_phys_par
IO_RADIUS_KM       = 1821.6     #<const> The mean radius of Jupiter's moon Io, expressed in kilometers.
EUROPA_RADIUS_KM   = 1560.8     #<const> The mean radius of Jupiter's moon Europa, expressed in kilometers.
GANYMEDE_RADIUS_KM = 2631.2     #<const> The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
CALLISTO_RADIUS_KM = 2410.3     #<const> The mean radius of Jupiter's moon Callisto, expressed in kilometers.

_CalcMoonCount = 0

_RAD2HOUR  =  3.819718634205488         # 12/pi = factor to convert radians to sidereal hours
_HOUR2RAD  =  0.2617993877991494365     # pi/12 = factor to convert sidereal hours to radians
_DAYS_PER_TROPICAL_YEAR = 365.24217
_PI2 = 2.0 * math.pi
_EPOCH = datetime.datetime(2000, 1, 1, 12)
_ASEC360 = 1296000.0
_ASEC2RAD = 4.848136811095359935899141e-6
_ARC = 3600.0 * 180.0 / math.pi     # arcseconds per radian
_METERS_PER_AU = KM_PER_AU * 1000.0
_ANGVEL = 7.2921150e-5
_SECONDS_PER_DAY = 24.0 * 3600.0
_SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592
_MEAN_SYNODIC_MONTH = 29.530588
_EARTH_ORBITAL_PERIOD = 365.256
_NEPTUNE_ORBITAL_PERIOD = 60189.0
_REFRACTION_NEAR_HORIZON = 34.0 / 60.0

_SUN_RADIUS_KM = 695700.0
_SUN_RADIUS_AU  = _SUN_RADIUS_KM / KM_PER_AU

_EARTH_FLATTENING = 0.996647180302104
_EARTH_EQUATORIAL_RADIUS_KM = 6378.1366
_EARTH_POLAR_RADIUS_KM = _EARTH_EQUATORIAL_RADIUS_KM * _EARTH_FLATTENING
_EARTH_EQUATORIAL_RADIUS_AU = _EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU
_EARTH_MEAN_RADIUS_KM = 6371.0      # mean radius of the Earth's geoid, without atmosphere
_EARTH_ATMOSPHERE_KM = 88.0         # effective atmosphere thickness for lunar eclipses
_EARTH_ECLIPSE_RADIUS_KM = _EARTH_MEAN_RADIUS_KM + _EARTH_ATMOSPHERE_KM

_MOON_EQUATORIAL_RADIUS_KM = 1738.1
_MOON_MEAN_RADIUS_KM       = 1737.4
_MOON_POLAR_RADIUS_KM      = 1736.0
_MOON_EQUATORIAL_RADIUS_AU = (_MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU)

_ASEC180 = 180.0 * 60.0 * 60.0
_AU_PER_PARSEC = _ASEC180 / math.pi
_EARTH_MOON_MASS_RATIO = 81.30056

#
#    Masses of the Sun and outer planets, used for:
#    (1) Calculating the Solar System Barycenter
#    (2) Integrating the movement of Pluto
#
#    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
#
#    Page 10 in the above document describes the constants used in the DE405 ephemeris.
#    The following are G*M values (gravity constant * mass) in [au^3 / day^2].
#    This side-steps issues of not knowing the exact values of G and masses M[i];
#    the products GM[i] are known extremely accurately.
#
_SUN_GM     = 0.2959122082855911e-03
_JUPITER_GM = 0.2825345909524226e-06
_SATURN_GM  = 0.8459715185680659e-07
_URANUS_GM  = 0.1292024916781969e-07
_NEPTUNE_GM = 0.1524358900784276e-07

@enum.unique
class _PrecessDir(enum.Enum):
    From2000 = 0
    Into2000 = 1

def _LongitudeOffset(diff):
    offset = diff
    while offset <= -180.0:
        offset += 360.0
    while offset > 180.0:
        offset -= 360.0
    return offset

def _NormalizeLongitude(lon):
    while lon < 0.0:
        lon += 360.0
    while lon >= 360.0:
        lon -= 360.0
    return lon

class Vector:
    """A Cartesian vector with 3 space coordinates and 1 time coordinate.

    The vector's space coordinates are measured in astronomical units (AU).
    The coordinate system varies and depends on context.
    The vector also includes a time stamp.

    Attributes
    ----------
    x : float
        The x-coordinate of the vector, measured in AU.
    y : float
        The y-coordinate of the vector, measured in AU.
    z : float
        The z-coordinate of the vector, measured in AU.
    t : Time
        The date and time at which the coordinate is valid.
    """
    def __init__(self, x, y, z, t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t

    def __repr__(self):
        return 'Vector({}, {}, {}, {})'.format(self.x, self.y, self.z, str(self.t))

    def __str__(self):
        return '({}, {}, {}, {})'.format(self.x, self.y, self.z, str(self.t))

    def Length(self):
        """Returns the length of the vector in AU."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z, self.t)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z, self.t)

class StateVector:
    """A combination of a position vector, a velocity vector, and a time.

    The position (x, y, z) is measured in astronomical units (AU).
    The velocity (vx, vy, vz) is measured in AU/day.
    The coordinate system varies and depends on context.
    The state vector also includes a time stamp.

    Attributes
    ----------
    x : float
        The x-coordinate of the position, measured in AU.
    y : float
        The y-coordinate of the position, measured in AU.
    z : float
        The z-coordinate of the position, measured in AU.
    vx : float
        The x-component of the velocity, measured in AU/day.
    vy : float
        The y-component of the velocity, measured in AU/day.
    vz : float
        The z-component of the velocity, measured in AU/day.
    t : Time
        The date and time at which the position and velocity vectors are valid.
    """
    def __init__(self, x, y, z, vx, vy, vz, t):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.t = t

    def __repr__(self):
        return 'StateVector[pos=({}, {}, {}), vel=({}, {}, {}), t={}]'.format(
            self.x, self.y, self.z,
            self.vx, self.vy, self.vz,
            str(self.t))

    def __str__(self):
        return '({}, {}, {}, {}, {}, {}, {})'.format(self.x, self.y, self.z, self.vx, self.vy, self.vz, str(self.t))


@enum.unique
class Body(enum.Enum):
    """The celestial bodies supported by Astronomy Engine calculations.

    Values
    ------
    Invalid: An unknown, invalid, or undefined celestial body.
    Mercury: The planet Mercury.
    Venus: The planet Venus.
    Earth: The planet Earth.
    Mars: The planet Mars.
    Jupiter: The planet Jupiter.
    Saturn: The planet Saturn.
    Uranus: The planet Uranus.
    Neptune: The planet Neptune.
    Pluto: The planet Pluto.
    Sun: The Sun.
    Moon: The Earth's moon.
    EMB: The Earth/Moon Barycenter.
    SSB: The Solar System Barycenter.
    """
    Invalid = -1
    Mercury = 0
    Venus = 1
    Earth = 2
    Mars = 3
    Jupiter = 4
    Saturn = 5
    Uranus = 6
    Neptune = 7
    Pluto = 8
    Sun = 9
    Moon = 10
    EMB = 11
    SSB = 12

def BodyCode(name):
    """Finds the Body enumeration value, given the name of a body.

    Parameters
    ----------
    name: str
        The common English name of a supported celestial body.

    Returns
    -------
    Body
        If `name` is a valid body name, returns the enumeration
        value associated with that body.
        Otherwise, returns `Body.Invalid`.

    Example
    -------

    >>> astronomy.BodyCode('Mars')
    <Body.Mars: 3>

    """
    if name not in Body.__members__:
        return Body.Invalid
    return Body[name]

def _IsSuperiorPlanet(body):
    return body in [Body.Mars, Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto]

_PlanetOrbitalPeriod = [
    87.969,
    224.701,
    _EARTH_ORBITAL_PERIOD,
    686.980,
    4332.589,
    10759.22,
    30685.4,
    _NEPTUNE_ORBITAL_PERIOD,
    90560.0
]

class Error(Exception):
    """Indicates an error in an astronomical calculation."""
    def __init__(self, message):
        Exception.__init__(self, message)

class DateTimeFormatError(Error):
    """The syntax of a UTC date/time string was not valid, or it contains invalid values."""
    def __init__(self, text):
        Error.__init__(self, 'The date/time string is not valid: "{}"'.format(text))

class EarthNotAllowedError(Error):
    """The Earth is not allowed as the celestial body in this calculation."""
    def __init__(self):
        Error.__init__(self, 'The Earth is not allowed as the body.')

class InvalidBodyError(Error):
    """The celestial body is not allowed for this calculation."""
    def __init__(self):
        Error.__init__(self, 'Invalid astronomical body.')

class BadVectorError(Error):
    """A vector magnitude is too small to have a direction in space."""
    def __init__(self):
        Error.__init__(self, 'Vector is too small to have a direction.')

class InternalError(Error):
    """An internal error occured that should be reported as a bug.

    Indicates an unexpected and unrecoverable condition occurred.
    If you encounter this error using Astronomy Engine, it would be very
    helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
    page on GitHub. Please include a copy of the stack trace, along with a description
    of how to reproduce the error. This will help improve the quality of
    Astronomy Engine for everyone! (Thank you in advance from the author.)
    """
    def __init__(self):
        Error.__init__(self, 'Internal error - please report issue at https://github.com/cosinekitty/astronomy/issues')

class NoConvergeError(Error):
    """A numeric solver did not converge.

    Indicates that there was a failure of a numeric solver to converge.
    If you encounter this error using Astronomy Engine, it would be very
    helpful to report it at the [Issues](https://github.com/cosinekitty/astronomy/issues)
    page on GitHub. Please include a copy of the stack trace, along with a description
    of how to reproduce the error. This will help improve the quality of
    Astronomy Engine for everyone! (Thank you in advance from the author.)
    """
    def __init__(self):
        Error.__init__(self, 'Numeric solver did not converge - please report issue at https://github.com/cosinekitty/astronomy/issues')

def _SynodicPeriod(body):
    if body == Body.Earth:
        raise EarthNotAllowedError()
    if body.value < 0 or body.value >= len(_PlanetOrbitalPeriod):
        raise InvalidBodyError()
    if body == Body.Moon:
        return _MEAN_SYNODIC_MONTH
    return abs(_EARTH_ORBITAL_PERIOD / (_EARTH_ORBITAL_PERIOD/_PlanetOrbitalPeriod[body.value] - 1.0))

def AngleBetween(a, b):
    """Calculates the angle in degrees between two vectors.

    Given a pair of vectors, this function returns the angle in degrees
    between the two vectors in 3D space.
    The angle is measured in the plane that contains both vectors.

    Parameters
    ----------
    a : Vector
        The first of a pair of vectors between which to measure an angle.
    b : Vector
        The second of a pair of vectors between which to measure an angle.

    Returns
    -------
    float
        The angle between the two vectors expressed in degrees.
        The value is in the range [0, 180].
    """
    r = a.Length() * b.Length()
    if r < 1.0e-8:
        return BadVectorError()
    dot = (a.x*b.x + a.y*b.y + a.z*b.z) / r
    if dot <= -1.0:
        return 180.0
    if dot >= +1.0:
        return 0.0
    return math.degrees(math.acos(dot))


def DeltaT_EspenakMeeus(ut):
    """The default Delta T function used by Astronomy Engine.

    Espenak and Meeus use a series of piecewise polynomials to
    approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
    See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    This is the default Delta T function used by Astronomy Engine.

    Parameters
    ----------
    ut: float
        The floating point number of days since noon UTC on January 1, 2000.

    Returns
    -------
    float
        The estimated difference TT-UT on the given date, expressed in seconds.
    """
    # Fred Espenak writes about Delta-T generically here:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html
    # He provides polynomial approximations for distant years here:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    # They start with a year value 'y' such that y=2000 corresponds
    # to the UTC Date 15-January-2000. Convert difference in days
    # to mean tropical years.

    y = 2000 + ((ut - 14) / _DAYS_PER_TROPICAL_YEAR)

    if y < -500:
        u = (y - 1820) / 100
        return -20 + (32 * u*u)

    if y < 500:
        u = y / 100
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3
        return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6

    if y < 1600:
        u = (y - 1000) / 100
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3
        return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6

    if y < 1700:
        u = y - 1600
        u2 = u*u; u3 = u*u2
        return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0

    if y < 1800:
        u = y - 1700
        u2 = u*u; u3 = u*u2; u4 = u2*u2
        return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000

    if y < 1860:
        u = y - 1800
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3; u7 = u3*u4
        return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7

    if y < 1900:
        u = y - 1860
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3
        return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174

    if y < 1920:
        u = y - 1900
        u2 = u*u; u3 = u*u2; u4 = u2*u2
        return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4

    if y < 1941:
        u = y - 1920
        u2 = u*u; u3 = u*u2
        return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3

    if y < 1961:
        u = y - 1950
        u2 = u*u; u3 = u*u2
        return 29.07 + 0.407*u - u2/233 + u3/2547

    if y < 1986:
        u = y - 1975
        u2 = u*u; u3 = u*u2
        return 45.45 + 1.067*u - u2/260 - u3/718

    if y < 2005:
        u = y - 2000
        u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3
        return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5

    if y < 2050:
        u = y - 2000
        return 62.92 + 0.32217*u + 0.005589*u*u

    if y < 2150:
        u = (y-1820)/100
        return -20 + 32*u*u - 0.5628*(2150 - y)

    # all years after 2150
    u = (y - 1820) / 100
    return -20 + (32 * u*u)


_DeltaT = DeltaT_EspenakMeeus


def _TerrestrialTime(ut):
    return ut + _DeltaT(ut) / 86400.0

def _UniversalTime(tt):
    # This is the inverse function of _TerrestrialTime.
    # This is an iterative numerical solver, but because
    # the relationship between UT and TT is almost perfectly linear,
    # it converges extremely fast (never more than 3 iterations).
    dt = _TerrestrialTime(tt) - tt      # first approximation of dt = tt - ut
    while True:
        ut = tt - dt
        tt_check = _TerrestrialTime(ut)
        err = tt_check - tt
        if abs(err) < 1.0e-12:
            return ut
        dt += err

_TimeRegex = re.compile(r'^([0-9]{1,4})-([0-9]{2})-([0-9]{2})(T([0-9]{2}):([0-9]{2})(:([0-9]{2}(\.[0-9]+)?))?Z)?$')

class Time:
    """Represents a date and time used for performing astronomy calculations.

    All calculations performed by Astronomy Engine are based on
    dates and times represented by `Time` objects.

    Parameters
    ----------
    ut : float
        UT1/UTC number of days since noon on January 1, 2000.
        See the `ut` attribute of this class for more details.

    Attributes
    ----------
    ut : float
        The floating point number of days of Universal Time since noon UTC January 1, 2000.
        Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
        not exactly equivalent; UTC and UT1 can disagree by up to 0.9 seconds.
        This approximation is sufficient for the accuracy requirements of Astronomy Engine.
        Universal Time Coordinate (UTC) is the international standard for legal and civil
        timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
        UTC is kept in sync with unpredictable observed changes in the Earth's rotation
        by occasionally adding leap seconds as needed.
        UT1 is an idealized time scale based on observed rotation of the Earth, which
        gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
        large scale weather events like hurricanes, and internal seismic and convection effects.
        Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
        is adjusted by a scheduled whole number of leap seconds as needed.
        The value in `ut` is appropriate for any calculation involving the Earth's rotation,
        such as calculating rise/set times, culumination, and anything involving apparent
        sidereal time.
        Before the era of atomic timekeeping, days based on the Earth's rotation
        were often known as *mean solar days*.
    tt : float
        Terrestrial Time days since noon on January 1, 2000.
        Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
        In this system, days are not based on Earth rotations, but instead by
        the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
        divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
        for changes in the Earth's rotation.
        The value in `tt` is used for calculations of movements not involving the Earth's rotation,
        such as the orbits of planets around the Sun, or the Moon around the Earth.
        Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
    """
    def __init__(self, ut, tt = None):
        self.ut = ut
        if tt is None:
            self.tt = _TerrestrialTime(ut)
        else:
            self.tt = tt
        self.etilt = None

    @staticmethod
    def FromTerrestrialTime(tt):
        """Creates a #Time object from a Terrestrial Time day value.

        Parameters
        ----------
        tt : The number of days after the J2000 epoch.

        Returns
        -------
        Time
        """
        return Time(_UniversalTime(tt), tt)

    @staticmethod
    def Parse(text):
        """Creates a #Time object from a string of the form 'yyyy-mm-ddThh:mm:ss.sssZ'

        Parses a UTC date and time from a string and returns a #Time object.
        Permits a subset of ISO 8601 format.
        The year, month, and day are required.
        Hours, minutes, seconds, and fractions of a second are optional.
        If time is specified, there must be a 'T' between the date and the time
        and a 'Z' at the end of the time.

        Parameters
        ----------
        text : string
            A string of the following formats:
            `yyyy-mm-dd`
            `yyyy-mm-ddThh:mmZ`
            `yyyy-mm-ddThh:mm:ssZ`
            `yyyy-mm-ddThh:mm:ss.sssZ`

        Returns
        -------
        Time
        """
        m = _TimeRegex.match(text)
        if m is None:
            raise DateTimeFormatError(text)
        year = int(m.group(1))
        month = int(m.group(2))
        if not (1 <= month <= 12):
            raise DateTimeFormatError(text)
        day = int(m.group(3))
        if not (1 <= day <= 31):
            raise DateTimeFormatError(text)
        hour = int(m.group(5) or '0')
        if not (0 <= hour <= 23):
            raise DateTimeFormatError(text)
        minute = int(m.group(6) or '0')
        if not (0 <= minute <= 59):
            raise DateTimeFormatError(text)
        second = float(m.group(8) or '0')
        if not (0.0 <= second < 60.0):
            raise DateTimeFormatError(text)
        return Time.Make(year, month, day, hour, minute, second)

    @staticmethod
    def Make(year, month, day, hour, minute, second):
        """Creates a #Time object from a UTC calendar date and time.

        Parameters
        ----------
        year : int
            The UTC 4-digit year value, e.g. 2019.
        month : int
            The UTC month in the range 1..12.
        day : int
            The UTC day of the month, in the range 1..31.
        hour : int
            The UTC hour, in the range 0..23.
        minute : int
            The UTC minute, in the range 0..59.
        second : float
            The real-valued UTC second, in the range [0, 60).

        Returns
        -------
        Time
        """
        micro = round(math.fmod(second, 1.0) * 1000000)
        second = math.floor(second - micro/1000000)
        d = datetime.datetime(year, month, day, hour, minute, second, micro)
        ut = (d - _EPOCH).total_seconds() / 86400
        return Time(ut)

    @staticmethod
    def Now():
        """Returns the computer's current date and time in the form of a #Time object.

        Uses the computer's system clock to find the current UTC date and time.
        Converts that date and time to a #Time value and returns the result.
        Callers can pass this value to other Astronomy Engine functions to
        calculate current observational conditions.

        Returns
        -------
        Time
        """
        ut = (datetime.datetime.utcnow() - _EPOCH).total_seconds() / 86400.0
        return Time(ut)

    def AddDays(self, days):
        """Calculates the sum or difference of a #Time with a specified real-valued number of days.

        Sometimes we need to adjust a given #Time value by a certain amount of time.
        This function adds the given real number of days in `days` to the date and time
        in the calling object.

        More precisely, the result's Universal Time field `ut` is exactly adjusted by `days`
        and the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
        using a best-fit piecewise polynomial model devised by
        [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).

        The value of the calling object is not modified. This function creates a brand new
        #Time object and returns it.

        Parameters
        ----------
        days : float
            A floating point number of days by which to adjust `time`.
            May be negative, 0, or positive.

        Returns
        -------
        Time
        """
        return Time(self.ut + days)

    def __repr__(self):
        return 'Time(' + str(self) + ')'

    def __str__(self):
        millis = round(self.ut * 86400000.0)
        n = _EPOCH + datetime.timedelta(milliseconds=millis)
        return '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}Z'.format(n.year, n.month, n.day, n.hour, n.minute, n.second, math.floor(n.microsecond / 1000))

    def Utc(self):
        """Returns the UTC date and time as a `datetime` object.

        Uses the standard [`datetime`](https://docs.python.org/3/library/datetime.html) class
        to represent the date and time in this Time object.

        Returns
        -------
        datetime
        """
        return _EPOCH + datetime.timedelta(days=self.ut)

    def _etilt(self):
        # Calculates precession and nutation of the Earth's axis.
        # The calculations are very expensive, so lazy-evaluate and cache
        # the result inside this Time object.
        if self.etilt is None:
            self.etilt = _e_tilt(self)
        return self.etilt

    def __lt__(self, other):
        return self.tt < other.tt

    def __eq__(self, other):
        return self.tt == other.tt

    def __le__(self, other):
        return self.tt <= other.tt

    def __ne__(self, other):
        return self.tt != other.tt

    def __gt__(self, other):
        return self.tt > other.tt

    def __ge__(self, other):
        return self.tt >= other.tt


class Observer:
    """Represents the geographic location of an observer on the surface of the Earth.

    Parameters
    ----------
    latitude : float
        Geographic latitude in degrees north of the equator.
    longitude : float
        Geographic longitude in degrees east of the prime meridian at Greenwich, England.
    height : float
        Elevation above sea level in meters.
    """
    def __init__(self, latitude, longitude, height=0.0):
        self.latitude = latitude
        self.longitude = longitude
        self.height = height

    def __repr__(self):
        return 'Observer({}, {}, {})'.format(self.latitude, self.longitude, self.height)

    def __str__(self):
        text = '('
        text += 'S' if (self.latitude < 0) else 'N'
        text += '{}, '.format(abs(self.latitude))
        text += 'W' if (self.longitude < 0) else 'E'
        text += '{}, '.format(abs(self.longitude))
        text += '{}m'.format(self.height)
        text += ')'
        return text

class RotationMatrix:
    """Contains a rotation matrix that can be used to transform one
    coordinate system into another.

    Parameters
    ----------
    rot : float[3][3]
        A normalized 3x3 rotation matrix.
    """
    def __init__(self, rot):
        self.rot = rot

class Spherical:
    """Holds spherical coordinates: latitude, longitude, distance.

    Parameters
    ----------
    lat : float
        The latitude angle: -90..+90 degrees.
    lon : float
        The longitude angle: 0..360 degrees.
    dist : float
        Distance in AU.
    """
    def __init__(self, lat, lon, dist):
        self.lat = lat
        self.lon = lon
        self.dist = dist

class _iau2000b:
    def __init__(self, time):
        t = time.tt / 36525.0
        el  = math.fmod((485868.249036 + t*1717915923.2178), _ASEC360) * _ASEC2RAD
        elp = math.fmod((1287104.79305 + t*129596581.0481),  _ASEC360) * _ASEC2RAD
        f   = math.fmod((335779.526232 + t*1739527262.8478), _ASEC360) * _ASEC2RAD
        d   = math.fmod((1072260.70369 + t*1602961601.2090), _ASEC360) * _ASEC2RAD
        om  = math.fmod((450160.398036 - t*6962890.5431),    _ASEC360) * _ASEC2RAD
        dp = 0
        de = 0
$ASTRO_IAU_DATA()
        self.dpsi = -0.000135 + (dp * 1.0e-7)
        self.deps = +0.000388 + (de * 1.0e-7)

def _mean_obliq(tt):
    t = tt / 36525
    asec = (
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406
    )
    return asec / 3600.0

class _e_tilt:
    def __init__(self, time):
        e = _iau2000b(time)
        self.dpsi = e.dpsi
        self.deps = e.deps
        self.mobl = _mean_obliq(time.tt)
        self.tobl = self.mobl + (e.deps / 3600.0)
        self.tt = time.tt
        self.ee = e.dpsi * math.cos(math.radians(self.mobl)) / 15.0

def _ecl2equ_vec(time, ecl):
    obl = math.radians(_mean_obliq(time.tt))
    cos_obl = math.cos(obl)
    sin_obl = math.sin(obl)
    return [
        ecl[0],
        ecl[1]*cos_obl - ecl[2]*sin_obl,
        ecl[1]*sin_obl + ecl[2]*cos_obl
    ]

def _precession_rot(time, direction):
    eps0 = 84381.406
    t = time.tt / 36525

    psia  = (((((-    0.0000000951  * t
                 +    0.000132851 ) * t
                 -    0.00114045  ) * t
                 -    1.0790069   ) * t
                 + 5038.481507    ) * t)

    omegaa = (((((+   0.0000003337  * t
                 -    0.000000467 ) * t
                 -    0.00772503  ) * t
                 +    0.0512623   ) * t
                 -    0.025754    ) * t + eps0)

    chia  = (((((-    0.0000000560  * t
                 +    0.000170663 ) * t
                 -    0.00121197  ) * t
                 -    2.3814292   ) * t
                 +   10.556403    ) * t)

    eps0 *= _ASEC2RAD
    psia *= _ASEC2RAD
    omegaa *= _ASEC2RAD
    chia *= _ASEC2RAD

    sa = math.sin(eps0)
    ca = math.cos(eps0)
    sb = math.sin(-psia)
    cb = math.cos(-psia)
    sc = math.sin(-omegaa)
    cc = math.cos(-omegaa)
    sd = math.sin(chia)
    cd = math.cos(chia)

    xx =  cd * cb - sb * sd * cc
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc
    xy = -sd * cb - sb * cd * cc
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc
    xz =  sb * sc
    yz = -sc * cb * ca - sa * cc
    zz = -sc * cb * sa + cc * ca
    if direction == _PrecessDir.Into2000:
        # Perform rotation from other epoch to J2000.0.
        return RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ])

    if direction == _PrecessDir.From2000:
        # Perform rotation from J2000.0 to other epoch.
        return RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])

    raise Error('Inalid precession direction')

def _precession(pos, time, direction):
    r = _precession_rot(time, direction)
    return [
        r.rot[0][0]*pos[0] + r.rot[1][0]*pos[1] + r.rot[2][0]*pos[2],
        r.rot[0][1]*pos[0] + r.rot[1][1]*pos[1] + r.rot[2][1]*pos[2],
        r.rot[0][2]*pos[0] + r.rot[1][2]*pos[1] + r.rot[2][2]*pos[2]
    ]

class Equatorial:
    """Equatorial angular coordinates

    Coordinates of a celestial body as seen from the Earth.
    Can be geocentric or topocentric, depending on context.
    The coordinates are oriented with respect to the Earth's
    equator projected onto the sky.

    Attributes
    ----------
    ra : float
        Right ascension in sidereal hours.
    dec : float
        Declination in degrees.
    dist : float
        Distance to the celestial body in AU.
    vec : Vector
        The equatorial coordinates in cartesian form, using AU distance units.
        x = direction of the March equinox,
        y = direction of the June solstice,
        z = north.
    """
    def __init__(self, ra, dec, dist, vec):
        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.vec = vec


def _vector2radec(pos, time):
    xyproj = pos[0]*pos[0] + pos[1]*pos[1]
    dist = math.sqrt(xyproj + pos[2]*pos[2])
    if xyproj == 0.0:
        if pos[2] == 0.0:
            # Indeterminate coordinates: pos vector has zero length.
            raise BadVectorError()
        ra = 0.0
        if pos[2] < 0.0:
            dec = -90.0
        else:
            dec = +90.0
    else:
        ra = _RAD2HOUR * math.atan2(pos[1], pos[0])
        if ra < 0:
            ra += 24
        dec = math.degrees(math.atan2(pos[2], math.sqrt(xyproj)))
    vec = Vector(pos[0], pos[1], pos[2], time)
    return Equatorial(ra, dec, dist, vec)


def _nutation_rot(time, direction):
    tilt = time._etilt()
    oblm = math.radians(tilt.mobl)
    oblt = math.radians(tilt.tobl)
    psi = tilt.dpsi * _ASEC2RAD
    cobm = math.cos(oblm)
    sobm = math.sin(oblm)
    cobt = math.cos(oblt)
    sobt = math.sin(oblt)
    cpsi = math.cos(psi)
    spsi = math.sin(psi)

    xx = cpsi
    yx = -spsi * cobm
    zx = -spsi * sobm
    xy = spsi * cobt
    yy = cpsi * cobm * cobt + sobm * sobt
    zy = cpsi * sobm * cobt - cobm * sobt
    xz = spsi * sobt
    yz = cpsi * cobm * sobt - sobm * cobt
    zz = cpsi * sobm * sobt + cobm * cobt

    if direction == _PrecessDir.From2000:
        # convert J2000 to of-date
        return RotationMatrix([
            [xx, xy, xz],
            [yx, yy, yz],
            [zx, zy, zz]
        ])

    if direction == _PrecessDir.Into2000:
        # convert of-date to J2000
        return RotationMatrix([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz]
        ])

    raise Error('Invalid nutation direction')


def _nutation(pos, time, direction):
    r = _nutation_rot(time, direction)
    return [
        r.rot[0][0]*pos[0] + r.rot[1][0]*pos[1] + r.rot[2][0]*pos[2],
        r.rot[0][1]*pos[0] + r.rot[1][1]*pos[1] + r.rot[2][1]*pos[2],
        r.rot[0][2]*pos[0] + r.rot[1][2]*pos[1] + r.rot[2][2]*pos[2]
    ]

def _era(time):        # Earth Rotation Angle
    thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut
    thet3 = math.fmod(time.ut, 1.0)
    theta = 360.0 * math.fmod((thet1 + thet3), 1.0)
    if theta < 0.0:
        theta += 360.0
    return theta

def _sidereal_time(time):
    t = time.tt / 36525.0
    eqeq = 15.0 * time._etilt().ee    # Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
    theta = _era(time)
    st = (eqeq + 0.014506 +
        (((( -    0.0000000368   * t
            -    0.000029956  ) * t
            -    0.00000044   ) * t
            +    1.3915817    ) * t
            + 4612.156534     ) * t)
    gst = math.fmod((st/3600.0 + theta), 360.0) / 15.0
    if gst < 0.0:
        gst += 24.0
    # return sidereal hours in the half-open range [0, 24).
    return gst

def _inverse_terra(ovec, st):
    # Convert from AU to kilometers
    x = ovec[0] * KM_PER_AU
    y = ovec[1] * KM_PER_AU
    z = ovec[2] * KM_PER_AU
    p = math.sqrt(x*x + y*y)
    if p < 1.0e-6:
        # Special case: within 1 millimeter of a pole!
        # Use arbitrary longitude, and latitude determined by polarity of z.
        lon_deg = 0.0
        if z > 0.0:
            lat_deg = +90.0
        else:
            lat_deg = -90.0
        # Elevation is calculated directly from z
        height_km = abs(z) - _EARTH_POLAR_RADIUS_KM
    else:
        stlocl = math.atan2(y, x)
        # Calculate exact longitude.
        lon_deg = math.degrees(stlocl) - (15.0 * st)
        # Normalize longitude to the range (-180, +180].
        while lon_deg <= -180.0:
            lon_deg += 360.0
        while lon_deg > +180.0:
            lon_deg -= 360.0
        # Numerically solve for exact latitude, using Newton's Method.
        F = _EARTH_FLATTENING ** 2
        # Start with initial latitude estimate, based on a spherical Earth.
        lat = math.atan2(z, p)
        while True:
            # Calculate the error function W(lat).
            # We try to find the root of W, meaning where the error is 0.
            cos = math.cos(lat)
            sin = math.sin(lat)
            factor = (F-1)*_EARTH_EQUATORIAL_RADIUS_KM
            cos2 = cos*cos
            sin2 = sin*sin
            radicand = cos2 + F*sin2
            denom = math.sqrt(radicand)
            W = (factor*sin*cos)/denom - z*cos + p*sin
            if abs(W) < 1.0e-12:
                # The error is now negligible
                break
            # Error is still too large. Find the next estimate.
            # Calculate D = the derivative of W with respect to lat.
            D = factor*((cos2 - sin2)/denom - sin2*cos2*(F-1)/(factor*radicand)) + z*sin + p*cos
            lat -= W/D
        # We now have a solution for the latitude in radians.
        lat_deg = math.degrees(lat)
        # Solve for exact height in meters.
        # There are two formulas I can use. Use whichever has the less risky denominator.
        adjust = _EARTH_EQUATORIAL_RADIUS_KM / denom
        if abs(sin) > abs(cos):
            height_km = z/sin - F*adjust
        else:
            height_km = p/cos - adjust
    return Observer(lat_deg, lon_deg, 1000*height_km)

def _terra(observer, st):
    df2 = _EARTH_FLATTENING ** 2
    phi = math.radians(observer.latitude)
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    c = 1.0 / math.sqrt(cosphi*cosphi + df2*sinphi*sinphi)
    s = df2 * c
    ht_km = observer.height / 1000.0
    ach = _EARTH_EQUATORIAL_RADIUS_KM*c + ht_km
    ash = _EARTH_EQUATORIAL_RADIUS_KM*s + ht_km
    stlocl = math.radians(15.0*st + observer.longitude)
    sinst = math.sin(stlocl)
    cosst = math.cos(stlocl)
    return [
        ach * cosphi * cosst / KM_PER_AU,
        ach * cosphi * sinst / KM_PER_AU,
        ash * sinphi / KM_PER_AU
    ]

def _geo_pos(time, observer):
    gast = _sidereal_time(time)
    pos1 = _terra(observer, gast)
    pos2 = _nutation(pos1, time, _PrecessDir.Into2000)
    outpos = _precession(pos2, time, _PrecessDir.Into2000)
    return outpos

def _spin(angle, pos1):
    angr = math.radians(angle)
    cosang = math.cos(angr)
    sinang = math.sin(angr)
    return [
        +cosang*pos1[0] + sinang*pos1[1],
        -sinang*pos1[0] + cosang*pos1[1],
        pos1[2]
    ]

#----------------------------------------------------------------------------
# BEGIN CalcMoon

def _Array1(xmin, xmax):
    return dict((key, 0j) for key in range(xmin, 1+xmax))

def _Array2(xmin, xmax, ymin, ymax):
    return dict((key, _Array1(ymin, ymax)) for key in range(xmin, 1+xmax))

class _moonpos:
    def __init__(self, lon, lat, dist):
        self.geo_eclip_lon = lon
        self.geo_eclip_lat = lat
        self.distance_au = dist

def _CalcMoon(time):
    global _CalcMoonCount
    _CalcMoonCount += 1

    T = time.tt / 36525
    ex = _Array2(-6, 6, 1, 4)

    def Sine(phi):
        return math.sin(_PI2 * phi)

    def Frac(x):
        return x - math.floor(x)

    T2 = T*T
    DLAM = 0
    DS = 0
    GAM1C = 0
    SINPI = 3422.7000
    S1 = Sine(0.19833+0.05611*T)
    S2 = Sine(0.27869+0.04508*T)
    S3 = Sine(0.16827-0.36903*T)
    S4 = Sine(0.34734-5.37261*T)
    S5 = Sine(0.10498-5.37899*T)
    S6 = Sine(0.42681-0.41855*T)
    S7 = Sine(0.14943-5.37511*T)
    DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6
    DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6
    DLS =-6.40*S1                                   -1.89*S6
    DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7
    DD  = DL0-DLS
    DGAM  = ((-3332E-9 * Sine(0.59734-5.37261*T)
               -539E-9 * Sine(0.35498-5.37899*T)
                -64E-9 * Sine(0.39943-5.37511*T)))

    L0 = _PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/_ARC
    L  = _PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /_ARC
    LS = _PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/_ARC
    F  = _PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /_ARC
    D  = _PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /_ARC

    I = 1
    while I <= 4:
        if I == 1:
            ARG=L; MAX=4; FAC=1.000002208
        elif I == 2:
            ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T
        elif I == 3:
            ARG=F; MAX=4; FAC=1.000002708+139.978*DGAM
        else:
            ARG=D; MAX=6; FAC=1.0

        ex[0][I] = complex(1, 0)
        ex[1][I] = complex(FAC * math.cos(ARG), FAC * math.sin(ARG))

        J = 2
        while J <= MAX:
            ex[J][I] = ex[J-1][I] * ex[1][I]
            J += 1

        J = 1
        while J <= MAX:
            ex[-J][I] = ex[J][I].conjugate()
            J += 1

        I += 1

$ASTRO_ADDSOL()

    def ADDN(coeffn, p, q, r, s):
        return coeffn * (ex[p][1] * ex[q][2] * ex[r][3] * ex[s][4]).imag

    N = 0
    N += ADDN(-526.069, 0, 0,1,-2)
    N += ADDN(  -3.352, 0, 0,1,-4)
    N += ADDN( +44.297,+1, 0,1,-2)
    N += ADDN(  -6.000,+1, 0,1,-4)
    N += ADDN( +20.599,-1, 0,1, 0)
    N += ADDN( -30.598,-1, 0,1,-2)
    N += ADDN( -24.649,-2, 0,1, 0)
    N += ADDN(  -2.000,-2, 0,1,-2)
    N += ADDN( -22.571, 0,+1,1,-2)
    N += ADDN( +10.985, 0,-1,1,-2)

    DLAM += (
        +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
        +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
        +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
        +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
        +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
        +0.33*Sine(0.3132   +6.3368*T)
    )
    S = F + DS/_ARC
    lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*math.sin(S) - 6.24*math.sin(3*S) + N
    return _moonpos(
        _PI2 * Frac((L0+DLAM/_ARC) / _PI2),
        (math.pi / (180 * 3600)) * lat_seconds,
        (_ARC * _EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
    )

def GeoMoon(time):
    """Calculates the geocentric position of the Moon at a given time.

    Given a time of observation, calculates the Moon's position as a vector.
    The vector gives the location of the Moon's center relative to the Earth's center
    with x-, y-, and z-components measured in astronomical units.

    This algorithm is based on Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
    which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
    It is adapted from Turbo Pascal code from the book
    [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
    by Montenbruck and Pfleger.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's position.

    Returns
    -------
    Vector
        The Moon's position as a vector in J2000 Cartesian equatorial coordinates.

    """
    m = _CalcMoon(time)

    # Convert geocentric ecliptic spherical coordinates to Cartesian coordinates.
    dist_cos_lat = m.distance_au * math.cos(m.geo_eclip_lat)
    gepos = [
        dist_cos_lat * math.cos(m.geo_eclip_lon),
        dist_cos_lat * math.sin(m.geo_eclip_lon),
        m.distance_au * math.sin(m.geo_eclip_lat)
    ]

    # Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
    mpos1 = _ecl2equ_vec(time, gepos)

    # Convert from mean equinox of date to J2000.
    mpos2 = _precession(mpos1, time, _PrecessDir.Into2000)
    return Vector(mpos2[0], mpos2[1], mpos2[2], time)

# END CalcMoon
#----------------------------------------------------------------------------
# BEGIN VSOP

_vsop = [
    # Mercury
    $ASTRO_LIST_VSOP(Mercury),

    # Venus
    $ASTRO_LIST_VSOP(Venus),

    # Earth
    $ASTRO_LIST_VSOP(Earth),

    # Mars
    $ASTRO_LIST_VSOP(Mars),

    # Jupiter
    $ASTRO_LIST_VSOP(Jupiter),

    # Saturn
    $ASTRO_LIST_VSOP(Saturn),

    # Uranus
    $ASTRO_LIST_VSOP(Uranus),

    # Neptune
    $ASTRO_LIST_VSOP(Neptune),
]

def _VsopFormula(formula, t, clamp_angle):
    tpower = 1.0
    coord = 0.0
    for series in formula:
        incr = tpower * sum(A * math.cos(B + C*t) for (A, B, C) in series)
        if clamp_angle:
            # Longitude angles can be hundreds of radians.
            # Improve precision by keeping each increment within [-2*pi, +2*pi].
            incr = math.fmod(incr, _PI2)
        coord += incr
        tpower *= t
    return coord

def _VsopDeriv(formula, t):
    tpower = 1      # t**s
    dpower = 0      # t**(s-1)
    deriv = 0
    s = 0
    for series in formula:
        sin_sum = 0
        cos_sum = 0
        for (ampl, phas, freq) in series:
            angle = phas + (t * freq)
            sin_sum += ampl * freq * math.sin(angle)
            if s > 0:
                cos_sum += ampl * math.cos(angle)
        deriv += (s * dpower * cos_sum) - (tpower * sin_sum)
        dpower = tpower
        tpower *= t
        s += 1
    return deriv

_DAYS_PER_MILLENNIUM = 365250.0
_LON_INDEX = 0
_LAT_INDEX = 1
_RAD_INDEX = 2


def _VsopRotate(eclip):
    # Convert ecliptic cartesian coordinates to equatorial cartesian coordinates.
    x = eclip.x + 0.000000440360*eclip.y - 0.000000190919*eclip.z
    y = -0.000000479966*eclip.x + 0.917482137087*eclip.y - 0.397776982902*eclip.z
    z = 0.397776982902*eclip.y + 0.917482137087*eclip.z
    return _TerseVector(x, y, z)

def _VsopSphereToRect(lon, lat, rad):
    # Convert spherical coordinates to cartesian coordinates.
    r_coslat = rad * math.cos(lat)
    return _TerseVector(
        r_coslat * math.cos(lon),
        r_coslat * math.sin(lon),
        rad * math.sin(lat)
    )

def _CalcVsop(model, time):
    t = time.tt / _DAYS_PER_MILLENNIUM
    lon = _VsopFormula(model[0], t, True)
    lat = _VsopFormula(model[1], t, False)
    rad = _VsopFormula(model[2], t, False)
    eclip = _VsopSphereToRect(lon, lat, rad)
    return _VsopRotate(eclip).ToAstroVector(time)

class _body_state_t:
    def __init__(self, tt, r, v):
        self.tt  = tt
        self.r = r
        self.v = v

def _CalcVsopPosVel(model, tt):
    t = tt / _DAYS_PER_MILLENNIUM

    lon = _VsopFormula(model[0], t, True)
    lat = _VsopFormula(model[1], t, False)
    rad = _VsopFormula(model[2], t, False)

    (dlon_dt, dlat_dt, drad_dt) = [_VsopDeriv(formula, t) for formula in model]

    # Use spherical coords and spherical derivatives to calculate
    # the velocity vector in rectangular coordinates.

    coslon = math.cos(lon)
    sinlon = math.sin(lon)
    coslat = math.cos(lat)
    sinlat = math.sin(lat)

    vx = (
        + (drad_dt * coslat * coslon)
        - (rad * sinlat * coslon * dlat_dt)
        - (rad * coslat * sinlon * dlon_dt)
    )

    vy = (
        + (drad_dt * coslat * sinlon)
        - (rad * sinlat * sinlon * dlat_dt)
        + (rad * coslat * coslon * dlon_dt)
    )

    vz = (
        + (drad_dt * sinlat)
        + (rad * coslat * dlat_dt)
    )

    eclip_pos = _VsopSphereToRect(lon, lat, rad)

    # Convert speed units from [AU/millennium] to [AU/day].
    eclip_vel = _TerseVector(vx, vy, vz) / _DAYS_PER_MILLENNIUM

    # Rotate the vectors from ecliptic to equatorial coordinates.
    equ_pos = _VsopRotate(eclip_pos)
    equ_vel = _VsopRotate(eclip_vel)
    return _body_state_t(tt, equ_pos, equ_vel)


def _AdjustBarycenter(ssb, time, body, pmass):
    shift = pmass / (pmass + _SUN_GM)
    planet = _CalcVsop(_vsop[body.value], time)
    ssb.x += shift * planet.x
    ssb.y += shift * planet.y
    ssb.z += shift * planet.z


def _CalcSolarSystemBarycenter(time):
    ssb = Vector(0.0, 0.0, 0.0, time)
    _AdjustBarycenter(ssb, time, Body.Jupiter, _JUPITER_GM)
    _AdjustBarycenter(ssb, time, Body.Saturn,  _SATURN_GM)
    _AdjustBarycenter(ssb, time, Body.Uranus,  _URANUS_GM)
    _AdjustBarycenter(ssb, time, Body.Neptune, _NEPTUNE_GM)
    return ssb


def _VsopHelioDistance(model, time):
    # The caller only wants to know the distance between the planet and the Sun.
    # So we only need to calculate the radial component of the spherical coordinates.
    # There is no need to translate coordinates.
    return _VsopFormula(model[2], time.tt / _DAYS_PER_MILLENNIUM, False)

def _CalcEarth(time):
    return _CalcVsop(_vsop[Body.Earth.value], time)

# END VSOP
#----------------------------------------------------------------------------
# BEGIN Pluto Integrator

$ASTRO_PLUTO_TABLE()

class _TerseVector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def ToAstroVector(self, time):
        '''Convert _TerseVector object to Vector object.'''
        return Vector(self.x, self.y, self.z, time)

    def quadrature(self):
        '''Return magnitude squared of this vector.'''
        return self.x**2 + self.y**2 + self.z**2

    def mean(self, other):
        '''Return the average of this vector and another vector.'''
        return _TerseVector((self.x + other.x)/2.0, (self.y + other.y)/2.0, (self.z + other.z)/2.0)

    def __add__(self, other):
        return _TerseVector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return _TerseVector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar):
        return _TerseVector(scalar * self.x, scalar * self.y, scalar * self.z)

    def __rmul__(self, scalar):
        return _TerseVector(scalar * self.x, scalar * self.y, scalar * self.z)

    def __truediv__(self, scalar):
        return _TerseVector(self.x / scalar, self.y / scalar, self.z / scalar)


def _BodyStateFromTable(entry):
    [ tt, [rx, ry, rz], [vx, vy, vz] ] = entry
    return _body_state_t(tt, _TerseVector(rx, ry, rz), _TerseVector(vx, vy, vz))


def _AdjustBarycenterPosVel(ssb, tt, body, planet_gm):
    shift = planet_gm / (planet_gm + _SUN_GM)
    planet = _CalcVsopPosVel(_vsop[body.value], tt)
    ssb.r += shift * planet.r
    ssb.v += shift * planet.v
    return planet

def _AccelerationIncrement(small_pos, gm, major_pos):
    delta = major_pos - small_pos
    r2 = delta.quadrature()
    return (gm / (r2 * math.sqrt(r2))) * delta


class _major_bodies_t:
    def __init__(self, tt):
        # Accumulate the Solar System Barycenter position.
        ssb = _body_state_t(tt, _TerseVector(0,0,0), _TerseVector(0,0,0))
        # Calculate the position and velocity vectors of the 4 major planets.
        self.Jupiter = _AdjustBarycenterPosVel(ssb, tt, Body.Jupiter, _JUPITER_GM)
        self.Saturn  = _AdjustBarycenterPosVel(ssb, tt, Body.Saturn,  _SATURN_GM)
        self.Uranus  = _AdjustBarycenterPosVel(ssb, tt, Body.Uranus,  _URANUS_GM)
        self.Neptune = _AdjustBarycenterPosVel(ssb, tt, Body.Neptune, _NEPTUNE_GM)
        # Convert the planets' vectors from heliocentric to barycentric.
        self.Jupiter.r -= ssb.r
        self.Jupiter.v -= ssb.v
        self.Saturn.r  -= ssb.r
        self.Saturn.v  -= ssb.v
        self.Uranus.r  -= ssb.r
        self.Uranus.v  -= ssb.v
        self.Neptune.r -= ssb.r
        self.Neptune.v -= ssb.v
        # Convert heliocentric SSB to barycentric Sun.
        self.Sun = _body_state_t(tt, -1*ssb.r, -1*ssb.v)

    def Acceleration(self, pos):
        '''Use barycentric coordinates of the Sun and major planets to calculate
        the gravitational acceleration vector experienced at location 'pos'.'''
        acc  = _AccelerationIncrement(pos, _SUN_GM,     self.Sun.r)
        acc += _AccelerationIncrement(pos, _JUPITER_GM, self.Jupiter.r)
        acc += _AccelerationIncrement(pos, _SATURN_GM,  self.Saturn.r)
        acc += _AccelerationIncrement(pos, _URANUS_GM,  self.Uranus.r)
        acc += _AccelerationIncrement(pos, _NEPTUNE_GM, self.Neptune.r)
        return acc


class _body_grav_calc_t:
    def __init__(self, tt, r, v, a):
        self.tt = tt    # J2000 terrestrial time [days]
        self.r = r      # position [au]
        self.v = v      # velocity [au/day]
        self.a = a      # acceleration [au/day^2]


class _grav_sim_t:
    def __init__(self, bary, grav):
        self.bary = bary
        self.grav = grav


def _UpdatePosition(dt, r, v, a):
    return _TerseVector(
        r.x + dt*(v.x + dt*a.x/2.0),
        r.y + dt*(v.y + dt*a.y/2.0),
        r.z + dt*(v.z + dt*a.z/2.0)
    )

def _UpdateVelocity(dt, v, a):
    return _TerseVector(
        v.x + dt*a.x,
        v.y + dt*a.y,
        v.z + dt*a.z
    )


def _GravSim(tt2, calc1):
    dt = tt2 - calc1.tt

    # Calculate where the major bodies (Sun, Jupiter...Neptune) will be at tt2.
    bary2 = _major_bodies_t(tt2)

    # Estimate position of small body as if current acceleration applies across the whole time interval.
    approx_pos = _UpdatePosition(dt, calc1.r, calc1.v, calc1.a)

    # Calculate the average acceleration of the endpoints.
    # This becomes our estimate of the mean effective acceleration over the whole interval.
    mean_acc = bary2.Acceleration(approx_pos).mean(calc1.a)

    # Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration.
    pos = _UpdatePosition(dt, calc1.r, calc1.v, mean_acc)
    vel = calc1.v + dt*mean_acc
    acc = bary2.Acceleration(pos)
    grav = _body_grav_calc_t(tt2, pos, vel, acc)
    return _grav_sim_t(bary2, grav)


_pluto_cache = [None] * (_PLUTO_NUM_STATES - 1)


def _ClampIndex(frac, nsteps):
    index = math.floor(frac)
    if index < 0:
        return 0
    if index >= nsteps:
        return nsteps-1
    return index


def _GravFromState(entry):
    state = _BodyStateFromTable(entry)
    bary = _major_bodies_t(state.tt)
    r = state.r + bary.Sun.r
    v = state.v + bary.Sun.v
    a = bary.Acceleration(r)
    grav = _body_grav_calc_t(state.tt, r, v, a)
    return _grav_sim_t(bary, grav)


def _GetSegment(cache, tt):
    if (tt < _PlutoStateTable[0][0]) or (tt > _PlutoStateTable[_PLUTO_NUM_STATES-1][0]):
        # Don't bother calculating a segment. Let the caller crawl backward/forward to this time.
        return None

    seg_index = _ClampIndex((tt - _PlutoStateTable[0][0]) / _PLUTO_TIME_STEP, _PLUTO_NUM_STATES-1)
    if cache[seg_index] is None:
        seg = cache[seg_index] = [None] * _PLUTO_NSTEPS

        # Each endpoint is exact.
        seg[0] = _GravFromState(_PlutoStateTable[seg_index]).grav
        seg[_PLUTO_NSTEPS-1] = _GravFromState(_PlutoStateTable[seg_index + 1]).grav

        # Simulate forwards from the lower time bound.
        step_tt = seg[0].tt
        i = 1
        while i < _PLUTO_NSTEPS-1:
            step_tt += _PLUTO_DT
            seg[i] = _GravSim(step_tt, seg[i-1]).grav
            i += 1

        # Simulate backwards from the upper time bound.
        step_tt = seg[_PLUTO_NSTEPS-1].tt
        reverse = [None] * _PLUTO_NSTEPS
        reverse[_PLUTO_NSTEPS-1] = seg[_PLUTO_NSTEPS-1]
        i = _PLUTO_NSTEPS - 2
        while i > 0:
            step_tt -= _PLUTO_DT
            reverse[i] = _GravSim(step_tt, reverse[i+1]).grav
            i -= 1

        # Fade-mix the two series so that there are no discontinuities.
        i = _PLUTO_NSTEPS - 2
        while i > 0:
            ramp = i / (_PLUTO_NSTEPS-1)
            seg[i].r = seg[i].r*(1 - ramp) + reverse[i].r*ramp
            seg[i].v = seg[i].v*(1 - ramp) + reverse[i].v*ramp
            seg[i].a = seg[i].a*(1 - ramp) + reverse[i].a*ramp
            i -= 1

    return cache[seg_index]


def _CalcPlutoOneWay(entry, target_tt, dt):
    sim = _GravFromState(entry)
    n = math.ceil((target_tt - sim.grav.tt) / dt)
    for i in range(n):
        sim = _GravSim(target_tt if (i+1 == n) else (sim.grav.tt + dt), sim.grav)
    return sim


def _CalcPluto(time, helio):
    bary = None
    seg = _GetSegment(_pluto_cache, time.tt)
    if seg is None:
        # The target time is outside the year range 0000..4000.
        # Calculate it by crawling backward from 0000 or forward from 4000.
        # FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed.
        if time.tt < _PlutoStateTable[0][0]:
            sim = _CalcPlutoOneWay(_PlutoStateTable[0], time.tt, -_PLUTO_DT)
        else:
            sim = _CalcPlutoOneWay(_PlutoStateTable[_PLUTO_NUM_STATES-1], time.tt, +_PLUTO_DT)
        r = sim.grav.r
        v = sim.grav.v
        bary = sim.bary
    else:
        left = _ClampIndex((time.tt - seg[0].tt) / _PLUTO_DT, _PLUTO_NSTEPS-1)
        s1 = seg[left]
        s2 = seg[left+1]

        # Find mean acceleration vector over the interval.
        acc = s1.a.mean(s2.a)

        # Use Newtonian mechanics to extrapolate away from t1 in the positive time direction.
        ra = _UpdatePosition(time.tt - s1.tt, s1.r, s1.v, acc)
        va = _UpdateVelocity(time.tt - s1.tt, s1.v, acc)

        # Use Newtonian mechanics to extrapolate away from t2 in the negative time direction.
        rb = _UpdatePosition(time.tt - s2.tt, s2.r, s2.v, acc)
        vb = _UpdateVelocity(time.tt - s2.tt, s2.v, acc)

        # Use fade in/out idea to blend the two position estimates.
        ramp = (time.tt - s1.tt)/_PLUTO_DT
        r = ra*(1 - ramp) + rb*ramp
        v = va*(1 - ramp) + vb*ramp

    if helio:
        # Convert barycentric vectors to heliocentric vectors.
        if bary is None:
            bary = _major_bodies_t(time.tt)
        r -= bary.Sun.r
        v -= bary.Sun.v

    return StateVector(r.x, r.y, r.z, v.x, v.y, v.z, time)


# END Pluto Integrator
#----------------------------------------------------------------------------
# BEGIN Jupiter Moons

$ASTRO_JUPITER_MOONS()

class JupiterMoonsInfo:
    """Holds the positions and velocities of Jupiter's major 4 moons.

    The #JupiterMoons function returns an object of this type
    to report position and velocity vectors for Jupiter's largest 4 moons
    Io, Europa, Ganymede, and Callisto. Each position vector is relative
    to the center of Jupiter. Both position and velocity are oriented in
    the EQJ system (that is, using Earth's equator at the J2000 epoch).
    The positions are expressed in astronomical units (AU),
    and the velocities in AU/day.

    Attributes
    ----------
    moon : StateVector[4]
        An array of state vectors, one for each of the four major moons
        of Jupiter, in the following order: 0=Io, 1=Europa, 2=Ganymede, 3=Callisto.
    """
    def __init__(self, moon):
        self.moon = moon


def _JupiterMoon_elem2pv(time, mu, A, AL, K, H, Q, P):
    # Translation of FORTRAN subroutine ELEM2PV from:
    # https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/
    AN = math.sqrt(mu / (A*A*A))
    EE = AL + K*math.sin(AL) - H*math.cos(AL)
    DE = 1
    while abs(DE) >= 1.0e-12:
        CE = math.cos(EE)
        SE = math.sin(EE)
        DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE)
        EE += DE
    CE = math.cos(EE)
    SE = math.sin(EE)
    DLE = H*CE - K*SE
    RSAM1 = -K*CE - H*SE
    ASR = 1.0/(1.0 + RSAM1)
    PHI = math.sqrt(1.0 - K*K - H*H)
    PSI = 1.0/(1.0 + PHI)
    X1 = A*(CE - K - PSI*H*DLE)
    Y1 = A*(SE - H + PSI*K*DLE)
    VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1)
    VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1)
    F2 = 2.0*math.sqrt(1.0 - Q*Q - P*P)
    P2 = 1.0 - 2.0*P*P
    Q2 = 1.0 - 2.0*Q*Q
    PQ = 2.0*P*Q
    return StateVector(
        X1*P2 + Y1*PQ,
        X1*PQ + Y1*Q2,
        (Q*Y1 - X1*P)*F2,
        VX1*P2 + VY1*PQ,
        VX1*PQ + VY1*Q2,
        (Q*VY1 - VX1*P)*F2,
        time
    )


def _CalcJupiterMoon(time, mu, al0, al1, a, l, z, zeta):
    # This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
    # https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    t = time.tt + 18262.5    # number of days since 1950-01-01T00:00:00Z

    # Calculate 6 orbital elements at the given time t.
    elem0 = 0.0
    for (amplitude, phase, frequency) in a:
        elem0 += amplitude * math.cos(phase + (t * frequency))

    elem1 = al0 + (t * al1)
    for (amplitude, phase, frequency) in l:
        elem1 += amplitude * math.sin(phase + (t * frequency))

    elem1 = math.fmod(elem1, _PI2)
    if elem1 < 0:
        elem1 += _PI2

    elem2 = 0.0
    elem3 = 0.0
    for (amplitude, phase, frequency) in z:
        arg = phase + (t * frequency)
        elem2 += amplitude * math.cos(arg)
        elem3 += amplitude * math.sin(arg)

    elem4 = 0.0
    elem5 = 0.0
    for (amplitude, phase, frequency) in zeta:
        arg = phase + (t * frequency)
        elem4 += amplitude * math.cos(arg)
        elem5 += amplitude * math.sin(arg)

    # Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
    state = _JupiterMoon_elem2pv(time, mu, elem0, elem1, elem2, elem3, elem4, elem5)

    # Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
    return RotateState(_Rotation_JUP_EQJ, state)


def JupiterMoons(time):
    """Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.

    Calculates position and velocity vectors for Jupiter's moons
    Io, Europa, Ganymede, and Callisto, at the given date and time.
    The vectors are jovicentric (relative to the center of Jupiter).
    Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
    The position components are expressed in astronomical units (AU), and the
    velocity components are in AU/day.

    To convert to heliocentric vectors, call #HelioVector
    with `Body.Jupiter` to get Jupiter's heliocentric position, then
    add the jovicentric vectors. Likewise, you can call #GeoVector
    to convert to geocentric vectors.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate Jupiter's moons.

    Returns
    -------
    JupiterMoonsInfo
        The positions and velocities of Jupiter's 4 largest moons.
    """
    infolist = []
    for (mu, al0, al1, a, l, z, zeta) in _JupiterMoonModel:
        infolist.append(_CalcJupiterMoon(time, mu, al0, al1, a, l, z, zeta))
    return JupiterMoonsInfo(infolist)

# END Jupiter Moons
#----------------------------------------------------------------------------
# BEGIN Search

def _QuadInterp(tm, dt, fa, fm, fb):
    Q = (fb + fa)/2 - fm
    R = (fb - fa)/2
    S = fm

    if Q == 0:
        # This is a line, not a parabola.
        if R == 0:
            # This is a HORIZONTAL line... can't make progress!
            return None
        x = -S / R
        if not (-1 <= x <= +1):
            return None  # out of bounds
    else:
        # It really is a parabola. Find roots x1, x2.
        u = R*R - 4*Q*S
        if u <= 0:
            return None
        ru = math.sqrt(u)
        x1 = (-R + ru) / (2 * Q)
        x2 = (-R - ru) / (2 * Q)

        if -1 <= x1 <= +1:
            if -1 <= x2 <= +1:
                # Two solutions... so parabola intersects twice.
                return None
            x = x1
        elif -1 <= x2 <= +1:
            x = x2
        else:
            return None

    t = tm + x*dt
    df_dt = (2*Q*x + R) / dt
    return (t, df_dt)

def Search(func, context, t1, t2, dt_tolerance_seconds):
    """Searches for a time at which a function's value increases through zero.

    Certain astronomy calculations involve finding a time when an event occurs.
    Often such events can be defined as the root of a function:
    the time at which the function's value becomes zero.

    `Search` finds the *ascending root* of a function: the time at which
    the function's value becomes zero while having a positive slope. That is, as time increases,
    the function transitions from a negative value, through zero at a specific moment,
    to a positive value later. The goal of the search is to find that specific moment.

    The search function is specified by two parameters: `func` and `context`.
    The `func` parameter is a function itself that accepts a time
    and a context containing any other arguments needed to evaluate the function.
    The `context` parameter supplies that context for the given search.
    As an example, a caller may wish to find the moment a celestial body reaches a certain
    ecliptic longitude. In that case, the caller might create a type (class, tuple, whatever)
    that contains a #Body member to specify the body and a numeric value to hold the target longitude.
    A different function might use a completely different context type.

    Every time it is called, `func` returns a `float` value or it raises an exception.
    If `func` raises an exception, the search immediately fails and the exception is
    propagated back to the caller.
    Otherwise, the search proceeds until it either finds the ascending root or fails for some reason.

    The search calls `func` repeatedly to rapidly narrow in on any ascending
    root within the time window specified by `t1` and `t2`. The search never
    reports a solution outside this time window.

    `Search` uses a combination of bisection and quadratic interpolation
    to minimize the number of function calls. However, it is critical that the
    supplied time window be small enough that there cannot be more than one root
    (ascedning or descending) within it; otherwise the search can fail.
    Beyond that, it helps to make the time window as small as possible, ideally
    such that the function itself resembles a smooth parabolic curve within that window.

    If an ascending root is not found, or more than one root
    (ascending and/or descending) exists within the window `t1`..`t2`,
    `Search` will return `None` to indicate a normal search failure.

    If the search does not converge within 20 iterations, it will raise
    an #Error exception.

    Parameters
    ----------
    func : function(context, Time)
        A function that takes an arbitrary context parameter and a #Time parameter.
        Returns a float value.  See remarks above for more details.

    context : object
        An arbitrary data structure needed to be passed to the function `func`
        every time it is called.

    t1 : float
        The lower time bound of the search window.
        See remarks above for more details.

    t2 : float
        The upper time bound of the search window.
        See remarks above for more details.

    dt_tolerance_seconds : float
        Specifies an amount of time in seconds within which a bounded ascending root
        is considered accurate enough to stop. A typical value is 1 second.

    Returns
    -------
    Time or `None`
        If the search is successful, returns a #Time object that is within
        `dt_tolerance_seconds` of an ascending root.
        In this case, the returned time value will always be within the
        inclusive range [`t1`, `t2`].
        If there is no ascending root, or there is more than one ascending root,
        the function returns `None`.

    """
    dt_days = abs(dt_tolerance_seconds / _SECONDS_PER_DAY)
    f1 = func(context, t1)
    f2 = func(context, t2)
    iter_count = 0
    iter_limit = 20
    calc_fmid = True
    while True:
        iter_count += 1
        if iter_count > iter_limit:
            raise Error('Excessive iteration in Search')

        dt = (t2.tt - t1.tt) / 2.0
        tmid = t1.AddDays(dt)
        if abs(dt) < dt_days:
            # We are close enough to the event to stop the search.
            return tmid

        if calc_fmid:
            fmid = func(context, tmid)
        else:
            # We already have the correct value of fmid from the previous loop.
            calc_fmid = True

        # Quadratic interpolation:
        # Try to find a parabola that passes through the 3 points we have sampled:
        # (t1,f1), (tmid,fmid), (t2,f2).
        q = _QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2)
        if q:
            (q_ut, q_df_dt) = q
            tq = Time(q_ut)
            fq = func(context, tq)
            if q_df_dt != 0.0:
                dt_guess = abs(fq / q_df_dt)
                if dt_guess < dt_days:
                    # The estimated time error is small enough that we can quit now.
                    return tq

                # Try guessing a tighter boundary with the interpolated root at the center.
                dt_guess *= 1.2
                if dt_guess < dt / 10.0:
                    tleft = tq.AddDays(-dt_guess)
                    tright = tq.AddDays(+dt_guess)
                    if (tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0.0:
                        if (tright.ut - t1.ut)*(tright.ut - t2.ut) < 0.0:
                            fleft = func(context, tleft)
                            fright = func(context, tright)
                            if fleft < 0.0 and fright >= 0.0:
                                f1 = fleft
                                f2 = fright
                                t1 = tleft
                                t2 = tright
                                fmid = fq
                                calc_fmid = False
                                continue

        # Quadratic interpolation attempt did not work out.
        # Just divide the region in two parts and pick whichever one appears to contain a root.
        if f1 < 0.0 and fmid >= 0.0:
            t2 = tmid
            f2 = fmid
            continue

        if fmid < 0.0 and f2 >= 0.0:
            t1 = tmid
            f1 = fmid
            continue

        # Either there is no ascending zero-crossing in this range
        # or the search window is too wide (more than one zero-crossing).
        return None

# END Search
#----------------------------------------------------------------------------


def HelioVector(body, time):
    """Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Sun as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    The position is not corrected for light travel time or aberration.
    This is different from the behavior of #GeoVector.

    If given an invalid value for `body`, this function raises an exception.

    Parameters
    ----------
    body : Body
        The celestial body whose heliocentric position is to be calculated:
        The Sun, Moon, EMB, SSB, or any of the planets.
    time : Time
        The time at which to calculate the heliocentric position.

    Returns
    -------
    Vector
        A heliocentric position vector of the center of the given body
        at the given time.
    """
    if body == Body.Pluto:
        planet = _CalcPluto(time, True)
        return Vector(planet.x, planet.y, planet.z, time)

    if 0 <= body.value < len(_vsop):
        return _CalcVsop(_vsop[body.value], time)

    if body == Body.Sun:
        return Vector(0.0, 0.0, 0.0, time)

    if body == Body.Moon:
        e = _CalcEarth(time)
        m = GeoMoon(time)
        return Vector(e.x+m.x, e.y+m.y, e.z+m.z, time)

    if body == Body.EMB:
        e = _CalcEarth(time)
        m = GeoMoon(time)
        d = 1.0 + _EARTH_MOON_MASS_RATIO
        return Vector(e.x+(m.x/d), e.y+(m.y/d), e.z+(m.z/d), time)

    if body == Body.SSB:
        return _CalcSolarSystemBarycenter(time)

    raise InvalidBodyError()


def HelioDistance(body, time):
    """Calculates the distance between a body and the Sun at a given time.

    Given a date and time, this function calculates the distance between
    the center of `body` and the center of the Sun.
    For the planets Mercury through Neptune, this function is significantly
    more efficient than calling #HelioVector followed by taking the length
    of the resulting vector.

    Parameters
    ----------
    body : Body
        A body for which to calculate a heliocentric distance:
        the Sun, Moon, or any of the planets.
    time : Time
        The date and time for which to calculate the heliocentric distance.

    Returns
    -------
    float
        The heliocentric distance in AU.
    """
    if body == Body.Sun:
        return 0.0

    if 0 <= body.value < len(_vsop):
        return _VsopHelioDistance(_vsop[body.value], time)

    return HelioVector(body, time).Length()


def GeoVector(body, time, aberration):
    """Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.

    This function calculates the position of the given celestial body as a vector,
    using the center of the Earth as the origin.  The result is expressed as a Cartesian
    vector in the J2000 equatorial system: the coordinates are based on the mean equator
    of the Earth at noon UTC on 1 January 2000.

    If given an invalid value for `body`, this function will raise an exception.

    Unlike #HelioVector, this function always corrects for light travel time.
    This means the position of the body is "back-dated" by the amount of time it takes
    light to travel from that body to an observer on the Earth.

    Also, the position can optionally be corrected for
    [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
    causing the apparent direction of the body to be shifted due to transverse
    movement of the Earth with respect to the rays of light coming from that body.

    Parameters
    ----------
    body : Body
        A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
    time : Time
        The date and time for which to calculate the position.
    aberration : bool
        A boolean value indicating whether to correct for aberration.

    Returns
    -------
    Vector
        A geocentric position vector of the center of the given body.
    """
    if body == Body.Moon:
        return GeoMoon(time)

    if body == Body.Earth:
        return Vector(0.0, 0.0, 0.0, time)

    if not aberration:
        # No aberration, so calculate Earth's position once, at the time of observation.
        earth = _CalcEarth(time)

    # Correct for light-travel time, to get position of body as seen from Earth's center.
    ltime = time
    for _ in range(10):
        h = HelioVector(body, ltime)
        if aberration:
            # Include aberration, so make a good first-order approximation
            # by backdating the Earth's position also.
            # This is confusing, but it works for objects within the Solar System
            # because the distance the Earth moves in that small amount of light
            # travel time (a few minutes to a few hours) is well approximated
            # by a line segment that substends the angle seen from the remote
            # body viewing Earth. That angle is pretty close to the aberration
            # angle of the moving Earth viewing the remote body.
            # In other words, both of the following approximate the aberration angle:
            #    (transverse distance Earth moves) / (distance to body)
            #    (transverse speed of Earth) / (speed of light).
            earth = _CalcEarth(ltime)

        geo = Vector(h.x-earth.x, h.y-earth.y, h.z-earth.z, time)
        ltime2 = time.AddDays(-geo.Length() / C_AUDAY)
        dt = abs(ltime2.tt - ltime.tt)
        if dt < 1.0e-9:
            return geo

        ltime = ltime2

    raise Error('Light-travel time solver did not converge: dt={}'.format(dt))


def _ExportState(terse, time):
    return StateVector(
        terse.r.x, terse.r.y, terse.r.z,
        terse.v.x, terse.v.y, terse.v.z,
        time
    )


def BaryState(body, time):
    """Calculates barycentric position and velocity vectors for the given body.

    Given a body and a time, calculates the barycentric position and velocity
    vectors for the center of that body at that time.
    The vectors are expressed in equatorial J2000 coordinates (EQJ).

    Parameters
    ----------
    body : Body
        The celestial body whose barycentric state vector is to be calculated.
        Supported values are `Body.Sun`, `Body.SSB`, and all planets:
        `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
        `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
    time : Time
        The date and time for which to calculate position and velocity.

    Returns
    -------
    StateVector
        An object that contains barycentric position and velocity vectors.
    """
    # Trivial case: the solar sytem barycenter itself.
    if body == Body.SSB:
        return StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)

    if body == Body.Pluto:
        return _CalcPluto(time, False)

    # Find the barycentric positions and velocities for the 5 major bodies.
    bary = _major_bodies_t(time.tt)

    # If the caller is asking for one of the major bodies,
    # we can immediately return the answer.
    if body == Body.Sun:
        return _ExportState(bary.Sun, time)

    if body == Body.Jupiter:
        return _ExportState(bary.Jupiter, time)

    if body == Body.Saturn:
        return _ExportState(bary.Saturn, time)

    if body == Body.Uranus:
        return _ExportState(bary.Uranus, time)

    if body == Body.Neptune:
        return _ExportState(bary.Neptune, time)

    if 0 <= body.value < len(_vsop):
        # Handle the remaining VSOP bodies: Mercury, Venus, Earth, Mars.
        # Calculate the heliocentric state of the given body
        # and add the Sun's heliocentric state to obtain the
        # body's barycentric state.
        # BarySun + HelioBody = BaryBody
        planet = _CalcVsopPosVel(_vsop[body.value], time.tt)
        return StateVector(
            bary.Sun.r.x + planet.r.x,
            bary.Sun.r.y + planet.r.y,
            bary.Sun.r.z + planet.r.z,
            bary.Sun.v.x + planet.v.x,
            bary.Sun.v.y + planet.v.y,
            bary.Sun.v.z + planet.v.z,
            time
        )

    # FIXFIXFIX: later, we can add support for Moon, EMB, etc.
    raise InvalidBodyError()


def Equator(body, time, observer, ofdate, aberration):
    """Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.

    Calculates topocentric equatorial coordinates in one of two different systems:
    J2000 or true-equator-of-date, depending on the value of the `ofdate` parameter.
    Equatorial coordinates include right ascension, declination, and distance in astronomical units.

    This function corrects for light travel time: it adjusts the apparent location
    of the observed body based on how long it takes for light to travel from the body to the Earth.

    This function corrects for *topocentric parallax*, meaning that it adjusts for the
    angular shift depending on where the observer is located on the Earth. This is most
    significant for the Moon, because it is so close to the Earth. However, parallax corection
    has a small effect on the apparent positions of other bodies.

    Correction for aberration is optional, using the `aberration` parameter.

    Parameters
    ----------
    body : Body
        The celestial body to be observed. Not allowed to be `Body.Earth`.
    time : Time
        The date and time at which the observation takes place.
    observer : Observer
        A location on or near the surface of the Earth.
    ofdate : bool
        Selects the date of the Earth's equator in which to express the equatorial coordinates.
        If `True`, returns coordinates using the equator and equinox of date.
        If `False`, returns coordinates converted to the J2000 system.
    aberration : bool
        If `True`, corrects for aberration of light based on the motion of the Earth
        with respect to the heliocentric origin.
        If `False`, does not correct for aberration.

    Returns
    -------
    Equatorial
        Equatorial coordinates in the specified frame of reference.
    """
    gc_observer = _geo_pos(time, observer)
    gc = GeoVector(body, time, aberration)
    j2000 = [
        gc.x - gc_observer[0],
        gc.y - gc_observer[1],
        gc.z - gc_observer[2]
    ]
    if not ofdate:
        return _vector2radec(j2000, time)
    temp = _precession(j2000, time, _PrecessDir.From2000)
    datevect = _nutation(temp, time, _PrecessDir.From2000)
    return _vector2radec(datevect, time)


def ObserverVector(time, observer, ofdate):
    """Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.

    This function calculates a vector from the center of the Earth to
    a point on or near the surface of the Earth, expressed in equatorial
    coordinates. It takes into account the rotation of the Earth at the given
    time, along with the given latitude, longitude, and elevation of the observer.

    The caller may pass `ofdate` as `True` to return coordinates relative to the Earth's
    equator at the specified time, or `False` to use the J2000 equator.

    The returned vector has components expressed in astronomical units (AU).
    To convert to kilometers, multiply the `x`, `y`, and `z` values by
    the constant value #KM_PER_AU.

    The inverse of this function is also available: #VectorObserver.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the observer's position vector.
    observer : Observer
        The geographic location of a point on or near the surface of the Earth.
    ofdate : bool
        Selects the date of the Earth's equator in which to express the equatorial coordinates.
        The caller may pass `False` to use the orientation of the Earth's equator
        at noon UTC on January 1, 2000, in which case this function corrects for precession
        and nutation of the Earth as it was at the moment specified by the `time` parameter.
        Or the caller may pass `true` to use the Earth's equator at `time`
        as the orientation.

    Returns
    -------
    Vector
        An equatorial vector from the center of the Earth to the specified location
        on (or near) the Earth's surface.
    """
    gast = _sidereal_time(time)
    ovec = _terra(observer, gast)
    if not ofdate:
        ovec = _nutation(ovec, time, _PrecessDir.Into2000)
        ovec = _precession(ovec, time, _PrecessDir.Into2000)
    return Vector(ovec[0], ovec[1], ovec[2], time)

def VectorObserver(vector, ofdate):
    """Calculates the geographic location corresponding to an equatorial vector.

    This is the inverse function of #ObserverVector.
    Given a geocentric equatorial vector, it returns the geographic
    latitude, longitude, and elevation for that vector.

    Parameters
    ----------
    vector : Vector
        The geocentric equatorial position vector for which to find geographic coordinates.
        The components are expressed in astronomical units (AU).
        The time `vector.t` determines the Earth's rotation.
    ofdate : bool
        Selects the date of the Earth's equator in which `vector` is expressed.
        The caller may pass `False` to use the orientation of the Earth's equator
        at noon UTC on January 1, 2000, in which case this function corrects for precession
        and nutation of the Earth as it was at the moment specified by the the time `vector.t`.
        Or the caller may pass `True` to use the Earth's equator at `vector.t` as the orientation.

    Returns
    -------
    Observer
        The geographic latitude, longitude, and elevation above sea level
        that corresponds to the given equatorial vector.
    """
    gast = _sidereal_time(vector.t)
    ovec = [vector.x, vector.y, vector.z]
    if not ofdate:
        ovec = _precession(ovec, vector.t, _PrecessDir.From2000)
        ovec = _nutation(ovec, vector.t, _PrecessDir.From2000)
    return _inverse_terra(ovec, gast)

def ObserverGravity(latitude, height):
    """Calculates the gravitational acceleration experienced by an observer on the Earth.

    This function implements the WGS 84 Ellipsoidal Gravity Formula.
    The result is a combination of inward gravitational acceleration
    with outward centrifugal acceleration, as experienced by an observer
    in the Earth's rotating frame of reference.
    The resulting value increases toward the Earth's poles and decreases
    toward the equator, consistent with changes of the weight measured
    by a spring scale of a fixed mass moved to different latitudes and heights
    on the Earth.

    Parameters
    ----------
    latitude : float
        The latitude of the observer in degrees north or south of the equator.
        By formula symmetry, positive latitudes give the same answer as negative
        latitudes, so the sign does not matter.
    height : float
        The height above the sea level geoid in meters.
        No range checking is done; however, accuracy is only valid in the
        range 0 to 100000 meters.

    Returns
    -------
    float
        The effective gravitational acceleration expressed in meters per second squared [m/s^2].
    """
    s2 = math.sin(math.radians(latitude)) ** 2
    g0 = 9.7803253359 * (1.0 + 0.00193185265241*s2) / math.sqrt(1.0 - 0.00669437999013*s2)
    return g0 * (1.0 - (3.15704e-07 - 2.10269e-09*s2)*height + 7.37452e-14*height*height)

@enum.unique
class Refraction(enum.Enum):
    """Selects if/how to correct for atmospheric refraction.

    Some functions allow enabling or disabling atmospheric refraction
    for the calculated apparent position of a celestial body
    as seen by an observer on the surface of the Earth.

    Values
    ------
    Airless:      No atmospheric refraction correction.
    Normal:       Recommended correction for standard atmospheric refraction.
    JplHorizons:  Used only for compatibility testing with JPL Horizons online tool.
    """
    Airless = 0
    Normal = 1
    JplHorizons = 2

class HorizontalCoordinates:
    """Coordinates of a celestial body as seen by a topocentric observer.

    Contains horizontal and equatorial coordinates as seen by an observer
    on or near the surface of the Earth (a topocentric observer).
    All coordinates are optionally corrected for atmospheric refraction.

    Attributes
    ----------
    azimuth : float
        The compass direction laterally around the observer's horizon,
        measured in degrees.
        North is 0 degrees, east is 90 degrees, south is 180 degrees, etc.
    altitude : float
        The angle in degrees above (positive) or below (negative) the observer's horizon.
    ra : float
        The right ascension in sidereal hours.
    dec : float
        The declination in degrees.
    """
    def __init__(self, azimuth, altitude, ra, dec):
        self.azimuth = azimuth
        self.altitude = altitude
        self.ra = ra
        self.dec = dec

def Horizon(time, observer, ra, dec, refraction):
    """Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.

    Given a date and time, the geographic location of an observer on the Earth, and
    equatorial coordinates (right ascension and declination) of a celestial body,
    this function returns horizontal coordinates (azimuth and altitude angles) for the body
    relative to the horizon at the geographic location.

    The right ascension `ra` and declination `dec` passed in must be *equator of date*
    coordinates, based on the Earth's true equator at the date and time of the observation.
    Otherwise the resulting horizontal coordinates will be inaccurate.
    Equator of date coordinates can be obtained by calling #Equator, passing in
    `True` as its `ofdate` parameter. It is also recommended to enable
    aberration correction by passing in `True` for the `aberration` parameter.

    This function optionally corrects for atmospheric refraction.
    For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
    correct for optical lensing of the Earth's atmosphere that causes objects
    to appear somewhat higher above the horizon than they actually are.
    However, callers may choose to avoid this correction by passing in `Refraction.Airless`.
    If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
    in the #HorizontalCoordinates object returned by this function will all be corrected for refraction.
    If refraction is disabled, none of these four coordinates will be corrected; in that case,
    the right ascension and declination in the returned object will be numerically identical
    to the respective `ra` and `dec` values passed in.

    Parameters
    ----------
    time : Time
        The date and time for which to find horizontal coordinates.
    observer : Observer
        The location of the observer for which to find horizontal coordinates.
    ra : float
        Right ascension in sidereal hours of the celestial object,
        referred to the mean equinox of date for the J2000 epoch.
    dec : float
        Declination in degrees of the celestial object,
        referred to the mean equator of date for the J2000 epoch.
        Positive values are north of the celestial equator
        and negative values are south of it.
    refraction : Refraction
        The option for selecting whether to correct for atmospheric lensing.
        If `Refraction.Normal`, a well-behaved refraction model is used.
        If `Refraction.None`, no refraction correct is performed.
        `Refraction.JplHorizons` is used only for compatibility testing
        with the JPL Horizons online tool.

    Returns
    -------
    HorizontalCoordinates
        The horizontal coordinates (altitude and azimuth), along with
        equatorial coordinates (right ascension and declination), all
        optionally corrected for atmospheric refraction. See remarks above
        for more details.
    """
    if not (Refraction.Airless.value <= refraction.value <= Refraction.JplHorizons.value):
        raise Error('Invalid refraction type')

    latrad = math.radians(observer.latitude)
    lonrad = math.radians(observer.longitude)
    decrad = math.radians(dec)
    rarad = ra * _HOUR2RAD

    sinlat = math.sin(latrad)
    coslat = math.cos(latrad)
    sinlon = math.sin(lonrad)
    coslon = math.cos(lonrad)
    sindc = math.sin(decrad)
    cosdc = math.cos(decrad)
    sinra = math.sin(rarad)
    cosra = math.cos(rarad)

    # Calculate three mutually perpendicular unit vectors
    # in equatorial coordinates: uze, une, uwe.
    #
    # uze = The direction of the observer's local zenith (straight up).
    # une = The direction toward due north on the observer's horizon.
    # uwe = The direction toward due west on the observer's horizon.
    #
    # HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.
    #
    # The components of these 3 vectors are as follows:
    # [0] = x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
    # [1] = y = direction from center of Earth toward 90 degrees west longitude on equator.
    # [2] = z = direction from center of Earth toward the north pole.

    uze = [coslat*coslon, coslat*sinlon, sinlat]
    une = [-sinlat*coslon, -sinlat*sinlon, coslat]
    uwe = [sinlon, -coslon, 0.0]

    # Correct the vectors uze, une, uwe for the Earth's rotation by calculating
    # sideral time. Call spin() for each uncorrected vector to rotate about
    # the Earth's axis to yield corrected unit vectors uz, un, uw.
    # Multiply sidereal hours by -15 to convert to degrees and flip eastward
    # rotation of the Earth to westward apparent movement of objects with time.

    angle = -15.0 * _sidereal_time(time)
    uz = _spin(angle, uze)
    un = _spin(angle, une)
    uw = _spin(angle, uwe)

    # Convert angular equatorial coordinates (RA, DEC) to
    # cartesian equatorial coordinates in 'p', using the
    # same orientation system as uze, une, uwe.

    p = [cosdc*cosra, cosdc*sinra, sindc]

    # Use dot products of p with the zenith, north, and west
    # vectors to obtain the cartesian coordinates of the body in
    # the observer's horizontal orientation system.
    #
    # pz = zenith component [-1, +1]
    # pn = north  component [-1, +1]
    # pw = west   component [-1, +1]

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2]
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2]
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2]

    # proj is the "shadow" of the body vector along the observer's flat ground.
    proj = math.sqrt(pn*pn + pw*pw)

    # Calculate az = azimuth (compass direction clockwise from East.)
    if proj > 0.0:
        # If the body is not exactly straight up/down, it has an azimuth.
        # Invert the angle to produce degrees eastward from north.
        az = math.degrees(-math.atan2(pw, pn))
        if az < 0:
            az += 360
    else:
        # The body is straight up/down, so it does not have an azimuth.
        # Report an arbitrary but reasonable value.
        az = 0.0

    # zd = the angle of the body away from the observer's zenith.
    zd = math.degrees(math.atan2(proj, pz))
    hor_ra = ra
    hor_dec = dec

    if refraction != Refraction.Airless:
        zd0 = zd
        refr = RefractionAngle(refraction, 90.0 - zd)
        zd -= refr
        if refr > 0.0 and zd > 3.0e-4:
            zdrad = math.radians(zd)
            sinzd = math.sin(zdrad)
            coszd = math.cos(zdrad)
            zd0rad = math.radians(zd0)
            sinzd0 = math.sin(zd0rad)
            coszd0 = math.cos(zd0rad)

            pr = [(((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd) for j in range(3)]
            proj = math.sqrt(pr[0]*pr[0] + pr[1]*pr[1])
            if proj > 0:
                hor_ra = _RAD2HOUR * math.atan2(pr[1], pr[0])
                if hor_ra < 0:
                    hor_ra += 24
            else:
                hor_ra = 0
            hor_dec = math.degrees(math.atan2(pr[2], proj))

    return HorizontalCoordinates(az, 90.0 - zd, hor_ra, hor_dec)

def RefractionAngle(refraction, altitude):
    """Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.

    Given an altitude angle and a refraction option, calculates
    the amount of "lift" caused by atmospheric refraction.
    This is the number of degrees higher in the sky an object appears
    due to lensing of the Earth's atmosphere.

    Parameters
    ----------
    refraction : Refraction
        The option for selecting whether to correct for atmospheric lensing.
        If `Refraction.Normal`, a well-behaved refraction model is used.
        If `Refraction.Airless`, no refraction correct is performed.
        `Refraction.JplHorizons` is used only for compatibility testing
        with the JPL Horizons online tool.
        Any other value raises an exception.
    altitude : float
        The number of degrees above (positive) or below (negative) the
        horizon an object is, before being corrected for refraction.

    Returns
    -------
    float
        The number of additional degrees of altitude an object appears
        to have, due to atmospheric refraction, depending on the
        option selected by the `refraction` parameter.
    """
    if altitude < -90.0 or altitude > +90.0:
        return 0.0      # No attempt to correct an invalid altitude

    if refraction in (Refraction.Normal, Refraction.JplHorizons):
        # http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        # JPL Horizons says it uses refraction algorithm from
        # Meeus "Astronomical Algorithms", 1991, p. 101-102.
        # I found the following Go implementation:
        # https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        # This is a translation from the function "Saemundsson" there.
        # I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        # This is important because the 'refr' formula below goes crazy near hd = -5.11.
        hd = max(altitude, -1.0)
        refr = (1.02 / math.tan(math.radians((hd+10.3/(hd+5.11))))) / 60.0

        if refraction == Refraction.Normal and altitude < -1.0:
            # In "normal" mode we gradually reduce refraction toward the nadir
            # so that we never get an altitude angle less than -90 degrees.
            # When horizon angle is -1 degrees, the factor is exactly 1.
            # As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            refr *= (altitude + 90.0) / 89.0
    elif refraction == Refraction.Airless:
        # The caller does not want refraction correction.
        refr = 0.0
    else:
        raise Error('Inalid refraction option')
    return refr

def InverseRefractionAngle(refraction, bent_altitude):
    """Calculates the inverse of an atmospheric refraction angle.

    Given an observed altitude angle that includes atmospheric refraction,
    calculate the negative angular correction to obtain the unrefracted
    altitude. This is useful for cases where observed horizontal
    coordinates are to be converted to another orientation system,
    but refraction first must be removed from the observed position.

    Parameters
    ----------
    refraction : Refraction
        `Refraction.Normal` - corrects for atmospheric refraction (recommended).
        `Refraction.Airless` - no correction is performed.
        `Refraction.JplHorizons` - For JPL Horizons compatibility testing only.
    bent_altitude : float
        The apparent altitude that includes atmospheric refraction.

    Returns
    -------
    float
        The angular adjustment in degrees, to be added to the
        altitude angle to correct for atmospheric lensing.
        This will be less than or equal to zero.
    """
    if bent_altitude < -90.0 or bent_altitude > +90.0:
        return 0.0      # No attempt to correct an invalid altitude
    # Find the pre-adjusted altitude whose refraction correction leads to 'altitude'.
    altitude = bent_altitude - RefractionAngle(refraction, bent_altitude)
    while True:
        # See how close we got. Keep iterating until the solution converges.
        diff = (altitude + RefractionAngle(refraction, altitude)) - bent_altitude
        if abs(diff) < 1.0e-14:
            return altitude - bent_altitude
        altitude -= diff

class EclipticCoordinates:
    """Ecliptic angular and Cartesian coordinates.

    Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
    oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

    Attributes
    ----------
    vec : Vector
        Ecliptic cartesian vector with the following components:
        x: in the direction of the equinox along the ecliptic plane.
        y: Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox.
        z: Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north.
    elat : float
        Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
    elon : float
        Longitude in degrees around the ecliptic plane prograde from the equinox.
    """
    def __init__(self, vec, elat, elon):
        self.vec = vec
        self.elat = elat
        self.elon = elon

def _RotateEquatorialToEcliptic(pos, obliq_radians, time):
    cos_ob = math.cos(obliq_radians)
    sin_ob = math.sin(obliq_radians)
    ex = +pos[0]
    ey = +pos[1]*cos_ob + pos[2]*sin_ob
    ez = -pos[1]*sin_ob + pos[2]*cos_ob
    xyproj = math.sqrt(ex*ex + ey*ey)
    if xyproj > 0.0:
        elon = math.degrees(math.atan2(ey, ex))
        if elon < 0.0:
            elon += 360.0
    else:
        elon = 0.0
    elat = math.degrees(math.atan2(ez, xyproj))
    vec = Vector(ex, ey, ez, time)
    return EclipticCoordinates(vec, elat, elon)

def SunPosition(time):
    """Calculates geocentric ecliptic coordinates for the Sun.

    This function calculates the position of the Sun as seen from the Earth.
    The returned value includes both Cartesian and spherical coordinates.
    The x-coordinate and longitude values in the returned object are based
    on the *true equinox of date*: one of two points in the sky where the instantaneous
    plane of the Earth's equator at the given date and time (the *equatorial plane*)
    intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
    By convention, the apparent location of the Sun at the March equinox is chosen
    as the longitude origin and x-axis direction, instead of the one for September.

    `SunPosition` corrects for precession and nutation of the Earth's axis
    in order to obtain the exact equatorial plane at the given time.

    This function can be used for calculating changes of seasons: equinoxes and solstices.
    In fact, the function #Seasons does use this function for that purpose.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Sun's position.

    Returns
    -------
    EclipticCoordinates
        The ecliptic coordinates of the Sun using the Earth's true equator of date.
    """
    # Correct for light travel time from the Sun.
    # Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes!
    adjusted_time = time.AddDays(-1.0 / C_AUDAY)
    earth2000 = _CalcEarth(adjusted_time)
    sun2000 = [-earth2000.x, -earth2000.y, -earth2000.z]

    # Convert to equatorial Cartesian coordinates of date.
    stemp = _precession(sun2000, adjusted_time, _PrecessDir.From2000)
    sun_ofdate = _nutation(stemp, adjusted_time, _PrecessDir.From2000)

    # Convert equatorial coordinates to ecliptic coordinates.
    true_obliq = math.radians(adjusted_time._etilt().tobl)
    return _RotateEquatorialToEcliptic(sun_ofdate, true_obliq, time)

def Ecliptic(equ):
    """Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.

    Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
    on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
    which are relative to the plane of the Earth's orbit around the Sun.

    Parameters
    ----------
    equ : Equatorial
        Equatorial coordinates in the J2000 frame of reference.

    Returns
    -------
    EclipticCoordinates
        Ecliptic coordinates in the J2000 frame of reference.
    """
    # Based on NOVAS functions equ2ecl() and equ2ecl_vec().
    ob2000 = 0.40909260059599012   # mean obliquity of the J2000 ecliptic in radians
    return _RotateEquatorialToEcliptic([equ.x, equ.y, equ.z], ob2000, equ.t)

def EclipticLongitude(body, time):
    """Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.

    This function calculates the angle around the plane of the Earth's orbit
    of a celestial body, as seen from the center of the Sun.
    The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
    in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).

    Parameters
    ----------
    body : Body
        A body other than the Sun.

    time : Time
        The date and time at which the body's ecliptic longitude is to be calculated.

    Returns
    -------
    float
        An angular value in degrees indicating the ecliptic longitude of the body.
    """
    if body == Body.Sun:
        raise InvalidBodyError()
    hv = HelioVector(body, time)
    eclip = Ecliptic(hv)
    return eclip.elon

def AngleFromSun(body, time):
    """Returns the angle between the given body and the Sun, as seen from the Earth.

    This function calculates the angular separation between the given body and the Sun,
    as seen from the center of the Earth. This angle is helpful for determining how
    easy it is to see the body away from the glare of the Sun.

    Parameters
    ----------
    body : Body
        The celestial body whose angle from the Sun is to be measured.
        Not allowed to be `Body.Earth`.
    time : Time
        The time at which the observation is made.

    Returns
    -------
    float
        A numeric value indicating the angle in degrees between the Sun
        and the specified body as seen from the center of the Earth.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    sv = GeoVector(Body.Sun, time, True)
    bv = GeoVector(body, time, True)
    return AngleBetween(sv, bv)

def PairLongitude(body1, body2, time):
    """Returns one body's ecliptic longitude with respect to another, as seen from the Earth.

    This function determines where one body appears around the ecliptic plane
    (the plane of the Earth's orbit around the Sun) as seen from the Earth,
    relative to the another body's apparent position.
    The function returns an angle in the half-open range [0, 360) degrees.
    The value is the ecliptic longitude of `body1` relative to the ecliptic
    longitude of `body2`.

    The angle is 0 when the two bodies are at the same ecliptic longitude
    as seen from the Earth. The angle increases in the prograde direction
    (the direction that the planets orbit the Sun and the Moon orbits the Earth).

    When the angle is 180 degrees, it means the two bodies appear on opposite sides
    of the sky for an Earthly observer.

    Neither `body1` nor `body2` is allowed to be `Body.Earth`.
    If this happens, the function throws an exception.

    Parameters
    ----------
    body1 : Body
        The first body, whose longitude is to be found relative to the second body.
    body2 : Body
        The second body, relative to which the longitude of the first body is to be found.
    time : Time
        The date and time of the observation.

    Returns
    -------
    float
        An angle in degrees in the range [0, 360).
    """
    if body1 == Body.Earth or body2 == Body.Earth:
        raise EarthNotAllowedError()
    vector1 = GeoVector(body1, time, False)
    eclip1 = Ecliptic(vector1)
    vector2 = GeoVector(body2, time, False)
    eclip2 = Ecliptic(vector2)
    return _NormalizeLongitude(eclip1.elon - eclip2.elon)

class ElongationEvent:
    """Contains information about the visibility of a celestial body at a given date and time.

    See the #Elongation function for more detailed information about the members of this class.
    See also #SearchMaxElongation for how to search for maximum elongation events.

    Attributes
    ----------
    time : Time
        The date and time of the observation.
    visibility : Visibility
        Whether the body is best seen in the morning or the evening.
    elongation : float
        The angle in degrees between the body and the Sun, as seen from the Earth.
    ecliptic_separation : float
        The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.
    """
    def __init__(self, time, visibility, elongation, ecliptic_separation):
        self.time = time
        self.visibility = visibility
        self.elongation = elongation
        self.ecliptic_separation = ecliptic_separation

@enum.unique
class Visibility(enum.Enum):
    """Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.

    Values
    ------
    Morning : The body is best visible in the morning, before sunrise.
    Evening : The body is best visible in the evening, after sunset.
    """
    Morning = 0
    Evening = 1

def Elongation(body, time):
    """Determines visibility of a celestial body relative to the Sun, as seen from the Earth.

    This function returns an #ElongationEvent object, which provides the following
    information about the given celestial body at the given time:

    - `visibility` is an enumerated type that specifies whether the body is more
      easily seen in the morning before sunrise, or in the evening after sunset.

    - `elongation` is the angle in degrees between two vectors: one from the center
      of the Earth to the center of the Sun, the other from the center of the Earth
      to the center of the specified body. This angle indicates how far away the body
      is from the glare of the Sun. The elongation angle is always in the range [0, 180].

    - `ecliptic_separation` is the absolute value of the difference between the body's
      ecliptic longitude and the Sun's ecliptic longitude, both as seen from the center
      of the Earth. This angle measures around the plane of the Earth's orbit, and ignores
      how far above or below that plane the body is.
      The ecliptic separation is measured in degrees and is always in the range [0, 180].

    Parameters
    ----------
    body : Body
        The celestial body whose visibility is to be calculated.

    time : Time
        The date and time of the observation.

    Returns
    -------
    ElongationEvent
    """
    angle = PairLongitude(body, Body.Sun, time)
    if angle > 180.0:
        visibility = Visibility.Morning
        esep = 360.0 - angle
    else:
        visibility = Visibility.Evening
        esep = angle
    angle = AngleFromSun(body, time)
    return ElongationEvent(time, visibility, angle, esep)

def _rlon_offset(body, time, direction, targetRelLon):
    plon = EclipticLongitude(body, time)
    elon = EclipticLongitude(Body.Earth, time)
    diff = direction * (elon - plon)
    return _LongitudeOffset(diff - targetRelLon)

def SearchRelativeLongitude(body, targetRelLon, startTime):
    """Searches for when the Earth and another planet are separated by a certain ecliptic longitude.

    Searches for the time when the Earth and another planet are separated by a specified angle
    in ecliptic longitude, as seen from the Sun.

    A relative longitude is the angle between two bodies measured in the plane of the
    Earth's orbit (the ecliptic plane). The distance of the bodies above or below the ecliptic
    plane is ignored. If you imagine the shadow of the body cast onto the ecliptic plane,
    and the angle measured around that plane from one body to the other in the direction
    the planets orbit the Sun, you will get an angle somewhere between 0 and 360 degrees.
    This is the relative longitude.

    Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
    this function searches for the next time that the relative longitude measured from the
    planet to the Earth is `targetRelLon`.

    Certain astronomical events are defined in terms of relative longitude between
    the Earth and another planet:

    - When the relative longitude is 0 degrees, it means both planets are in the same
      direction from the Sun. For planets that orbit closer to the Sun (Mercury and Venus),
      this is known as *inferior conjunction*, a time when the other planet becomes very
      difficult to see because of being lost in the Sun's glare.
      (The only exception is in the rare event of a transit, when we see the silhouette
      of the planet passing between the Earth and the Sun.)

    - When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
      this is known as *opposition*. Opposition is when the planet is closest to the Earth,
      and also when it is visible for most of the night, so it is considered the best time
      to observe the planet.

    - When the relative longitude is 180 degrees, it means the other planet is on the opposite
      side of the Sun from the Earth.  This is called *superior conjunction*.  Like inferior
      conjunction, the planet is very difficult to see from the Earth.
      Superior conjunction is possible for any planet other than the Earth.

    Parameters
    ----------
    body : Body
        A planet other than the Earth. If `body` is not a planet, or if it is `Body.Earth`, an error occurs.
    targetRelLon : float
        The desired relative longitude, expressed in degrees. Must be in the range [0, 360).
    startTime : Time
        The date and time at which to begin the search.

    Returns
    -------
    Time
        The date and time of the relative longitude event.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    if body in (Body.Moon, Body.Sun):
        raise InvalidBodyError()
    syn = _SynodicPeriod(body)
    direction = +1 if _IsSuperiorPlanet(body) else -1
    # Iterate until we converge on the desired event.
    # Calculate the error angle, which will be a negative number of degrees,
    # meaning we are "behind" the target relative longitude.
    error_angle = _rlon_offset(body, startTime, direction, targetRelLon)
    if error_angle > 0.0:
        error_angle -= 360.0    # force searching forward in time
    time = startTime
    iter_count = 0
    while iter_count < 100:
        # Estimate how many days in the future (positive) or past (negative)
        # we have to go to get closer to the target relative longitude.
        day_adjust = (-error_angle/360.0) * syn
        time = time.AddDays(day_adjust)
        if abs(day_adjust) * _SECONDS_PER_DAY < 1.0:
            return time
        prev_angle = error_angle
        error_angle = _rlon_offset(body, time, direction, targetRelLon)
        if abs(prev_angle) < 30.0 and prev_angle != error_angle:
            # Improve convergence for Mercury/Mars (eccentric orbits)
            # by adjusting the synodic period to more closely match the
            # variable speed of both planets in this part of their respective orbits.
            ratio = prev_angle / (prev_angle - error_angle)
            if 0.5 < ratio < 2.0:
                syn *= ratio
        iter_count += 1
    raise NoConvergeError()

def _neg_elong_slope(body, time):
    dt = 0.1
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    e1 = AngleFromSun(body, t1)
    e2 = AngleFromSun(body, t2)
    return (e1 - e2)/dt

def SearchMaxElongation(body, startTime):
    """Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.

    Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
    Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
    The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
    a telescope without atmospheric interference, are when these planets reach maximum elongation.
    These are events where the planets reach the maximum angle from the Sun as seen from the Earth.

    This function solves for those times, reporting the next maximum elongation event's date and time,
    the elongation value itself, the relative longitude with the Sun, and whether the planet is best
    observed in the morning or evening. See #ElongationEvent for more details about the returned object.

    Parameters
    ----------
    body : Body
        Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception.
        To find the best viewing opportunities for planets farther from the Sun than the
        Earth is (Mars through Pluto), use #SearchRelativeLongitude to find the next opposition event.
    startTime : Time
        The date and time at which to begin the search. The maximum elongation event
        found will always be the first one that occurs after this date and time.

    Returns
    -------
    ElongationEvent
    """
    if body == Body.Mercury:
        s1 = 50.0
        s2 = 85.0
    elif body == Body.Venus:
        s1 = 40.0
        s2 = 50.0
    else:
        raise InvalidBodyError()
    syn = _SynodicPeriod(body)
    iter_count = 1
    while iter_count <= 2:
        plon = EclipticLongitude(body, startTime)
        elon = EclipticLongitude(Body.Earth, startTime)
        rlon = _LongitudeOffset(plon - elon)    # clamp to (-180, +180]

        # The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        # because there is a cusp there that causes a discontinuity in the derivative.
        # So we need to guard against searching near such times.
        if -s1 <= rlon < +s1:
            # Seek to the window [+s1, +s2].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +s1
            # Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +s2
        elif rlon > +s2 or rlon < -s2:
            # Seek to the next search window at [-s2, -s1].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -s2
            # Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -s1
        elif rlon >= 0.0:
            # rlon must be in the middle of the window [+s1, +s2].
            # Search BACKWARD for the time t1 when rel lon = +s1.
            adjust_days = -syn / 4.0
            rlon_lo = +s1
            rlon_hi = +s2
            # Search forward from t1 to find t2 such that rel lon = +s2.
        else:
            # rlon must be in the middle of the window [-s2, -s1].
            # Search BACKWARD for the time t1 when rel lon = -s2.
            adjust_days = -syn / 4.0
            rlon_lo = -s2
            # Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -s1

        t_start = startTime.AddDays(adjust_days)
        t1 = SearchRelativeLongitude(body, rlon_lo, t_start)
        if t1 is None:
            return None

        t2 = SearchRelativeLongitude(body, rlon_hi, t1)
        if t2 is None:
            return None

        # Now we have a time range [t1,t2] that brackets a maximum elongation event.
        # Confirm the bracketing.
        m1 = _neg_elong_slope(body, t1)
        if m1 >= 0.0:
            raise InternalError()   # there is a bug in the bracketing algorithm!

        m2 = _neg_elong_slope(body, t2)
        if m2 <= 0.0:
            raise InternalError()   # there is a bug in the bracketing algorithm!

        # Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        tx = Search(_neg_elong_slope, body, t1, t2, 10.0)
        if tx is None:
            return None

        if tx.tt >= startTime.tt:
            return Elongation(body, tx)

        # This event is in the past (earlier than startTime).
        # We need to search forward from t2 to find the next possible window.
        # We never need to search more than twice.
        startTime = t2.AddDays(1.0)
        iter_count += 1


def _sun_offset(targetLon, time):
    ecl = SunPosition(time)
    return _LongitudeOffset(ecl.elon - targetLon)

def SearchSunLongitude(targetLon, startTime, limitDays):
    """Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.

    This function finds the moment in time, if any exists in the given time window,
    that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

    This function can be used to determine equinoxes and solstices.
    However, it is usually more convenient and efficient to call #Seasons
    to calculate all equinoxes and solstices for a given calendar year.

    The function searches the window of time specified by `startTime` and `startTime+limitDays`.
    The search will return `None` if the Sun never reaches the longitude `targetLon` or
    if the window is so large that the longitude ranges more than 180 degrees within it.
    It is recommended to keep the window smaller than 10 days when possible.

    Parameters
    ----------

    targetLon : float
         The desired ecliptic longitude in degrees, relative to the true equinox of date.
         This may be any value in the range [0, 360), although certain values have
         conventional meanings:
         0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice.
    startTime : Time
         The date and time for starting the search for the desired longitude event.
    limitDays : float
         The real-valued number of days, which when added to `startTime`, limits the
         range of time over which the search looks.
         It is recommended to keep this value between 1 and 10 days.
         See remarks above for more details.

    Returns
    -------
    Time or `None`
    """
    t2 = startTime.AddDays(limitDays)
    return Search(_sun_offset, targetLon, startTime, t2, 1.0)

def MoonPhase(time):
    """Returns the Moon's phase as an angle from 0 to 360 degrees.

    This function determines the phase of the Moon using its apparent
    ecliptic longitude relative to the Sun, as seen from the center of the Earth.
    Certain values of the angle have conventional definitions:

    - 0 = new moon
    - 90 = first quarter
    - 180 = full moon
    - 270 = third quarter

    Parameters
    ----------
    time : Time
         The date and time of the observation.

    Returns
    -------
    float
    """
    return PairLongitude(Body.Moon, Body.Sun, time)

def _moon_offset(targetLon, time):
    angle = MoonPhase(time)
    return _LongitudeOffset(angle - targetLon)

def SearchMoonPhase(targetLon, startTime, limitDays):
    """Searches for the time that the Moon reaches a specified phase.

    Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
    longitude with respect to the Sun's geocentric ecliptic longitude.
    When the Moon and the Sun have the same longitude, that is defined as a new moon.
    When their longitudes are 180 degrees apart, that is defined as a full moon.

    This function searches for any value of the lunar phase expressed as an
    angle in degrees in the range [0, 360).

    If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
    it is much easier to call the functions #SearchMoonQuarter and #NextMoonQuarter.
    This function is useful for finding general phase angles outside those four quarters.

    Parameters
    ----------
    targetLon : float
         The difference in geocentric longitude between the Sun and Moon
         that specifies the lunar phase being sought. This can be any value
         in the range [0, 360).  Certain values have conventional names:
         0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter.
    startTime : Time
         The beginning of the time window in which to search for the Moon reaching the specified phase.
    limitDays : float
         The number of days after `startTime` that limits the time window for the search.

    Returns
    -------
    Time or `None`
    """
    # To avoid discontinuities in the _moon_offset function causing problems,
    # we need to approximate when that function will next return 0.
    # We probe it with the start time and take advantage of the fact
    # that every lunar phase repeats roughly every 29.5 days.
    # There is a surprising uncertainty in the quarter timing,
    # due to the eccentricity of the moon's orbit.
    # I have seen more than 0.9 days away from the simple prediction.
    # To be safe, we take the predicted time of the event and search
    # +/-1.5 days around it (a 3-day wide window).
    # But we must return None if the final result goes beyond limitDays after startTime.
    uncertainty = 1.5
    ya = _moon_offset(targetLon, startTime)
    if ya > 0.0:
        ya -= 360.0     # force searching forward in time, not backward
    est_dt = -(_MEAN_SYNODIC_MONTH * ya) / 360.0
    dt1 = est_dt - uncertainty
    if dt1 > limitDays:
        return None     # not possible for moon phase to occur within the specified window
    dt2 = min(limitDays, est_dt + uncertainty)
    t1 = startTime.AddDays(dt1)
    t2 = startTime.AddDays(dt2)
    return Search(_moon_offset, targetLon, t1, t2, 1.0)

class MoonQuarter:
    """A lunar quarter event along with its date and time.

    An object of this type represents one of the four major
    lunar phases that appear on calendars:
    new moon, first quarter, full moon, or third quarter.
    Along with the `quarter` attribute that specifies the
    type of quarter, it contains a `time` field that indicates
    when the lunar quarter event happens.

    Attributes
    ----------
    quarter : int
        0=new moon, 1=first quarter, 2=full moon, 3=third quarter.
    time : Time
        The date and time of the lunar quarter.
    """
    def __init__(self, quarter, time):
        self.quarter = quarter
        self.time = time

def SearchMoonQuarter(startTime):
    """Finds the first lunar quarter after the specified date and time.

    A lunar quarter is one of the following four lunar phase events:
    new moon, first quarter, full moon, third quarter.
    This function finds the lunar quarter that happens soonest
    after the specified date and time.

    To continue iterating through consecutive lunar quarters, call this function once,
    followed by calls to #NextMoonQuarter as many times as desired.

    Parameters
    ----------
    startTime : Time
        The date and time at which to start the search.

    Returns
    -------
    MoonQuarter
    """
    angle = MoonPhase(startTime)
    quarter = int(1 + math.floor(angle / 90.0)) % 4
    time = SearchMoonPhase(90.0 * quarter, startTime, 10.0)
    if time is None:
        # The search should never fail. We should always find another lunar quarter.
        raise InternalError()
    return MoonQuarter(quarter, time)

def NextMoonQuarter(mq):
    """Continues searching for lunar quarters from a previous search.

    After calling #SearchMoonQuarter, this function can be called
    one or more times to continue finding consecutive lunar quarters.
    This function finds the next consecutive moon quarter event after
    the one passed in as the parameter `mq`.

    Parameters
    ----------
    mq : MoonQuarter
        A value returned by a prior call to #SearchMoonQuarter or #NextMoonQuarter.

    Returns
    -------
    MoonQuarter
    """
    # Skip 6 days past the previous found moon quarter to find the next one.
    # This is less than the minimum possible increment.
    # So far I have seen the interval well contained by the range (6.5, 8.3) days.
    time = mq.time.AddDays(6.0)
    next_mq = SearchMoonQuarter(time)
    # Verify that we found the expected moon quarter.
    if next_mq.quarter != (1 + mq.quarter) % 4:
        raise InternalError()
    return next_mq


class IlluminationInfo:
    """Information about the brightness and illuminated shape of a celestial body.

    Returned by functions #Illumination and #SearchPeakMagnitude
    to report the visual magnitude and illuminated fraction of a celestial
    body at a given date and time.

    Attributes
    ----------
    time : Time
        The date and time of the observation.
    mag : float
        The visual magnitude of the body. Smaller values are brighter.
    phase_angle : float
        The angle in degrees between the Sun and the Earth, as seen from the body.
        Indicates the body's phase as seen from the Earth.
    phase_fraction : float
        A value in the range [0.0, 1.0] indicating what fraction of the
        body's apparent disc is illuminated, as seen from the Earth.
    helio_dist : float
        The distance between the Sun and the body at the observation time, in AU.
    ring_tilt : float
        For Saturn, the tilt angle in degrees of its rings as seen from Earth.
        When the `ring_tilt` is very close to 0, it means the rings are edge-on
        as seen from observers on the Earth, and are thus very difficult to see.
        For bodies other than Saturn, `ring_tilt` is `None`.
    """
    def __init__(self, time, mag, phase, helio_dist, geo_dist, gc, hc, ring_tilt):
        self.time = time
        self.mag = mag
        self.phase_angle = phase
        self.phase_fraction = (1.0 + math.cos(math.radians(phase))) / 2.0
        self.helio_dist = helio_dist
        self.geo_dist = geo_dist
        self.gc = gc
        self.hc = hc
        self.ring_tilt = ring_tilt

def _MoonMagnitude(phase, helio_dist, geo_dist):
    # https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
    rad = math.radians(phase)
    mag = -12.717 + 1.49*abs(rad) + 0.0431*(rad**4)
    moon_mean_distance_au = 385000.6 / KM_PER_AU
    geo_au = geo_dist / moon_mean_distance_au
    mag += 5.0 * math.log10(helio_dist * geo_au)
    return mag

def _SaturnMagnitude(phase, helio_dist, geo_dist, gc, time):
    # Based on formulas by Paul Schlyter found here:
    # http://www.stjarnhimlen.se/comp/ppcomp.html#15

    # We must handle Saturn's rings as a major component of its visual magnitude.
    # Find geocentric ecliptic coordinates of Saturn.
    eclip = Ecliptic(gc)

    ir = math.radians(28.06)   # tilt of Saturn's rings to the ecliptic, in radians
    Nr = math.radians(169.51 + (3.82e-5 * time.tt))    # ascending node of Saturn's rings, in radians

    # Find tilt of Saturn's rings, as seen from Earth.
    lat = math.radians(eclip.elat)
    lon = math.radians(eclip.elon)
    tilt = math.asin(math.sin(lat)*math.cos(ir) - math.cos(lat)*math.sin(ir)*math.sin(lon-Nr))
    sin_tilt = math.sin(abs(tilt))

    mag = -9.0 + 0.044*phase
    mag += sin_tilt*(-2.6 + 1.2*sin_tilt)
    mag += 5.0 * math.log10(helio_dist * geo_dist)
    ring_tilt = math.degrees(tilt)
    return (mag, ring_tilt)

def _VisualMagnitude(body, phase, helio_dist, geo_dist):
    # For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
    c0 = c1 = c2 = c3 = 0
    if body == Body.Mercury:
        c0 = -0.60; c1 = +4.98; c2 = -4.88; c3 = +3.02
    elif body == Body.Venus:
        if phase < 163.6:
            c0 = -4.47; c1 = +1.03; c2 = +0.57; c3 = +0.13
        else:
            c0 = +0.98; c1 = -1.02
    elif body == Body.Mars:
        c0 = -1.52; c1 = +1.60
    elif body == Body.Jupiter:
        c0 = -9.40; c1 = +0.50
    elif body == Body.Uranus:
        c0 = -7.19; c1 = +0.25
    elif body == Body.Neptune:
        c0 = -6.87
    elif body == Body.Pluto:
        c0 = -1.00; c1 = +4.00
    else:
        raise InvalidBodyError()

    x = phase / 100.0
    mag = c0 + x*(c1 + x*(c2 + x*c3))
    mag += 5.0 * math.log10(helio_dist * geo_dist)
    return mag

def Illumination(body, time):
    """Finds visual magnitude, phase angle, and other illumination information about a celestial body.

    This function calculates information about how bright a celestial body appears from the Earth,
    reported as visual magnitude, which is a smaller (or even negative) number for brighter objects,
    and a larger number for dimmer objects.

    For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
    the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
    of the body appears illuminated as seen from the Earth. For example, when the phase angle is
    near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
    180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
    of 90 degrees means the body appears "half full".
    For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
    so it doesn't have a phase angle.

    When the body is Saturn, the returned object contains a field `ring_tilt` that holds
    the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
    the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
    0 for all bodies other than Saturn.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    time : Time
        The date and time of the observation.

    Returns
    -------
    IlluminationInfo
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()
    earth = _CalcEarth(time)
    if body == Body.Sun:
        gc = Vector(-earth.x, -earth.y, -earth.z, time)
        hc = Vector(0.0, 0.0, 0.0, time)
        phase = 0.0     # placeholder value; the Sun does not have a phase angle.
    else:
        if body == Body.Moon:
            # For extra numeric precision, use geocentric moon formula directly.
            gc = GeoMoon(time)
            hc = Vector(earth.x + gc.x, earth.y + gc.y, earth.z + gc.z, time)
        else:
            # For planets, heliocentric vector is most direct to calculate.
            hc = HelioVector(body, time)
            gc = Vector(hc.x - earth.x, hc.y - earth.y, hc.z - earth.z, time)
        phase = AngleBetween(gc, hc)

    geo_dist = gc.Length()      # distance from body to center of Earth
    helio_dist = hc.Length()    # distance from body to center of Sun
    ring_tilt = None            # only reported for Saturn
    if body == Body.Sun:
        mag = -0.17 + 5.0*math.log10(geo_dist / _AU_PER_PARSEC)
    elif body == Body.Moon:
        mag = _MoonMagnitude(phase, helio_dist, geo_dist)
    elif body == Body.Saturn:
        mag, ring_tilt = _SaturnMagnitude(phase, helio_dist, geo_dist, gc, time)
    else:
        mag = _VisualMagnitude(body, phase, helio_dist, geo_dist)
    return IlluminationInfo(time, mag, phase, helio_dist, geo_dist, gc, hc, ring_tilt)

def _mag_slope(body, time):
    # The Search() function finds a transition from negative to positive values.
    # The derivative of magnitude y with respect to time t (dy/dt)
    # is negative as an object gets brighter, because the magnitude numbers
    # get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
    # dy/dt > 0.
    dt = 0.01
    t1 = time.AddDays(-dt/2)
    t2 = time.AddDays(+dt/2)
    y1 = Illumination(body, t1)
    y2 = Illumination(body, t2)
    return (y2.mag - y1.mag) / dt

def SearchPeakMagnitude(body, startTime):
    """Searches for the date and time Venus will next appear brightest as seen from the Earth.

    This function searches for the date and time Venus appears brightest as seen from the Earth.
    Currently only Venus is supported for the `body` parameter, though this could change in the future.
    Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see
    from the Earth, so peak magnitude events have little practical value for that planet.
    Planets other than Venus and Mercury reach peak magnitude at opposition, which can
    be found using #SearchRelativeLongitude.
    The Moon reaches peak magnitude at full moon, which can be found using
    #SearchMoonQuarter or #SearchMoonPhase.
    The Sun reaches peak magnitude at perihelion, which occurs each year in January.
    However, the difference is minor and has little practical value.

    Parameters
    ----------
    body : Body
        Currently only `Body.Venus` is allowed. Any other value results in an exception.
        See remarks above for more details.
    startTime : Time
        The date and time to start searching for the next peak magnitude event.

    Returns
    -------
    IlluminationInfo
    """
    # s1 and s2 are relative longitudes within which peak magnitude of Venus can occur.
    s1 = 10.0
    s2 = 30.0
    if body != Body.Venus:
        raise InvalidBodyError()

    iter_count = 1
    while iter_count <= 2:
        # Find current heliocentric relative longitude between the
        # inferior planet and the Earth.
        plon = EclipticLongitude(body, startTime)
        elon = EclipticLongitude(Body.Earth, startTime)
        rlon = _LongitudeOffset(plon - elon)
        # The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
        # because there is a cusp there that causes a discontinuity in the derivative.
        # So we need to guard against searching near such times.
        if -s1 <= rlon < +s1:
            # Seek to the window [+s1, +s2].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = +s1.
            rlon_lo = +s1
            # Search forward for the time t2 when rel lon = +s2.
            rlon_hi = +s2
        elif rlon >= +s2 or rlon < -s2:
            # Seek to the next search window at [-s2, -s1].
            adjust_days = 0.0
            # Search forward for the time t1 when rel lon = -s2.
            rlon_lo = -s2
            # Search forward for the time t2 when rel lon = -s1.
            rlon_hi = -s1
        elif rlon >= 0:
            # rlon must be in the middle of the window [+s1, +s2].
            # Search BACKWARD for the time t1 when rel lon = +s1.
            syn = _SynodicPeriod(body)
            adjust_days = -syn / 4
            rlon_lo = +s1
            # Search forward from t1 to find t2 such that rel lon = +s2.
            rlon_hi = +s2
        else:
            # rlon must be in the middle of the window [-s2, -s1].
            # Search BACKWARD for the time t1 when rel lon = -s2.
            syn = _SynodicPeriod(body)
            adjust_days = -syn / 4
            rlon_lo = -s2
            # Search forward from t1 to find t2 such that rel lon = -s1.
            rlon_hi = -s1

        t_start = startTime.AddDays(adjust_days)
        t1 = SearchRelativeLongitude(body, rlon_lo, t_start)
        t2 = SearchRelativeLongitude(body, rlon_hi, t1)

        # Now we have a time range [t1,t2] that brackets a maximum magnitude event.
        # Confirm the bracketing.
        m1 = _mag_slope(body, t1)
        if m1 >= 0.0:
            raise InternalError()

        m2 = _mag_slope(body, t2)
        if m2 <= 0.0:
            raise InternalError()

        # Use the generic search algorithm to home in on where the slope crosses from negative to positive.
        tx = Search(_mag_slope, body, t1, t2, 10.0)
        if tx is None:
            # The search should have found the ascending root in the interval [t1, t2].
            raise InternalError()

        if tx.tt >= startTime.tt:
            return Illumination(body, tx)

        # This event is in the past (earlier than startTime).
        # We need to search forward from t2 to find the next possible window.
        # We never need to search more than twice.
        startTime = t2.AddDays(1.0)
        iter_count += 1

    # We should have found the peak magnitude in at most 2 iterations.
    raise InternalError()

class HourAngleEvent:
    """Information about a celestial body crossing a specific hour angle.

    Returned by the function #SearchHourAngle to report information about
    a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

    Attributes
    ----------
    time : Time
        The date and time when the body crosses the specified hour angle.
    hor : HorizontalCoordinates
        Apparent coordinates of the body at the time it crosses the specified hour angle.
    """
    def __init__(self, time, hor):
        self.time = time
        self.hor = hor

def SearchHourAngle(body, observer, hourAngle, startTime):
    """Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.

    The *hour angle* of a celestial body indicates its position in the sky with respect
    to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
    The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
    The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
    to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
    the number of hours that have passed since the most recent time that the body has culminated,
    or reached its highest point.

    This function searches for the next time a celestial body reaches the given hour angle
    after the date and time specified by `startTime`.
    To find when a body culminates, pass 0 for `hourAngle`.
    To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.

    Note that, especially close to the Earth's poles, a body as seen on a given day
    may always be above the horizon or always below the horizon, so the caller cannot
    assume that a culminating object is visible nor that an object is below the horizon
    at its minimum altitude.

    On success, the function reports the date and time, along with the horizontal coordinates
    of the body at that time, as seen by the given observer.

    Parameters
    ----------
    body : Body
         The celestial body, which can the Sun, the Moon, or any planet other than the Earth.
    observer : Observer
         Indicates a location on or near the surface of the Earth where the observer is located.
    hourAngle : float
         An hour angle value in the range [0.0, 24.0) indicating the number of sidereal hours after the
         body's most recent culmination.
    startTime : Time
         The date and time at which to start the search.

    Returns
    -------
    HourAngleEvent
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if hourAngle < 0.0 or hourAngle >= 24.0:
        raise Error('Invalid hour angle.')

    iter_count = 0
    time = startTime
    while True:
        iter_count += 1
        # Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
        gast = _sidereal_time(time)
        ofdate = Equator(body, time, observer, True, True)

        # Calculate the adjustment needed in sidereal time to bring
        # the hour angle to the desired value.
        delta_sidereal_hours = math.fmod(((hourAngle + ofdate.ra - observer.longitude/15) - gast), 24.0)
        if iter_count == 1:
            # On the first iteration, always search forward in time.
            if delta_sidereal_hours < 0.0:
                delta_sidereal_hours += 24.0
        else:
            # On subsequent iterations, we make the smallest possible adjustment,
            # either forward or backward in time.
            if delta_sidereal_hours < -12.0:
                delta_sidereal_hours += 24.0
            elif delta_sidereal_hours > +12.0:
                delta_sidereal_hours -= 24.0

        # If the error is tolerable (less than 0.1 seconds), stop searching.
        if abs(delta_sidereal_hours) * 3600.0 < 0.1:
            hor = Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.Normal)
            return HourAngleEvent(time, hor)

        # We need to loop another time to get more accuracy.
        # Update the terrestrial time (in solar days) adjusting by sidereal time.
        delta_days = (delta_sidereal_hours / 24.0) * _SOLAR_DAYS_PER_SIDEREAL_DAY
        time = time.AddDays(delta_days)

@enum.unique
class Direction(enum.Enum):
    """Indicates whether a body is rising above or setting below the horizon.

    Specifies the direction of a rising or setting event for a body.
    For example, `Direction.Rise` is used to find sunrise times,
    and `Direction.Set` is used to find sunset times.

    Values
    ------
    Rise:   First appearance of a body as it rises above the horizon.
    Set:    Last appearance of a body as it sinks below the horizon.
    """
    Rise = +1
    Set  = -1


def _InternalSearchAltitude(body, observer, direction, startTime, limitDays, altitude_error_func, altitude_error_context):
    if direction == Direction.Rise:
        ha_before = 12.0    # minimum altitude (bottom) happens BEFORE the body rises.
        ha_after  =  0.0    # maximum altitude (culmination) happens AFTER the body rises.
    elif direction == Direction.Set:
        ha_before =  0.0    # culmination happens BEFORE the body sets.
        ha_after  = 12.0    # bottom happens AFTER the body sets.
    else:
        raise Error('Invalid value for direction parameter')

    # See if the body is currently above/below the horizon.
    # If we are looking for next rise time and the body is below the horizon,
    # we use the current time as the lower time bound and the next culmination
    # as the upper bound.
    # If the body is above the horizon, we search for the next bottom and use it
    # as the lower bound and the next culmination after that bottom as the upper bound.
    # The same logic applies for finding set times, only we swap the hour angles.
    time_start = startTime
    alt_before = altitude_error_func(altitude_error_context, time_start)
    if alt_before > 0.0:
        # We are past the sought event, so we have to wait for the next "before" event (culm/bottom).
        evt_before = SearchHourAngle(body, observer, ha_before, time_start)
        time_before = evt_before.time
        alt_before = altitude_error_func(altitude_error_context, time_before)
    else:
        # We are before or at the sought ebvent, so we find the next "after" event (bottom/culm),
        # and use the current time as the "before" event.
        time_before = time_start

    evt_after = SearchHourAngle(body, observer, ha_after, time_before)
    alt_after = altitude_error_func(altitude_error_context, evt_after.time)

    while True:
        if alt_before <= 0.0 and alt_after > 0.0:
            # Search between the "before time" and the "after time" for the desired event.
            event_time = Search(altitude_error_func, altitude_error_context, time_before, evt_after.time, 1.0)
            if event_time is not None:
                return event_time
        # We didn't find the desired event, so use the "after" time to find the next "before" event.
        evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time)
        evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time)
        if evt_before.time.ut >= time_start.ut + limitDays:
            return None
        time_before = evt_before.time
        alt_before = altitude_error_func(altitude_error_context, evt_before.time)
        alt_after = altitude_error_func(altitude_error_context, evt_after.time)


class _peak_altitude_context:
    def __init__(self, body, direction, observer, body_radius_au):
        self.body = body
        self.direction = direction
        self.observer = observer
        self.body_radius_au = body_radius_au


def _peak_altitude(context, time):
    # Return the angular altitude above or below the horizon
    # of the highest part (the peak) of the given object.
    # This is defined as the apparent altitude of the center of the body plus
    # the body's angular radius.
    # The 'direction' parameter controls whether the angle is measured
    # positive above the horizon or positive below the horizon,
    # depending on whether the caller wants rise times or set times, respectively.

    ofdate = Equator(context.body, time, context.observer, True, True)

    # We calculate altitude without refraction, then add fixed refraction near the horizon.
    # This gives us the time of rise/set without the extra work.
    hor = Horizon(time, context.observer, ofdate.ra, ofdate.dec, Refraction.Airless)
    alt = hor.altitude + math.degrees(context.body_radius_au / ofdate.dist)
    return context.direction.value * (alt + _REFRACTION_NEAR_HORIZON)


def SearchRiseSet(body, observer, direction, startTime, limitDays):
    """Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.

    This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
    Rise time is when the body first starts to be visible above the horizon.
    For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
    Set time is the moment when the body appears to vanish below the horizon.

    This function corrects for typical atmospheric refraction, which causes celestial
    bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
    It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).

    Note that rise or set may not occur in every 24 hour period.
    For example, near the Earth's poles, there are long periods of time where
    the Sun stays below the horizon, never rising.
    Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
    This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
    significant amount during each rotation of the Earth.
    Therefore callers must not assume that the function will always succeed.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    observer : Observer
        The location where observation takes place.
    direction : Direction
        Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time.
    startTime : Time
        The date and time at which to start the search.
    limitDays : float
        Limits how many days to search for a rise or set time.
        To limit a rise or set time to the same day, you can use a value of 1 day.
        In cases where you want to find the next rise or set time no matter how far
        in the future (for example, for an observer near the south pole), you can pass
        in a larger value like 365.

    Returns
    -------
    Time or `None`
        If the rise or set time is found within the specified time window,
        this function returns that time. Otherwise, it returns `None`.
    """
    if body == Body.Earth:
        raise EarthNotAllowedError()

    if body == Body.Sun:
        body_radius = _SUN_RADIUS_AU
    elif body == Body.Moon:
        body_radius = _MOON_EQUATORIAL_RADIUS_AU
    else:
        body_radius = 0.0

    context = _peak_altitude_context(body, direction, observer, body_radius)
    return _InternalSearchAltitude(body, observer, direction, startTime, limitDays, _peak_altitude, context)


class _altitude_error_context:
    def __init__(self, body, direction, observer, altitude):
        self.body = body
        self.direction = direction
        self.observer = observer
        self.altitude = altitude

def _altitude_error_func(context, time):
    ofdate = Equator(context.body, time, context.observer, True, True)
    hor = Horizon(time, context.observer, ofdate.ra, ofdate.dec, Refraction.Airless)
    return context.direction.value * (hor.altitude - context.altitude)


def SearchAltitude(body, observer, direction, dateStart, limitDays, altitude):
    """Finds the next time a body reaches a given altitude.

    Finds when the given body ascends or descends through a given
    altitude angle, as seen by an observer at the specified location on the Earth.
    By using the appropriate combination of `direction` and `altitude` parameters,
    this function can be used to find when civil, nautical, or astronomical twilight
    begins (dawn) or ends (dusk).

    Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
    the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.

    Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
    To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.

    Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.

    Astronomical twilight uses -18 degrees as the `altitude` value.

    Parameters
    ----------
    body : Body
        The Sun, Moon, or any planet other than the Earth.
    observer : Observer
        The location where observation takes place.
    direction : Direction
        Either `Direction.Rise` to find an ascending altitude event
        or `Direction.Set` to find a descending altitude event.
    startTime : Time
        The date and time at which to start the search.
    limitDays : float
        The fractional number of days after `dateStart` that limits
        when the altitude event is to be found. Must be a positive number.

    Returns
    -------
    Time or `None`
        If the altitude event time is found within the specified time window,
        this function returns that time. Otherwise, it returns `None`.
    """
    if not (-90.0 <= altitude <= +90.0):
        raise Error('Invalid altitude: {}'.format(altitude))

    context = _altitude_error_context(body, direction, observer, altitude)
    return _InternalSearchAltitude(body, observer, direction, dateStart, limitDays, _altitude_error_func, context)

class SeasonInfo:
    """The dates and times of changes of season for a given calendar year.

    Call #Seasons to calculate this data structure for a given year.

    Attributes
    ----------
    mar_equinox : Time
        The date and time of the March equinox for the specified year.
    jun_solstice : Time
        The date and time of the June solstice for the specified year.
    sep_equinox : Time
        The date and time of the September equinox for the specified year.
    dec_solstice : Time
        The date and time of the December solstice for the specified year.
    """
    def __init__(self, mar_equinox, jun_solstice, sep_equinox, dec_solstice):
        self.mar_equinox = mar_equinox
        self.jun_solstice = jun_solstice
        self.sep_equinox = sep_equinox
        self.dec_solstice = dec_solstice

def _FindSeasonChange(targetLon, year, month, day):
    startTime = Time.Make(year, month, day, 0, 0, 0)
    time = SearchSunLongitude(targetLon, startTime, 4.0)
    if time is None:
        # We should always be able to find a season change.
        raise InternalError()
    return time

def Seasons(year):
    """Finds both equinoxes and both solstices for a given calendar year.

    The changes of seasons are defined by solstices and equinoxes.
    Given a calendar year number, this function calculates the
    March and September equinoxes and the June and December solstices.

    The equinoxes are the moments twice each year when the plane of the
    Earth's equator passes through the center of the Sun. In other words,
    the Sun's declination is zero at both equinoxes.
    The March equinox defines the beginning of spring in the northern hemisphere
    and the beginning of autumn in the southern hemisphere.
    The September equinox defines the beginning of autumn in the northern hemisphere
    and the beginning of spring in the southern hemisphere.

    The solstices are the moments twice each year when one of the Earth's poles
    is most tilted toward the Sun. More precisely, the Sun's declination reaches
    its minimum value at the December solstice, which defines the beginning of
    winter in the northern hemisphere and the beginning of summer in the southern
    hemisphere. The Sun's declination reaches its maximum value at the June solstice,
    which defines the beginning of summer in the northern hemisphere and the beginning
    of winter in the southern hemisphere.

    Parameters
    ----------
    year : int
        The calendar year number for which to calculate equinoxes and solstices.
        The value may be any integer, but only the years 1800 through 2100 have
        been validated for accuracy: unit testing against data from the
        United States Naval Observatory confirms that all equinoxes and solstices
        for that range of years are within 2 minutes of the correct time.

    Returns
    -------
    SeasonInfo
    """
    mar_equinox = _FindSeasonChange(0, year, 3, 19)
    jun_solstice = _FindSeasonChange(90, year, 6, 19)
    sep_equinox = _FindSeasonChange(180, year, 9, 21)
    dec_solstice = _FindSeasonChange(270, year, 12, 20)
    return SeasonInfo(mar_equinox, jun_solstice, sep_equinox, dec_solstice)

def _MoonDistance(time):
    return _CalcMoon(time).distance_au

def _moon_distance_slope(direction, time):
    dt = 0.001
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    dist1 = _MoonDistance(t1)
    dist2 = _MoonDistance(t2)
    return direction * (dist2 - dist1) / dt

@enum.unique
class ApsisKind(enum.Enum):
    """Represents whether a satellite is at a closest or farthest point in its orbit.

    An apsis is a point in a satellite's orbit that is closest to,
    or farthest from, the body it orbits (its primary).
    `ApsisKind` is an enumerated type that indicates which of these
    two cases applies to a particular apsis event.

    Values
    ------
    Pericenter: The satellite is at its closest point to its primary.
    Apocenter: The satellite is at its farthest point from its primary.
    Invalid: A placeholder for an undefined, unknown, or invalid apsis.
    """
    Pericenter = 0
    Apocenter  = 1
    Invalid    = 2

class Apsis:
    """An event where a satellite is closest to or farthest from the body it orbits.

    For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
    event where the orbiting body reaches its closest or farthest point from the primary body.
    The closest approach is called *pericenter* and the farthest point is *apocenter*.

    More specific terminology is common for particular orbiting bodies.
    The Moon's closest approach to the Earth is called *perigee* and its furthest
    point is called *apogee*. The closest approach of a planet to the Sun is called
    *perihelion* and the furthest point is called *aphelion*.

    This data structure is returned by #SearchLunarApsis and #NextLunarApsis
    to iterate through consecutive alternating perigees and apogees.

    Attributes
    ----------
    time : Time
        The date and time of the apsis.
    kind : ApsisKind
        Whether this is a pericenter or apocenter event.
    dist_au : float
        The distance between the centers of the bodies in astronomical units.
    dist_km : float
        The distance between the centers of the bodies in kilometers.
    """
    def __init__(self, time, kind, dist_au):
        self.time = time
        self.kind = kind
        self.dist_au = dist_au
        self.dist_km = dist_au * KM_PER_AU

def SearchLunarApsis(startTime):
    """Finds the time of the first lunar apogee or perigee after the given time.

    Given a date and time to start the search in `startTime`, this function finds
    the next date and time that the center of the Moon reaches the closest or
    farthest point in its orbit with respect to the center of the Earth, whichever
    comes first after `startTime`.  The return value (of type #Apsis) also
    contains an indicator of whether the event is apogee or perigee.

    The closest point is called *perigee* and the farthest point is called *apogee*.
    The word *apsis* refers to either event.

    To iterate through consecutive alternating perigee and apogee events,
    call #SearchLunarApsis once, then use the return value to call #NextLunarApsis.
    After that, keep feeding the previous return value from `NextLunarApsis` into
    another call of `NextLunarApsis` as many times as desired.

    Parameters
    ----------
    startTime : Time
        The date and time at which to start searching for the next perigee or apogee.

    Returns
    -------
    Apsis
    """
    increment = 5.0     # number of days to skip on each iteration
    t1 = startTime
    m1 = _moon_distance_slope(+1, t1)
    iter_count = 0
    while iter_count * increment < 2.0 * _MEAN_SYNODIC_MONTH:
        t2 = t1.AddDays(increment)
        m2 = _moon_distance_slope(+1, t2)
        if m1 * m2 <= 0.0:
            # There is a change of slope polarity within the time range [t1, t2].
            # Therefore this time range contains an apsis.
            # Figure out whether it is perigee or apogee.
            if m1 < 0.0 or m2 > 0.0:
                # We found a minimum-distance event: perigee.
                # Search the time range for the time when the slope goes from negative to positive.
                apsis_time = Search(_moon_distance_slope, +1, t1, t2, 1.0)
                kind = ApsisKind.Pericenter
            elif m1 > 0.0 or m2 < 0.0:
                # We found a maximum-distance event: apogee.
                # Search the time range for the time when the slope goes from positive to negative.
                apsis_time = Search(_moon_distance_slope, -1, t1, t2, 1.0)
                kind = ApsisKind.Apocenter
            else:
                # This should never happen. It should not be possible for both slopes to be zero.
                raise InternalError()

            if apsis_time is None:
                # We should always be able to find a lunar apsis!
                raise InternalError()

            dist = _MoonDistance(apsis_time)
            return Apsis(apsis_time, kind, dist)

        # We have not yet found a slope polarity change. Keep searching.
        t1 = t2
        m1 = m2
        iter_count += 1

    # It should not be possible to fail to find an apsis within 2 synodic months.
    raise InternalError()


def NextLunarApsis(apsis):
    """Finds the next lunar perigee or apogee in a series.

    This function requires an #Apsis value obtained from a call to
    #SearchLunarApsis or `NextLunarApsis`.
    Given an apogee event, this function finds the next perigee event,
    and vice versa.

    See #SearchLunarApsis for more details.

    Parameters
    ----------
    apsis : Apsis

    Returns
    -------
    Apsis
    """
    skip = 11.0     # number of days to skip to start looking for next apsis event
    time = apsis.time.AddDays(skip)
    next_apsis = SearchLunarApsis(time)
    # Verify that we found the opposite apsis from the previous one.
    if apsis.kind not in [ApsisKind.Apocenter, ApsisKind.Pericenter]:
        raise Error('Parameter "apsis" contains an invalid "kind" value.')
    if next_apsis.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite apsis kind
    return next_apsis


def _planet_distance_slope(context, time):
    (direction, body) = context
    dt = 0.001
    t1 = time.AddDays(-dt/2.0)
    t2 = time.AddDays(+dt/2.0)
    dist1 = HelioDistance(body, t1)
    dist2 = HelioDistance(body, t2)
    return direction * (dist2 - dist1) / dt


def SearchPlanetApsis(body, startTime):
    """Finds the next planet perihelion or aphelion, after a given time.

    Given a date and time to start the search in `startTime`, this function finds the
    next date and time that the center of the specified planet reaches the closest or farthest point
    in its orbit with respect to the center of the Sun, whichever comes first after `startTime`.

    The closest point is called *perihelion* and the farthest point is called *aphelion*.
    The word *apsis* refers to either event.

    To iterate through consecutive alternating perihelion and aphelion events,
    call `SearchPlanetApsis` once, then use the return value to call #NextPlanetApsis.
    After that, keep feeding the previous return value from `NextPlanetApsis`
    into another call of `NextPlanetApsis` as many times as desired.

    Parameters
    ----------
    body : Body
        The planet for which to find the next perihelion/aphelion event.
        Not allowed to be `Body.Sun` or `Body.Moon`.
    startTime : Time
        The date and time at which to start searching for the next perihelion or aphelion.

    Returns
    -------
    Apsis
    """
    if body == Body.Neptune or body == Body.Pluto:
        return _BruteSearchPlanetApsis(body, startTime)
    positive_slope = (+1.0, body)
    negative_slope = (-1.0, body)
    orbit_period_days = _PlanetOrbitalPeriod[body.value]
    increment = orbit_period_days / 6.0
    t1 = startTime
    m1 = _planet_distance_slope(positive_slope, t1)
    iter_count = 0
    while iter_count * increment < 2 * orbit_period_days:
        t2 = t1.AddDays(increment)
        m2 = _planet_distance_slope(positive_slope, t2)
        if m1 * m2 <= 0.0:
            # There is a change of slope polarity within the time range [t1, t2].
            # Therefore this time range contains an apsis.
            # Figure out whether it is perihelion or aphelion.
            if m1 < 0.0 or m2 > 0.0:
                # We found a minimum-distance event: perihelion.
                # Search the time range for the time when the slope goes from negative to positive.
                context = positive_slope
                kind = ApsisKind.Pericenter
            elif m1 > 0.0 or m2 < 0.0:
                # We found a maximum-distance event: aphelion.
                # Search the time range for the time when the slope goes from positive to negative.
                context = negative_slope
                kind = ApsisKind.Apocenter
            else:
                raise InternalError()   # at least one of the planet distance slopes should have been nonzero
            search = Search(_planet_distance_slope, context, t1, t2, 1.0)
            if search is None:
                raise InternalError()   # failed to find where planet distance slope passed through zero

            dist = HelioDistance(body, search)
            return Apsis(search, kind, dist)
        t1 = t2
        m1 = m2
        iter_count += 1
    raise InternalError()   # should have found planet apsis within 2 planet orbits


def NextPlanetApsis(body, apsis):
    """Finds the next planetary perihelion or aphelion event in a series.

    This function requires an #Apsis value obtained from a call
    to #SearchPlanetApsis or `NextPlanetApsis`.
    Given an aphelion event, this function finds the next perihelion event, and vice versa.
    See #SearchPlanetApsis for more details.

    Parameters
    ----------
    body : Body
        The planet for which to find the next perihelion/aphelion event.
        Not allowed to be `Body.Sun` or `Body.Moon`.
        Must match the body passed into the call that produced the `apsis` parameter.
    apsis : Apsis
        An apsis event obtained from a call to #SearchPlanetApsis or `NextPlanetApsis`.

    Returns
    -------
    Apsis
    """
    if apsis.kind not in [ApsisKind.Apocenter, ApsisKind.Pericenter]:
        raise Error('Parameter "apsis" contains an invalid "kind" value.')
    # Skip 1/4 of an orbit before starting search again.
    skip = 0.25 * _PlanetOrbitalPeriod[body.value]
    time = apsis.time.AddDays(skip)
    next_apsis = SearchPlanetApsis(body, time)
    # Verify that we found the opposite apsis from the previous one.
    if next_apsis.kind.value + apsis.kind.value != 1:
        raise InternalError()   # should have found opposite planetary apsis type
    return next_apsis


def _PlanetExtreme(body, kind, start_time, dayspan):
    direction = +1.0 if (kind == ApsisKind.Apocenter) else -1.0
    npoints = 10
    while True:
        interval = dayspan / (npoints - 1)
        # Iterate until uncertainty is less than one minute.
        if interval < 1/1440:
            apsis_time = start_time.AddDays(interval/2)
            dist_au = HelioDistance(body, apsis_time)
            return Apsis(apsis_time, kind, dist_au)
        for i in range(npoints):
            time = start_time.AddDays(i * interval)
            dist = direction * HelioDistance(body, time)
            if i==0 or dist > best_dist:
                best_i = i
                best_dist = dist
        # Narrow in on the extreme point.
        start_time = start_time.AddDays((best_i - 1) * interval)
        dayspan = 2 * interval


def _BruteSearchPlanetApsis(body, startTime):
    # Neptune is a special case for two reasons:
    # 1. Its orbit is nearly circular (low orbital eccentricity).
    # 2. It is so distant from the Sun that the orbital period is very long.
    # Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
    # to be so significant that there are 3 local minima in the distance-vs-time curve
    # near each apsis. Therefore, unlike for other planets, we can't use an optimized
    # algorithm for finding dr/dt = 0.
    # Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
    # heliocentric distance.
    #
    # There is a similar problem in the TOP2013 model for Pluto:
    # Its position vector has high-frequency oscillations that confuse the
    # slope-based determination of apsides.
    #
    # Rewind approximately 30 degrees in the orbit,
    # then search forward for 270 degrees.
    # This is a very cautious way to prevent missing an apsis.
    # Typically we will find two apsides, and we pick whichever
    # apsis is ealier, but after startTime.
    # Sample points around this orbital arc and find when the distance
    # is greatest and smallest.
    period = _PlanetOrbitalPeriod[body.value]
    t1 = startTime.AddDays(period * ( -30 / 360))
    t2 = startTime.AddDays(period * (+270 / 360))
    t_min = t_max = t1
    npoints = 100
    interval = (t2.ut - t1.ut) / (npoints - 1)

    for i in range(npoints):
        time = t1.AddDays(i * interval)
        dist = HelioDistance(body, time)
        if i == 0:
            min_dist = max_dist = dist
        else:
            if dist > max_dist:
                max_dist = dist
                t_max = time
            if dist < min_dist:
                min_dist = dist
                t_min = time

    perihelion = _PlanetExtreme(body, ApsisKind.Pericenter, t_min.AddDays(-2*interval), 4*interval)
    aphelion = _PlanetExtreme(body, ApsisKind.Apocenter, t_max.AddDays(-2*interval), 4*interval)
    if perihelion.time.tt >= startTime.tt:
        if startTime.tt <= aphelion.time.tt < perihelion.time.tt:
            return aphelion
        return perihelion
    if aphelion.time.tt >= startTime.tt:
        return aphelion
    raise InternalError()   # failed to find Neptune apsis


def VectorFromSphere(sphere, time):
    """Converts spherical coordinates to Cartesian coordinates.

    Given spherical coordinates and a time at which they are valid,
    returns a vector of Cartesian coordinates. The returned value
    includes the time, as required by all `Time` objects.

    Parameters
    ----------
    sphere : Spherical
        Spherical coordinates to be converted.
    time : Time
        The time that should be included in the returned vector.

    Returns
    -------
    Vector
        The vector form of the supplied spherical coordinates.
    """
    radlat = math.radians(sphere.lat)
    radlon = math.radians(sphere.lon)
    rcoslat = sphere.dist * math.cos(radlat)
    return Vector(
        rcoslat * math.cos(radlon),
        rcoslat * math.sin(radlon),
        sphere.dist * math.sin(radlat),
        time
    )


def EquatorFromVector(vec):
    """Given an equatorial vector, calculates equatorial angular coordinates.

    Parameters
    ----------
    vec : Vector
        A vector in an equatorial coordinate system.

    Returns
    -------
    Equatorial
        Angular coordinates expressed in the same equatorial system as `vec`.
    """
    sphere = SphereFromVector(vec)
    return Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist, vec)


def SphereFromVector(vector):
    """Converts Cartesian coordinates to spherical coordinates.

    Given a Cartesian vector, returns latitude, longitude, and distance.

    Parameters
    ----------
    vector : Vector
        Cartesian vector to be converted to spherical coordinates.

    Returns
    -------
    Spherical
        Spherical coordinates that are equivalent to the given vector.
    """
    xyproj = vector.x*vector.x + vector.y*vector.y
    dist = math.sqrt(xyproj + vector.z*vector.z)
    if xyproj == 0.0:
        if vector.z == 0.0:
            raise Exception('Zero-length vector not allowed.')
        lon = 0.0
        if vector.z < 0.0:
            lat = -90.0
        else:
            lat = +90.0
    else:
        lon = math.degrees(math.atan2(vector.y, vector.x))
        if lon < 0.0:
            lon += 360.0
        lat = math.degrees(math.atan2(vector.z, math.sqrt(xyproj)))
    return Spherical(lat, lon, dist)


def _ToggleAzimuthDirection(az):
    az = 360.0 - az
    if az >= 360.0:
        az -= 360.0
    elif az < 0.0:
        az += 360.0
    return az


def VectorFromHorizon(sphere, time, refraction):
    """Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.

    Parameters
    ----------
    sphere : Spherical
        A structure that contains apparent horizontal coordinates:
        `lat` holds the refracted azimuth angle,
        `lon` holds the azimuth in degrees clockwise from north,
        and `dist` holds the distance from the observer to the object in AU.
    time : Time
        The date and time of the observation. This is needed because the returned
        vector object requires a valid time value when passed to certain other functions.
    refraction : Refraction
        See remarks in function #RefractionAngle.

    Returns
    -------
    Vector
        A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
    """
    lon = _ToggleAzimuthDirection(sphere.lon)
    lat = sphere.lat + InverseRefractionAngle(refraction, sphere.lat)
    xsphere = Spherical(lat, lon, sphere.dist)
    return VectorFromSphere(xsphere, time)


def HorizonFromVector(vector, refraction):
    """Converts Cartesian coordinates to horizontal coordinates.

    Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.

    *IMPORTANT:* This function differs from `SphereFromVector` in two ways:
    - `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
      from north (e.g., west = +90), but this function represents a clockwise rotation
      (e.g., east = +90). The difference is because `SphereFromVector` is intended
      to preserve the vector "right-hand rule", while this function defines azimuth in a more
      traditional way as used in navigation and cartography.
    - This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.

    The returned object contains the azimuth in `lon`.
    It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

    The altitude is stored in `lat`.

    The distance to the observed object is stored in `dist`,
    and is expressed in astronomical units (AU).

    Parameters
    ----------
    vector : Vector
        Cartesian vector to be converted to horizontal angular coordinates.
    refraction : Refraction
        See comments in the #RefractionAngle function.
    """
    sphere = SphereFromVector(vector)
    sphere.lon = _ToggleAzimuthDirection(sphere.lon)
    sphere.lat += RefractionAngle(refraction, sphere.lat)
    return sphere


def InverseRotation(rotation):
    """Calculates the inverse of a rotation matrix.

    Given a rotation matrix that performs some coordinate transform,
    this function returns the matrix that reverses that trasnform.

    Parameters
    ----------
    rotation : RotationMatrix
        The rotation matrix to be inverted.

    Returns
    -------
    RotationMatrix
        The inverse rotation matrix.
    """
    return RotationMatrix([
        [rotation.rot[0][0], rotation.rot[1][0], rotation.rot[2][0]],
        [rotation.rot[0][1], rotation.rot[1][1], rotation.rot[2][1]],
        [rotation.rot[0][2], rotation.rot[1][2], rotation.rot[2][2]]
    ])


def CombineRotation(a, b):
    """Creates a rotation based on applying one rotation followed by another.

    Given two rotation matrices, returns a combined rotation matrix that is
    equivalent to rotating based on the first matrix, followed by the second.

    Parameters
    ----------
    a : RotationMatrix
        The first rotation to apply.

    b : RotationMatrix
        The second rotation to apply.

    Returns
    -------
    RotationMatrix
        The combined rotation matrix.
    """
    # Use matrix multiplication: c = b*a.
    # We put 'b' on the left and 'a' on the right because,
    # just like when you use a matrix M to rotate a vector V,
    # you put the M on the left in the product M*V.
    # We can think of this as 'b' rotating all the 3 column vectors in 'a'.

    return RotationMatrix([
        [
            b.rot[0][0]*a.rot[0][0] + b.rot[1][0]*a.rot[0][1] + b.rot[2][0]*a.rot[0][2],
            b.rot[0][1]*a.rot[0][0] + b.rot[1][1]*a.rot[0][1] + b.rot[2][1]*a.rot[0][2],
            b.rot[0][2]*a.rot[0][0] + b.rot[1][2]*a.rot[0][1] + b.rot[2][2]*a.rot[0][2]
        ],
        [
            b.rot[0][0]*a.rot[1][0] + b.rot[1][0]*a.rot[1][1] + b.rot[2][0]*a.rot[1][2],
            b.rot[0][1]*a.rot[1][0] + b.rot[1][1]*a.rot[1][1] + b.rot[2][1]*a.rot[1][2],
            b.rot[0][2]*a.rot[1][0] + b.rot[1][2]*a.rot[1][1] + b.rot[2][2]*a.rot[1][2]
        ],
        [
            b.rot[0][0]*a.rot[2][0] + b.rot[1][0]*a.rot[2][1] + b.rot[2][0]*a.rot[2][2],
            b.rot[0][1]*a.rot[2][0] + b.rot[1][1]*a.rot[2][1] + b.rot[2][1]*a.rot[2][2],
            b.rot[0][2]*a.rot[2][0] + b.rot[1][2]*a.rot[2][1] + b.rot[2][2]*a.rot[2][2]
        ]
    ])


def IdentityMatrix():
    """Creates an identity rotation matrix.

    Returns a rotation matrix that has no effect on orientation.
    This matrix can be the starting point for other operations,
    such as using a series of calls to #Pivot to
    create a custom rotation matrix.

    Returns
    -------
    RotationMatrix
        The identity rotation matrix.
    """
    return RotationMatrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])


def Pivot(rotation, axis, angle):
    """Re-orients a rotation matrix by pivoting it by an angle around one of its axes.

    Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
    this function pivots the rotation matrix by that angle around that coordinate axis.

    For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
    to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
    of a telescope camera pointed at a given body, you can use `Pivot` twice:
    (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
    western axis by the body's altitude angle. The resulting rotation matrix will then
    reorient ECL coordinates to the orientation of your telescope camera.

    Parameters
    ----------
    rotation : RotationMatrix
        The input rotation matrix.
    axis : int
        An integer that selects which coordinate axis to rotate around:
        0 = x, 1 = y, 2 = z. Any other value will cause an exception.
    angle : float
        An angle in degrees indicating the amount of rotation around the specified axis.
        Positive angles indicate rotation counterclockwise as seen from the positive
        direction along that axis, looking towards the origin point of the orientation system.
        Any finite number of degrees is allowed, but best precision will result from
        keeping `angle` in the range [-360, +360].

    Returns
    -------
    RotationMatrix
        A pivoted matrix object.
    """
    # Check for an invalid coordinate axis.
    if axis not in [0, 1, 2]:
        raise Error('Invalid axis {}. Must be [0, 1, 2].'.format(axis))

    radians = math.radians(angle)
    c = math.cos(radians)
    s = math.sin(radians)

    # We need to maintain the "right-hand" rule, no matter which
    # axis was selected. That means we pick (i, j, k) axis order
    # such that the following vector cross product is satisfied:
    # i x j = k
    i = (axis + 1) % 3
    j = (axis + 2) % 3
    k = axis

    rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    rot[i][i] = c*rotation.rot[i][i] - s*rotation.rot[i][j]
    rot[i][j] = s*rotation.rot[i][i] + c*rotation.rot[i][j]
    rot[i][k] = rotation.rot[i][k]

    rot[j][i] = c*rotation.rot[j][i] - s*rotation.rot[j][j]
    rot[j][j] = s*rotation.rot[j][i] + c*rotation.rot[j][j]
    rot[j][k] = rotation.rot[j][k]

    rot[k][i] = c*rotation.rot[k][i] - s*rotation.rot[k][j]
    rot[k][j] = s*rotation.rot[k][i] + c*rotation.rot[k][j]
    rot[k][k] = rotation.rot[k][k]

    return RotationMatrix(rot)


def RotateVector(rotation, vector):
    """Applies a rotation to a vector, yielding a rotated vector.

    This function transforms a vector in one orientation to a vector
    in another orientation.

    Parameters
    ----------
    rotation : RotationMatrix
        A rotation matrix that specifies how the orientation of the vector is to be changed.
    vector : Vector
        The vector whose orientation is to be changed.

    Returns
    -------
    Vector
        A vector in the orientation specified by `rotation`.
    """
    return Vector(
        rotation.rot[0][0]*vector.x + rotation.rot[1][0]*vector.y + rotation.rot[2][0]*vector.z,
        rotation.rot[0][1]*vector.x + rotation.rot[1][1]*vector.y + rotation.rot[2][1]*vector.z,
        rotation.rot[0][2]*vector.x + rotation.rot[1][2]*vector.y + rotation.rot[2][2]*vector.z,
        vector.t
    )


def RotateState(rotation, state):
    """Applies a rotation to a state vector, yielding a rotated state vector.

    This function transforms a state vector in one orientation to a
    state vector in another orientation. Both the position and velocity
    vectors are rotated the same way.

    Parameters
    ----------
    rotation : RotationMatrix
        A rotation matrix that specifies how the orientation of the vector is to be changed.
    state : StateVector
        The state vector whose orientation is to be changed.

    Returns
    -------
    StateVector
        A state vector in the orientation specified by `rotation`.
    """
    return StateVector(
        rotation.rot[0][0]*state.x + rotation.rot[1][0]*state.y + rotation.rot[2][0]*state.z,
        rotation.rot[0][1]*state.x + rotation.rot[1][1]*state.y + rotation.rot[2][1]*state.z,
        rotation.rot[0][2]*state.x + rotation.rot[1][2]*state.y + rotation.rot[2][2]*state.z,
        rotation.rot[0][0]*state.vx + rotation.rot[1][0]*state.vy + rotation.rot[2][0]*state.vz,
        rotation.rot[0][1]*state.vx + rotation.rot[1][1]*state.vy + rotation.rot[2][1]*state.vz,
        rotation.rot[0][2]*state.vx + rotation.rot[1][2]*state.vy + rotation.rot[2][2]*state.vz,
        state.t
    )


def Rotation_EQJ_ECL():
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using equator at J2000 epoch.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to ECL.
    """
    # ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    c = 0.9174821430670688    # cos(ob)
    s = 0.3977769691083922    # sin(ob)
    return RotationMatrix([
        [ 1,  0,  0],
        [ 0, +c, -s],
        [ 0, +s, +c]
    ])


def Rotation_ECL_EQJ():
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: EQJ = equatorial system, using equator at J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to EQJ.
    """
    # ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    c = 0.9174821430670688    # cos(ob)
    s = 0.3977769691083922    # sin(ob)
    return RotationMatrix([
        [ 1,  0,  0],
        [ 0, +c, +s],
        [ 0, -s, +c]
    ])

def Rotation_EQJ_EQD(time):
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using equator at J2000 epoch.
    Target: EQD = equatorial system, using equator of the specified date/time.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator defines the target orientation.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to EQD at `time`.
    """
    prec = _precession_rot(time, _PrecessDir.From2000)
    nut = _nutation_rot(time, _PrecessDir.From2000)
    return CombineRotation(prec, nut)


def Rotation_EQD_EQJ(time):
    """Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of the specified date/time.
    Target: EQJ = equatorial system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator defines the source orientation.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD at `time` to EQJ.
    """
    nut = _nutation_rot(time, _PrecessDir.Into2000)
    prec = _precession_rot(time, _PrecessDir.Into2000)
    return CombineRotation(nut, prec)


def Rotation_EQD_HOR(time, observer):
    """Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of the specified date/time.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value
    to a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator applies.
    observer: Observer
        A location near the Earth's mean sea level that defines the observer's location.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    sinlat = math.sin(math.radians(observer.latitude))
    coslat = math.cos(math.radians(observer.latitude))
    sinlon = math.sin(math.radians(observer.longitude))
    coslon = math.cos(math.radians(observer.longitude))
    uze = [coslat * coslon, coslat * sinlon, sinlat]
    une = [-sinlat * coslon, -sinlat * sinlon, coslat]
    uwe = [sinlon, -coslon, 0.0]
    spin_angle = -15.0 * _sidereal_time(time)
    uz = _spin(spin_angle, uze)
    un = _spin(spin_angle, une)
    uw = _spin(spin_angle, uwe)
    return RotationMatrix([
        [un[0], uw[0], uz[0]],
        [un[1], uw[1], uz[1]],
        [un[2], uw[2], uz[2]],
    ])


def Rotation_HOR_EQD(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system (x=North, y=West, z=Zenith).
    Target: EQD = equatorial system, using equator of the specified date/time.

    Parameters
    ----------
    time : Time
        The date and time at which the Earth's equator applies.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to EQD at `time` and for `observer`.
    """
    rot = Rotation_EQD_HOR(time, observer)
    return InverseRotation(rot)


def Rotation_HOR_EQJ(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system (x=North, y=West, z=Zenith).
    Target: EQJ = equatorial system, using equator at the J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the observation.
    observer : Observer
        A location near the Earth's mean sea level that define's the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to EQD at `time` and for `observer`.
    """
    hor_eqd = Rotation_HOR_EQD(time, observer)
    eqd_eqj = Rotation_EQD_EQJ(time)
    return CombineRotation(hor_eqd, eqd_eqj)


def Rotation_EQJ_HOR(time, observer):
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using the equator at the J2000 epoch.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value to
    a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time of the desired horizontal orientation.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    rot = Rotation_HOR_EQJ(time, observer)
    return InverseRotation(rot)


def Rotation_EQD_ECL(time):
    """Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQD = equatorial system, using equator of date.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the source equator.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQD to ECL.
    """
    eqd_eqj = Rotation_EQD_EQJ(time)
    eqj_ecl = Rotation_EQJ_ECL()
    return CombineRotation(eqd_eqj, eqj_ecl)


def Rotation_ECL_EQD(time):
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: EQD = equatorial system, using equator of date.

    Parameters
    ----------
    time : Time
        The date and time of the desired equator.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to EQD.
    """
    rot = Rotation_EQD_ECL(time)
    return InverseRotation(rot)


def Rotation_ECL_HOR(time, observer):
    """Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: ECL = ecliptic system, using equator at J2000 epoch.
    Target: HOR = horizontal system.

    Use #HorizonFromVector to convert the return value
    to a traditional altitude/azimuth pair.

    Parameters
    ----------
    time : Time
        The date and time of the desired horizontal orientation.
    observer : Observer
        A location near the Earth's mean sea level that defines the observer's horizon.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts ECL to HOR at `time` and for `observer`.
        The components of the horizontal vector are:
        x = north, y = west, z = zenith (straight up from the observer).
        These components are chosen so that the "right-hand rule" works for the vector
        and so that north represents the direction where azimuth = 0.
    """
    ecl_eqd = Rotation_ECL_EQD(time)
    eqd_hor = Rotation_EQD_HOR(time, observer)
    return CombineRotation(ecl_eqd, eqd_hor)


def Rotation_HOR_ECL(time, observer):
    """Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: HOR = horizontal system.
    Target: ECL = ecliptic system, using equator at J2000 epoch.

    Parameters
    ----------
    time : Time
        The date and time of the horizontal observation.
    observer : Observer
        The location of the horizontal observer.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts HOR to ECL.
    """
    rot = Rotation_ECL_HOR(time, observer)
    return InverseRotation(rot)

def Rotation_EQJ_GAL():
    """Calculates a rotation matrix from equatorial J2000 (EQJ) to galactic (GAL).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: EQJ = equatorial system, using the equator at the J2000 epoch.
    Target: GAL = galactic system (IAU 1958 definition).

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts EQJ to GAL.
    """
    # This rotation matrix was calculated by the following script
    # in this same source code repository:
    # demo/python/galeqj_matrix.py
    return RotationMatrix([
        [-0.0548624779711344, +0.4941095946388765, -0.8676668813529025],
        [-0.8734572784246782, -0.4447938112296831, -0.1980677870294097],
        [-0.4838000529948520, +0.7470034631630423, +0.4559861124470794]
    ])

def Rotation_GAL_EQJ():
    """Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).

    This is one of the family of functions that returns a rotation matrix
    for converting from one orientation to another.
    Source: GAL = galactic system (IAU 1958 definition).
    Target: EQJ = equatorial system, using the equator at the J2000 epoch.

    Returns
    -------
    RotationMatrix
        A rotation matrix that converts GAL to EQJ.
    """
    # This rotation matrix was calculated by the following script
    # in this same source code repository:
    # demo/python/galeqj_matrix.py
    return RotationMatrix([
        [-0.0548624779711344, -0.8734572784246782, -0.4838000529948520],
        [+0.4941095946388765, -0.4447938112296831, +0.7470034631630423],
        [-0.8676668813529025, -0.1980677870294097, +0.4559861124470794]
    ])

class ConstellationInfo:
    """Reports the constellation that a given celestial point lies within.

    The #Constellation function returns a `ConstellationInfo` object
    to report which constellation corresponds with a given point in the sky.
    Constellations are defined with respect to the B1875 equatorial system
    per IAU standard. Although the `Constellation` function requires J2000 equatorial
    coordinates as input, the returned object contains converted B1875 coordinates for reference.

    Attributes
    ----------
    symbol : string
        3-character mnemonic symbol for the constellation, e.g. "Ori".
    name : string
        Full name of constellation, e.g. "Orion".
    ra1875 : float
        Right ascension expressed in B1875 coordinates.
    dec1875 : float
        Declination expressed in B1875 coordinates.
    """
    def __init__(self, symbol, name, ra1875, dec1875):
        self.symbol = symbol
        self.name = name
        self.ra1875 = ra1875
        self.dec1875 = dec1875


_ConstelRot = None
_Epoch2000 = None
$ASTRO_CONSTEL()


def Constellation(ra, dec):
    """Determines the constellation that contains the given point in the sky.

    Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
    constellation that contains that point.

    Parameters
    ----------
    ra : float
        The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
    dec : float
        The declination (DEC) of a point in the sky, using the J2000 equatorial system.

    Returns
    -------
    ConstellationInfo
        A structure that contains the 3-letter abbreviation and full name
        of the constellation that contains the given (ra,dec), along with
        the converted B1875 (ra,dec) for that point.
    """
    global _ConstelRot, _Epoch2000

    if dec < -90.0 or dec > +90.0:
        raise Error('Invalid declination angle. Must be -90..+90.')

    # Clamp right ascension to [0, 24) sidereal hours.
    ra = math.fmod(ra, 24.0)
    if ra < 0.0:
        ra += 24.0

    # Lazy-initialize rotation matrix.
    if _ConstelRot is None:
        # Need to calculate the B1875 epoch. Based on this:
        # https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years
        # B = 1900 + (JD - 2415020.31352) / 365.242198781
        # I'm interested in using TT instead of JD, giving:
        # B = 1900 + ((TT+2451545) - 2415020.31352) / 365.242198781
        # B = 1900 + (TT + 36524.68648) / 365.242198781
        # TT = 365.242198781*(B - 1900) - 36524.68648 = -45655.741449525
        # But the Time constructor wants UT, not TT.
        # Near that date, I get a historical correction of ut-tt = 3.2 seconds.
        # That gives UT = -45655.74141261017 for the B1875 epoch,
        # or 1874-12-31T18:12:21.950Z.
        _ConstelRot = Rotation_EQJ_EQD(Time(-45655.74141261017))
        _Epoch2000 = Time(0.0)

    # Convert coordinates from J2000 to B1875.
    sph2000 = Spherical(dec, 15.0 * ra, 1.0)
    vec2000 = VectorFromSphere(sph2000, _Epoch2000)
    vec1875 = RotateVector(_ConstelRot, vec2000)
    equ1875 = EquatorFromVector(vec1875)

    # Convert DEC from degrees and RA from hours, into compact angle units used in the _ConstelBounds table.
    x_dec = 24.0 * equ1875.dec
    x_ra = (24.0 * 15.0) * equ1875.ra

    # Search for the constellation using the B1875 coordinates.
    for b in _ConstelBounds:
        index, ra_lo, ra_hi, b_dec = b
        if (b_dec <= x_dec) and (ra_lo <= x_ra < ra_hi):
            symbol, name = _ConstelNames[index]
            return ConstellationInfo(symbol, name, equ1875.ra, equ1875.dec)

    # This should never happen!
    raise Error('Unable to find constellation for given coordinates.')


@enum.unique
class EclipseKind(enum.Enum):
    """Selects if/how to correct for atmospheric refraction.

    Some functions allow enabling or disabling atmospheric refraction
    for the calculated apparent position of a celestial body
    as seen by an observer on the surface of the Earth.

    Values
    ------
    Invalid: No eclipse found.
    Penumbral: A penumbral lunar eclipse. (Never used for a solar eclipse.)
    Partial: A partial lunar/solar eclipse.
    Annular: An annular solar eclipse. (Never used for a lunar eclipse.)
    Total: A total lunar/solar eclipse.
    """
    Invalid = 0
    Penumbral = 1
    Partial = 2
    Annular = 3
    Total = 4


class _ShadowInfo:
    def __init__(self, time, u, r, k, p, target, direction):
        self.time = time
        self.u = u   # dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is
        self.r = r   # km distance between center of Moon and the line passing through the centers of the Sun and Earth.
        self.k = k   # umbra radius in km, at the shadow plane
        self.p = p   # penumbra radius in km, at the shadow plane
        self.target = target        # vector from center of shadow-casting body to target location that might receive the shadow
        self.dir = direction        # vector from center of Sun to center of shadow-casting body


def _CalcShadow(body_radius_km, time, target, sdir):
    u = (sdir.x*target.x + sdir.y*target.y + sdir.z*target.z) / (sdir.x*sdir.x + sdir.y*sdir.y + sdir.z*sdir.z)
    dx = (u * sdir.x) - target.x
    dy = (u * sdir.y) - target.y
    dz = (u * sdir.z) - target.z
    r = KM_PER_AU * math.sqrt(dx*dx + dy*dy + dz*dz)
    k = +_SUN_RADIUS_KM - (1.0 + u)*(_SUN_RADIUS_KM - body_radius_km)
    p = -_SUN_RADIUS_KM + (1.0 + u)*(_SUN_RADIUS_KM + body_radius_km)
    return _ShadowInfo(time, u, r, k, p, target, sdir)


def _EarthShadow(time):
    e = _CalcEarth(time)
    m = GeoMoon(time)
    return _CalcShadow(_EARTH_ECLIPSE_RADIUS_KM, time, m, e)


def _MoonShadow(time):
    # This is a variation on the logic in _EarthShadow().
    # Instead of a heliocentric Earth and a geocentric Moon,
    # we want a heliocentric Moon and a lunacentric Earth.
    h = _CalcEarth(time)    # heliocentric Earth
    m = GeoMoon(time)       # geocentric Moon
    # Calculate lunacentric Earth.
    e = Vector(-m.x, -m.y, -m.z, m.t)
    # Convert geocentric moon to heliocentric Moon.
    m.x += h.x
    m.y += h.y
    m.z += h.z
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, e, m)


def _LocalMoonShadow(time, observer):
    # Calculate observer's geocentric position.
    # For efficiency, do this first, to populate the earth rotation parameters in 'time'.
    # That way they can be recycled instead of recalculated.
    pos = _geo_pos(time, observer)
    h = _CalcEarth(time)     # heliocentric Earth
    m = GeoMoon(time)        # geocentric Moon
    # Calculate lunacentric location of an observer on the Earth's surface.
    o = Vector(pos[0] - m.x, pos[1] - m.y, pos[2] - m.z, time)
    # Convert geocentric moon to heliocentric Moon.
    m.x += h.x
    m.y += h.y
    m.z += h.z
    return _CalcShadow(_MOON_MEAN_RADIUS_KM, time, o, m)


def _PlanetShadow(body, planet_radius_km, time):
    # Calculate light-travel-corrected vector from Earth to planet.
    g = GeoVector(body, time, False)
    # Calculate light-travel-corrected vector from Earth to Sun.
    e = GeoVector(Body.Sun, time, False)
    # Deduce light-travel-corrected vector from Sun to planet.
    p = Vector(g.x - e.x, g.y - e.y, g.z - e.z, time)
    # Calcluate Earth's position from the planet's point of view.
    e.x = -g.x
    e.y = -g.y
    e.z = -g.z
    return _CalcShadow(planet_radius_km, time, e, p)


def _ShadowDistanceSlope(shadowfunc, time):
    dt = 1.0 / 86400.0
    t1 = time.AddDays(-dt)
    t2 = time.AddDays(+dt)
    shadow1 = shadowfunc(t1)
    shadow2 = shadowfunc(t2)
    return (shadow2.r - shadow1.r) / dt


def _PlanetShadowSlope(context, time):
    (body, planet_radius_km) = context
    dt = 1.0 / 86400.0
    shadow1 = _PlanetShadow(body, planet_radius_km, time.AddDays(-dt))
    shadow2 = _PlanetShadow(body, planet_radius_km, time.AddDays(+dt))
    return (shadow2.r - shadow1.r) / dt


def _PeakEarthShadow(search_center_time):
    window = 0.03        # initial search window, in days, before/after given time
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, _EarthShadow, t1, t2, 1.0)
    return _EarthShadow(tx)


def _PeakMoonShadow(search_center_time):
    window = 0.03        # initial search window, in days, before/after given time
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, _MoonShadow, t1, t2, 1.0)
    return _MoonShadow(tx)

def _PeakLocalMoonShadow(search_center_time, observer):
    # Search for the time near search_center_time that the Moon's shadow comes
    # closest to the given observer.
    window = 0.2
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_ShadowDistanceSlope, lambda time: _LocalMoonShadow(time, observer), t1, t2, 1.0)
    return _LocalMoonShadow(tx, observer)

def _PeakPlanetShadow(body, planet_radius_km, search_center_time):
    # Search for when the body's shadow is closest to the center of the Earth.
    window = 1.0     # days before/after inferior conjunction to search for minimum shadow distance.
    t1 = search_center_time.AddDays(-window)
    t2 = search_center_time.AddDays(+window)
    tx = Search(_PlanetShadowSlope, (body, planet_radius_km), t1, t2, 1.0)
    return _PlanetShadow(body, planet_radius_km, tx)

class _ShadowDiffContext:
    def __init__(self, radius_limit, direction):
        self.radius_limit = radius_limit
        self.direction = direction

def _ShadowDiff(context, time):
    return context.direction * (_EarthShadow(time).r - context.radius_limit)


class LunarEclipseInfo:
    """Returns information about a lunar eclipse.

    Returned by #SearchLunarEclipse or #NextLunarEclipse
    to report information about a lunar eclipse event.
    When a lunar eclipse is found, it is classified as penumbral, partial, or total.
    Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
    by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
    Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
    Total eclipses occur when the entire Moon passes into the Earth's umbra.

    The `kind` field thus holds one of the values `EclipseKind.Penumbral`, `EclipseKind.Partial`,
    or `EclipseKind.Total`, depending on the kind of lunar eclipse found.

    Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.

    Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
    of the eclipse, which is half of the amount of time the eclipse spends in each
    phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
    By converting from minutes to days, and subtracting/adding with `peak`, the caller
    may determine the date and time of the beginning/end of each eclipse phase.

    Attributes
    ----------
    kind : string
         The type of lunar eclipse found.
    peak : Time
         The time of the eclipse at its peak.
    sd_penum : float
         The semi-duration of the penumbral phase in minutes.
    sd_partial : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    sd_total : float
         The semi-duration of the penumbral phase in minutes, or 0.0 if none.
    """
    def __init__(self, kind, peak, sd_penum, sd_partial, sd_total):
        self.kind = kind
        self.peak = peak
        self.sd_penum = sd_penum
        self.sd_partial = sd_partial
        self.sd_total = sd_total


class GlobalSolarEclipseInfo:
    """Reports the time and geographic location of the peak of a solar eclipse.

    Returned by #SearchGlobalSolarEclipse or #NextGlobalSolarEclipse
    to report information about a solar eclipse event.

    Field `peak` holds the date and time of the peak of the eclipse, defined as
    the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.

    The eclipse is classified as partial, annular, or total, depending on the
    maximum amount of the Sun's disc obscured, as seen at the peak location
    on the surface of the Earth.

    The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
    An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
    to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
    A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
    observes either a total or annular eclipse.

    If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
    fields give the geographic coordinates of the center of the Moon's shadow projected
    onto the daytime side of the Earth at the instant of the eclipse's peak.
    If `kind` has any other value, `latitude` and `longitude` are undefined and should
    not be used.

    Attributes
    ----------
    kind : EclipseKind
        The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    peak : Time
        The date and time of the eclipse at its peak.
    distance : float
        The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.
    latitude : float
        The geographic latitude at the center of the peak eclipse shadow.
    longitude : float
        The geographic longitude at the center of the peak eclipse shadow.
    """
    def __init__(self, kind, peak, distance, latitude, longitude):
        self.kind = kind
        self.peak = peak
        self.distance = distance
        self.latitude = latitude
        self.longitude = longitude


class EclipseEvent:
    """Holds a time and the observed altitude of the Sun at that time.

    When reporting a solar eclipse observed at a specific location on the Earth
    (a "local" solar eclipse), a series of events occur. In addition
    to the time of each event, it is important to know the altitude of the Sun,
    because each event may be invisible to the observer if the Sun is below
    the horizon (i.e. it at night).

    If `altitude` is negative, the event is theoretical only; it would be
    visible if the Earth were transparent, but the observer cannot actually see it.
    If `altitude` is positive but less than a few degrees, visibility will be impaired by
    atmospheric interference (sunrise or sunset conditions).

    Attributes
    ----------
    time : Time
        The date and time of the event.
    altitude : float
        The angular altitude of the center of the Sun above/below the horizon, at `time`,
        corrected for atmospheric refraction and expressed in degrees.
    """
    def __init__(self, time, altitude):
        self.time = time
        self.altitude = altitude


class LocalSolarEclipseInfo:
    """Information about a solar eclipse as seen by an observer at a given time and geographic location.

    Returned by #SearchLocalSolarEclipse or #NextLocalSolarEclipse
    to report information about a solar eclipse as seen at a given geographic location.

    When a solar eclipse is found, it is classified as partial, annular, or total.
    The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    A partial solar eclipse is when the Moon does not line up directly enough with the Sun
    to completely block the Sun's light from reaching the observer.
    An annular eclipse occurs when the Moon's disc is completely visible against the Sun
    but the Moon is too far away to completely block the Sun's light; this leaves the
    Sun with a ring-like appearance.
    A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
    Sun just right to completely block all sunlight from reaching the observer.

    There are 5 "event" fields, each of which contains a time and a solar altitude.
    Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
    The fields `partial_begin` and `partial_end` are always set, and indicate when
    the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
    `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
    When an event field is valid, the caller must also check its `altitude` field to
    see whether the Sun is above the horizon at the time indicated by the `time` field.
    See #EclipseEvent for more information.

    Attributes
    ----------
    kind : EclipseKind
        The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    partial_begin : EclipseEvent
        The time and Sun altitude at the beginning of the eclipse.
    total_begin : EclipseEvent
        If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise `None`.
    peak : EclipseEvent
        The time and Sun altitude when the eclipse reaches its peak.
    total_end : EclipseEvent
        If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise `None`.
    partial_end : EclipseEvent
        The time and Sun altitude at the end of the eclipse.
    """
    def __init__(self, kind, partial_begin, total_begin, peak, total_end, partial_end):
        self.kind = kind
        self.partial_begin = partial_begin
        self.total_begin = total_begin
        self.peak = peak
        self.total_end = total_end
        self.partial_end = partial_end


def _EclipseKindFromUmbra(k):
    # The umbra radius tells us what kind of eclipse the observer sees.
    # If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular.
    # HACK: I added a tiny bias (14 meters) to match Espenak test data.
    if k > 0.014:
        return EclipseKind.Total
    return EclipseKind.Annular

class _LocalTransitionContext:
    def __init__(self, observer, direction, func):
        self.observer = observer
        self.direction = direction
        self.func = func


def _LocalTransitionFunc(context, time):
    shadow = _LocalMoonShadow(time, context.observer)
    return context.direction * context.func(shadow)


def _LocalEclipseTransition(observer, direction, func, t1, t2):
    context = _LocalTransitionContext(observer, direction, func)
    search = Search(_LocalTransitionFunc, context, t1, t2, 1.0)
    if search is None:
        raise Error('Local eclipse transition search failed')
    return _CalcEvent(observer, search)


def _CalcEvent(observer, time):
    altitude = _SunAltitude(time, observer)
    return EclipseEvent(time, altitude)


def _SunAltitude(time, observer):
    equ = Equator(Body.Sun, time, observer, True, True)
    hor = Horizon(time, observer, equ.ra, equ.dec, Refraction.Normal)
    return hor.altitude


def _local_partial_distance(shadow):
    # Must take the absolute value of the umbra radius 'k'
    # because it can be negative for an annular eclipse.
    return shadow.p - shadow.r


def _local_total_distance(shadow):
    return abs(shadow.k) - shadow.r


def _LocalEclipse(shadow, observer):
    PARTIAL_WINDOW = 0.2
    TOTAL_WINDOW = 0.01
    peak = _CalcEvent(observer, shadow.time)
    t1 = shadow.time.AddDays(-PARTIAL_WINDOW)
    t2 = shadow.time.AddDays(+PARTIAL_WINDOW)
    partial_begin = _LocalEclipseTransition(observer, +1.0, _local_partial_distance, t1, shadow.time)
    partial_end   = _LocalEclipseTransition(observer, -1.0, _local_partial_distance, shadow.time, t2)
    if shadow.r < abs(shadow.k):    # take absolute value of 'k' to handle annular eclipses too.
        t1 = shadow.time.AddDays(-TOTAL_WINDOW)
        t2 = shadow.time.AddDays(+TOTAL_WINDOW)
        total_begin = _LocalEclipseTransition(observer, +1.0, _local_total_distance, t1, shadow.time)
        total_end   = _LocalEclipseTransition(observer, -1.0, _local_total_distance, shadow.time, t2)
        kind = _EclipseKindFromUmbra(shadow.k)
    else:
        total_begin = None
        total_end = None
        kind = EclipseKind.Partial
    return LocalSolarEclipseInfo(kind, partial_begin, total_begin, peak, total_end, partial_end)


def _GeoidIntersect(shadow):
    kind = EclipseKind.Partial
    peak = shadow.time
    distance = shadow.r
    latitude = longitude = math.nan

    # We want to calculate the intersection of the shadow axis with the Earth's geoid.
    # First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
    # coordinates that are perfectly aligned with the Earth's equator at this
    # moment in time.
    rot = Rotation_EQJ_EQD(shadow.time)
    v = RotateVector(rot, shadow.dir)        # shadow-axis vector in equator-of-date coordinates
    e = RotateVector(rot, shadow.target)     # lunacentric Earth in equator-of-date coordinates

    # Convert all distances from AU to km.
    # But dilate the z-coordinates so that the Earth becomes a perfect sphere.
    # Then find the intersection of the vector with the sphere.
    # See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
    v.x *= KM_PER_AU
    v.y *= KM_PER_AU
    v.z *= KM_PER_AU / _EARTH_FLATTENING
    e.x *= KM_PER_AU
    e.y *= KM_PER_AU
    e.z *= KM_PER_AU / _EARTH_FLATTENING

    # Solve the quadratic equation that finds whether and where
    # the shadow axis intersects with the Earth in the dilated coordinate system.
    R = _EARTH_EQUATORIAL_RADIUS_KM
    A = v.x*v.x + v.y*v.y + v.z*v.z
    B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z)
    C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R
    radic = B*B - 4*A*C

    if radic > 0.0:
        # Calculate the closer of the two intersection points.
        # This will be on the day side of the Earth.
        u = (-B - math.sqrt(radic)) / (2 * A)

        # Convert lunacentric dilated coordinates to geocentric coordinates.
        px = u*v.x - e.x
        py = u*v.y - e.y
        pz = (u*v.z - e.z) * _EARTH_FLATTENING

        # Convert cartesian coordinates into geodetic latitude/longitude.
        proj = math.sqrt(px*px + py*py) * (_EARTH_FLATTENING * _EARTH_FLATTENING)
        if proj == 0.0:
            if pz > 0.0:
                latitude = +90.0
            else:
                latitude = -90.0
        else:
            latitude = math.degrees(math.atan(pz / proj))

        # Adjust longitude for Earth's rotation at the given UT.
        gast = _sidereal_time(peak)
        longitude = math.fmod(math.degrees(math.atan2(py, px)) - (15*gast), 360.0)
        if longitude <= -180.0:
            longitude += 360.0
        elif longitude > +180.0:
            longitude -= 360.0

        # We want to determine whether the observer sees a total eclipse or an annular eclipse.
        # We need to perform a series of vector calculations...
        # Calculate the inverse rotation matrix, so we can convert EQD to EQJ.
        inv = InverseRotation(rot)

        # Put the EQD geocentric coordinates of the observer into the vector 'o'.
        # Also convert back from kilometers to astronomical units.
        o = Vector(px / KM_PER_AU, py / KM_PER_AU, pz / KM_PER_AU, shadow.time)

        # Rotate the observer's geocentric EQD back to the EQJ system.
        o = RotateVector(inv, o)

        # Convert geocentric vector to lunacentric vector.
        o.x += shadow.target.x
        o.y += shadow.target.y
        o.z += shadow.target.z

        # Recalculate the shadow using a vector from the Moon's center toward the observer.
        surface = _CalcShadow(_MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir)

        # If we did everything right, the shadow distance should be very close to zero.
        # That's because we already determined the observer 'o' is on the shadow axis!
        if surface.r > 1.0e-9 or surface.r < 0.0:
            raise Error('Unexpected shadow distance from geoid intersection = {}'.format(surface.r))

        kind = _EclipseKindFromUmbra(surface.k)

    return GlobalSolarEclipseInfo(kind, peak, distance, latitude, longitude)


def _ShadowSemiDurationMinutes(center_time, radius_limit, window_minutes):
    # Search backwards and forwards from the center time until shadow axis distance crosses radius limit.
    window = window_minutes / (24.0 * 60.0)
    before = center_time.AddDays(-window)
    after  = center_time.AddDays(+window)
    t1 = Search(_ShadowDiff, _ShadowDiffContext(radius_limit, -1.0), before, center_time, 1.0)
    t2 = Search(_ShadowDiff, _ShadowDiffContext(radius_limit, +1.0), center_time, after, 1.0)
    if (t1 is None) or (t2 is None):
        raise Error('Failed to find shadow semiduration')
    return (t2.ut - t1.ut) * ((24.0 * 60.0) / 2.0)   # convert days to minutes and average the semi-durations.


def _MoonEclipticLatitudeDegrees(time):
    moon = _CalcMoon(time)
    return math.degrees(moon.geo_eclip_lat)


def SearchLunarEclipse(startTime):
    """Searches for a lunar eclipse.

    This function finds the first lunar eclipse that occurs after `startTime`.
    A lunar eclipse may be penumbral, partial, or total.
    See #LunarEclipseInfo for more information.
    To find a series of lunar eclipses, call this function once,
    then keep calling #NextLunarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    Parameters
    ----------
    startTime : Time
         The date and time for starting the search for a lunar eclipse.

    Returns
    -------
    LunarEclipseInfo
    """
    PruneLatitude = 1.8   # full Moon's ecliptic latitude above which eclipse is impossible
    fmtime = startTime
    for _ in range(12):
        # Search for the next full moon. Any eclipse will be near it.
        fullmoon = SearchMoonPhase(180, fmtime, 40)
        if fullmoon is None:
            raise Error('Cannot find full moon.')

        # Pruning: if the full Moon's ecliptic latitude is too large,
        # a lunar eclipse is not possible. Avoid needless work searching for
        # the minimum moon distance.
        eclip_lat = _MoonEclipticLatitudeDegrees(fullmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the full moon for the time when the center of the Moon
            # is closest to the line passing through the centers of the Sun and Earth.
            shadow = _PeakEarthShadow(fullmoon)
            if shadow.r < shadow.p + _MOON_MEAN_RADIUS_KM:
                # This is at least a penumbral eclipse. We will return a result.
                kind = EclipseKind.Penumbral
                sd_total = 0.0
                sd_partial = 0.0
                sd_penum = _ShadowSemiDurationMinutes(shadow.time, shadow.p + _MOON_MEAN_RADIUS_KM, 200.0)

                if shadow.r < shadow.k + _MOON_MEAN_RADIUS_KM:
                    # This is at least a partial eclipse.
                    kind = EclipseKind.Partial
                    sd_partial = _ShadowSemiDurationMinutes(shadow.time, shadow.k + _MOON_MEAN_RADIUS_KM, sd_penum)

                    if shadow.r + _MOON_MEAN_RADIUS_KM < shadow.k:
                        # This is a total eclipse.
                        kind = EclipseKind.Total
                        sd_total = _ShadowSemiDurationMinutes(shadow.time, shadow.k - _MOON_MEAN_RADIUS_KM, sd_partial)

                return LunarEclipseInfo(kind, shadow.time, sd_penum, sd_partial, sd_total)

        # We didn't find an eclipse on this full moon, so search for the next one.
        fmtime = fullmoon.AddDays(10)

    # This should never happen because there are always at least 2 full moons per year.
    raise Error('Failed to find lunar eclipse within 12 full moons.')


def NextLunarEclipse(prevEclipseTime):
    """Searches for the next lunar eclipse in a series.

     After using #SearchLunarEclipse to find the first lunar eclipse
     in a series, you can call this function to find the next consecutive lunar eclipse.
     Pass in the `peak` value from the #LunarEclipseInfo returned by the
     previous call to `SearchLunarEclipse` or `NextLunarEclipse`
     to find the next lunar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a full moon. Lunar eclipse search will start at the next full moon.

    Returns
    -------
    LunarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchLunarEclipse(startTime)


def SearchGlobalSolarEclipse(startTime):
    """Searches for a solar eclipse visible anywhere on the Earth's surface.

    This function finds the first solar eclipse that occurs after `startTime`.
    A solar eclipse may be partial, annular, or total.
    See #GlobalSolarEclipseInfo for more information.
    To find a series of solar eclipses, call this function once,
    then keep calling #NextGlobalSolarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    Parameters
    ----------
    startTime : Time
        The date and time for starting the search for a solar eclipse.

    Returns
    -------
    GlobalSolarEclipseInfo
    """
    PruneLatitude = 1.8     # Moon's ecliptic latitude beyond which eclipse is impossible
    # Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
    nmtime = startTime
    for _ in range(12):
        # Search for the next new moon. Any eclipse will be near it.
        newmoon = SearchMoonPhase(0.0, nmtime, 40.0)
        if newmoon is None:
            raise Error('Cannot find new moon')

        # Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
        eclip_lat = _MoonEclipticLatitudeDegrees(newmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the new moon for the time when the center of the Earth
            # is closest to the line passing through the centers of the Sun and Moon.
            shadow = _PeakMoonShadow(newmoon)
            if shadow.r < shadow.p + _EARTH_MEAN_RADIUS_KM:
                # This is at least a partial solar eclipse visible somewhere on Earth.
                # Try to find an intersection between the shadow axis and the Earth's oblate geoid.
                return _GeoidIntersect(shadow)

        # We didn't find an eclipse on this new moon, so search for the next one.
        nmtime = newmoon.AddDays(10.0)

    # Safety valve to prevent infinite loop.
    # This should never happen, because at least 2 solar eclipses happen per year.
    raise Error('Failed to find solar eclipse within 12 full moons.')


def NextGlobalSolarEclipse(prevEclipseTime):
    """Searches for the next global solar eclipse in a series.

    After using #SearchGlobalSolarEclipse to find the first solar eclipse
    in a series, you can call this function to find the next consecutive solar eclipse.
    Pass in the `peak` value from the #GlobalSolarEclipseInfo returned by the
    previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
    to find the next solar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a new moon. Solar eclipse search will start at the next new moon.

    Returns
    -------
    GlobalSolarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchGlobalSolarEclipse(startTime)


def SearchLocalSolarEclipse(startTime, observer):
    """Searches for a solar eclipse visible at a specific location on the Earth's surface.
    This function finds the first solar eclipse that occurs after `startTime`.
    A solar eclipse may be partial, annular, or total.
    See #LocalSolarEclipseInfo for more information.

    To find a series of solar eclipses, call this function once,
    then keep calling #NextLocalSolarEclipse as many times as desired,
    passing in the `peak` value returned from the previous call.

    IMPORTANT: An eclipse reported by this function might be partly or
    completely invisible to the observer due to the time of day.
    See #LocalSolarEclipseInfo for more information about this topic.

    Parameters
    ----------
    startTime : Time
        The date and time for starting the search for a solar eclipse.
    observer : Observer
        The geographic location of the observer.

    Returns
    -------
    LocalSolarEclipseInfo
    """
    PruneLatitude = 1.8   # Moon's ecliptic latitude beyond which eclipse is impossible

    # Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
    nmtime = startTime
    while True:
        # Search for the next new moon. Any eclipse will be near it.
        newmoon = SearchMoonPhase(0.0, nmtime, 40.0)

        # Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
        eclip_lat = _MoonEclipticLatitudeDegrees(newmoon)
        if abs(eclip_lat) < PruneLatitude:
            # Search near the new moon for the time when the observer
            # is closest to the line passing through the centers of the Sun and Moon.
            shadow = _PeakLocalMoonShadow(newmoon, observer)
            if shadow.r < shadow.p:
                # This is at least a partial solar eclipse for the observer.
                eclipse = _LocalEclipse(shadow, observer)

                # Ignore any eclipse that happens completely at night.
                # More precisely, the center of the Sun must be above the horizon
                # at the beginning or the end of the eclipse, or we skip the event.
                if eclipse.partial_begin.altitude > 0.0 or eclipse.partial_end.altitude > 0.0:
                    return eclipse

        # We didn't find an eclipse on this new moon, so search for the next one.
        nmtime = newmoon.AddDays(10.0)


def NextLocalSolarEclipse(prevEclipseTime, observer):
    """Searches for the next local solar eclipse in a series.

    After using #SearchLocalSolarEclipse to find the first solar eclipse
    in a series, you can call this function to find the next consecutive solar eclipse.
    Pass in the `peak` value from the #LocalSolarEclipseInfo returned by the
    previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
    to find the next solar eclipse.

    Parameters
    ----------
    prevEclipseTime : Time
        A date and time near a new moon. Solar eclipse search will start at the next new moon.
    observer : Observer
        The geographic location of the observer.

    Returns
    -------
    LocalSolarEclipseInfo
    """
    startTime = prevEclipseTime.AddDays(10.0)
    return SearchLocalSolarEclipse(startTime, observer)


class TransitInfo:
    """Information about a transit of Mercury or Venus, as seen from the Earth.

    Returned by #SearchTransit or #NextTransit to report
    information about a transit of Mercury or Venus.
    A transit is when Mercury or Venus passes between the Sun and Earth so that
    the other planet is seen in silhouette against the Sun.

    The calculations are performed from the point of view of a geocentric observer.

    Attributes
    ----------
    start : Time
        The date and time at the beginning of the transit.
        This is the moment the planet first becomes visible against the Sun in its background.
    peak : Time
        When the planet is most aligned with the Sun, as seen from the Earth.
    finish : Time
        The date and time at the end of the transit.
        This is the moment the planet is last seen against the Sun in its background.
    separation : float
        The minimum angular separation, in arcminutes, between the centers of the Sun and the planet.
        This angle pertains to the time stored in `peak`.
    """
    def __init__(self, start, peak, finish, separation):
        self.start = start
        self.peak = peak
        self.finish = finish
        self.separation = separation


def _PlanetShadowBoundary(context, time):
    (body, planet_radius_km, direction) = context
    shadow = _PlanetShadow(body, planet_radius_km, time)
    return direction * (shadow.r - shadow.p)


def _PlanetTransitBoundary(body, planet_radius_km, t1, t2, direction):
    # Search for the time the planet's penumbra begins/ends making contact with the center of the Earth.
    # context = new SearchContext_PlanetShadowBoundary(body, planet_radius_km, direction);
    tx = Search(_PlanetShadowBoundary, (body, planet_radius_km, direction), t1, t2, 1.0)
    if tx is None:
        raise Error('Planet transit boundary search failed')
    return tx


def SearchTransit(body, startTime):
    """Searches for the first transit of Mercury or Venus after a given date.

    Finds the first transit of Mercury or Venus after a specified date.
    A transit is when an inferior planet passes between the Sun and the Earth
    so that the silhouette of the planet is visible against the Sun in the background.
    To continue the search, pass the `finish` time in the returned structure to
    #NextTransit.

    Parameters
    ----------
    body : Body
        The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
    startTime : Time
        The date and time for starting the search for a transit.

    Returns
    -------
    TransitInfo
    """
    threshold_angle = 0.4     # maximum angular separation to attempt transit calculation
    dt_days = 1.0

    # Validate the planet and find its mean radius.
    if body == Body.Mercury:
        planet_radius_km = 2439.7
    elif body == Body.Venus:
        planet_radius_km = 6051.8
    else:
        raise InvalidBodyError()

    search_time = startTime
    while True:
        # Search for the next inferior conjunction of the given planet.
        # This is the next time the Earth and the other planet have the same
        # ecliptic longitude as seen from the Sun.
        conj = SearchRelativeLongitude(body, 0.0, search_time)

        # Calculate the angular separation between the body and the Sun at this time.
        conj_separation = AngleFromSun(body, conj)

        if conj_separation < threshold_angle:
            # The planet's angular separation from the Sun is small enough
            # to consider it a transit candidate.
            # Search for the moment when the line passing through the Sun
            # and planet are closest to the Earth's center.
            shadow = _PeakPlanetShadow(body, planet_radius_km, conj)
            if shadow.r < shadow.p:      # does the planet's penumbra touch the Earth's center?
                # Find the beginning and end of the penumbral contact.
                time_before = shadow.time.AddDays(-dt_days)
                start = _PlanetTransitBoundary(body, planet_radius_km, time_before, shadow.time, -1.0)
                time_after = shadow.time.AddDays(+dt_days)
                finish = _PlanetTransitBoundary(body, planet_radius_km, shadow.time, time_after, +1.0)
                min_separation = 60.0 * AngleFromSun(body, shadow.time)
                return TransitInfo(start, shadow.time, finish, min_separation)

        # This inferior conjunction was not a transit. Try the next inferior conjunction.
        search_time = conj.AddDays(10.0)


def NextTransit(body, prevTransitTime):
    """Searches for another transit of Mercury or Venus.

    After calling #SearchTransit to find a transit of Mercury or Venus,
    this function finds the next transit after that.
    Keep calling this function as many times as you want to keep finding more transits.

    Parameters
    ----------
    body : Body
        The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
    prevTransitTime : Time
        A date and time near the previous transit.

    Returns
    -------
    TransitInfo
    """
    startTime = prevTransitTime.AddDays(100.0)
    return SearchTransit(body, startTime)


class LibrationInfo:
    """Lunar libration angles, returned by #Libration.

    Contains lunar libration angles and lunar position information
    for a given moment in time. See #Libration for more details.

    Attributes
    ----------
    elat : float
        Sub-Earth libration ecliptic latitude angle, in degrees.
    elon : float
        Sub-Earth libration ecliptic longitude angle, in degrees.
    mlat : float
        Moon's geocentric ecliptic latitude.
    mlon : float
        Moon's geocentric ecliptic longitude.
    dist_km : float
        Distance between the centers of the Earth and Moon in kilometers.
    diam_deg : float
        The apparent angular diameter of the Moon as seen from the center of the Earth.
    """
    def __init__(self, elat, elon, mlat, mlon, dist_km, diam_deg):
        self.elat = elat
        self.elon = elon
        self.mlat = mlat
        self.mlon = mlon
        self.dist_km = dist_km
        self.diam_deg = diam_deg


def Libration(time):
    """Calculates the Moon's libration angles at a given moment in time.

    Libration is an observed back-and-forth wobble of the portion of the
    Moon visible from the Earth. It is caused by the imperfect tidal locking
    of the Moon's fixed rotation rate, compared to its variable angular speed
    of orbit around the Earth.

    This function calculates a pair of perpendicular libration angles,
    one representing rotation of the Moon in eclitpic longitude `elon`, the other
    in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.

    This function also returns the geocentric position of the Moon
    expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
    distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
    and the apparent angular diameter of the Moon `diam_deg`.

    Parameters
    ----------
    time : Time
        The date and time for which to calculate the Moon's libration angles.

    Returns
    -------
    LibrationInfo
    """
    t = time.tt / 36525.0
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    moon = _CalcMoon(time)
    mlon = moon.geo_eclip_lon
    mlat = moon.geo_eclip_lat
    dist_km = moon.distance_au * KM_PER_AU
    diam_deg = 2.0 * math.degrees(math.atan(_MOON_MEAN_RADIUS_KM / math.sqrt(dist_km*dist_km - _MOON_MEAN_RADIUS_KM*_MOON_MEAN_RADIUS_KM)))

    # Inclination angle
    I = math.radians(1.54242)

    # Moon's argument of latitude in radians.
    f = math.radians(_NormalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000))

    # Moon's ascending node's mean longitude in radians.
    omega = math.radians(_NormalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000))

    # Sun's mean anomaly.
    m = math.radians(_NormalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000))

    # Moon's mean anomaly.
    mdash = math.radians(_NormalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000))

    # Moon's mean elongation.
    d = math.radians(_NormalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000))

    # Eccentricity of the Earth's orbit.
    e = 1.0 - 0.002516*t - 0.0000074*t2

    # Optical librations
    w = mlon - omega
    a = math.atan2(math.sin(w)*math.cos(mlat)*math.cos(I) - math.sin(mlat)*math.sin(I), math.cos(w)*math.cos(mlat))
    ldash = _LongitudeOffset(math.degrees(a - f))
    bdash = math.asin(-math.sin(w)*math.cos(mlat)*math.sin(I) - math.sin(mlat)*math.cos(I))

    # Physical librations
    k1 = math.radians(119.75 + 131.849*t)
    k2 = math.radians(72.56 + 20.186*t)

    rho = (
        -0.02752*math.cos(mdash) +
        -0.02245*math.sin(f) +
        +0.00684*math.cos(mdash - 2*f) +
        -0.00293*math.cos(2*f) +
        -0.00085*math.cos(2*f - 2*d) +
        -0.00054*math.cos(mdash - 2*d) +
        -0.00020*math.sin(mdash + f) +
        -0.00020*math.cos(mdash + 2*f) +
        -0.00020*math.cos(mdash - f) +
        +0.00014*math.cos(mdash + 2*f - 2*d)
    )

    sigma = (
        -0.02816*math.sin(mdash) +
        +0.02244*math.cos(f) +
        -0.00682*math.sin(mdash - 2*f) +
        -0.00279*math.sin(2*f) +
        -0.00083*math.sin(2*f - 2*d) +
        +0.00069*math.sin(mdash - 2*d) +
        +0.00040*math.cos(mdash + f) +
        -0.00025*math.sin(2*mdash) +
        -0.00023*math.sin(mdash + 2*f) +
        +0.00020*math.cos(mdash - f) +
        +0.00019*math.sin(mdash - f) +
        +0.00013*math.sin(mdash + 2*f - 2*d) +
        -0.00010*math.cos(mdash - 3*f)
    )

    tau = (
        +0.02520*e*math.sin(m) +
        +0.00473*math.sin(2*mdash - 2*f) +
        -0.00467*math.sin(mdash) +
        +0.00396*math.sin(k1) +
        +0.00276*math.sin(2*mdash - 2*d) +
        +0.00196*math.sin(omega) +
        -0.00183*math.cos(mdash - f) +
        +0.00115*math.sin(mdash - 2*d) +
        -0.00096*math.sin(mdash - d) +
        +0.00046*math.sin(2*f - 2*d) +
        -0.00039*math.sin(mdash - f) +
        -0.00032*math.sin(mdash - m - d) +
        +0.00027*math.sin(2*mdash - m - 2*d) +
        +0.00023*math.sin(k2) +
        -0.00014*math.sin(2*d) +
        +0.00014*math.cos(2*mdash - 2*f) +
        -0.00012*math.sin(mdash - 2*f) +
        -0.00012*math.sin(2*mdash) +
        +0.00011*math.sin(2*mdash - 2*m - 2*d)
    )

    ldash2 = -tau + (rho*math.cos(a) + sigma*math.sin(a))*math.tan(bdash)
    bdash = math.degrees(bdash)
    bdash2 = sigma*math.cos(a) - rho*math.sin(a)
    return LibrationInfo(bdash + bdash2, ldash + ldash2, mlat, mlon, dist_km, diam_deg)
