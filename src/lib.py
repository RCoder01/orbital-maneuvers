import math
import constants

def period(
        semimajor_axis: float, # m
        combined_mass: int = constants.EARTH_MASS # kg
        ) -> float: # s
    """
    Period of an object in elliptical orbit, based on Kepler's third law

    https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Third_law
    """
    return ((semimajor_axis**3 * 4 * math.pi**2) / (constants.GRAVITATIONAL_CONSTANT * combined_mass)) ** 0.5

def angular_velocity(
        period: float, # s
        ) -> float: # rad/s
    return 2 * math.pi / period

def nodal_precession(
        inclination: float, # deg
        semimajor_axis: float, # m
        eccentricity: float = 0, # unitless
        /,
        equitorial_radius: int = constants.EARTH_EQUITORIAL_RADIUS, # m
        combined_mass: int = constants.EARTH_MASS, # kg
        j2: float = constants.EARTH_J2, # unitless
        ) -> float: # rad/s
    """
    Returns the approximate nodal precession rate of a satellite

    https://en.wikipedia.org/wiki/Nodal_precession#Rate_of_precession
    """
    return ((-3 * equitorial_radius**2 * j2 * angular_velocity(period(semimajor_axis, combined_mass)) * math.cos(math.radians(inclination)))
     / (2 * (semimajor_axis * (1 - eccentricity**2))**2))

def semimajor_axis(*, mean_motion: float, combined_mass: int = constants.EARTH_MASS) -> float:
    if mean_motion is not None:
        return (constants.GRAVITATIONAL_CONSTANT * combined_mass / (2 * math.pi * mean_motion) ** 2) ** (1 / 3)
    else:
        raise ValueError('mean_motion must be specified')

def apoapsis(semimajor_axis: float, eccentricity: float = 0) -> float:
    return semimajor_axis * (1 + eccentricity)

def periapsis(semimajor_axis: float, eccentricity: float = 0) -> float:
    return semimajor_axis * (1 - eccentricity)
