import math
from typing import Self
import constants


def gravitational_parameter(mass_1: float = constants.EARTH_MASS, mass_2: float = 0, /) -> float:
    return constants.GRAVITATIONAL_CONSTANT * (mass_1 + mass_2) # m^3 s^-2

class Orbit2d:
    def __init__(self, semimajor_axis: float, eccentricity: float = 0, gravitational_parameter: float = gravitational_parameter()):
        self.semimajor_axis = semimajor_axis # m
        self.eccentricity = eccentricity # unitless
        self.gravitational_parameter = gravitational_parameter # m^3 s^-2

    def __repr__(self):
        return f"Ellipse(semimajor_axis={self.semimajor_axis}, eccentricity={self.eccentricity})"

    @classmethod
    def from_dict(cls, orbit: dict, **kwargs) -> Self:
        return cls(semimajor_axis=float(orbit['SEMIMAJOR_AXIS']) * 1000, eccentricity=float(orbit['ECCENTRICITY']), **kwargs)

    @classmethod
    def from_mean_motion(cls, mean_motion: float, eccentricity: float, gravitational_parameter: float = gravitational_parameter()) -> Self:
        return cls(
            semimajor_axis=(gravitational_parameter / (2 * math.pi * mean_motion) ** 2) ** (1 / 3),
            eccentricity=eccentricity,
            gravitational_parameter=gravitational_parameter)

    @classmethod
    def from_apsis(cls, /, periapsis: float, apoapsis: float, **kwargs) -> Self:
        periapsis, apoapsis = sorted([periapsis, apoapsis])
        return cls(semimajor_axis=(periapsis + apoapsis) / 2, eccentricity=(apoapsis - periapsis) / (apoapsis + periapsis), **kwargs)

    @classmethod
    def from_velocities(cls, periapsis_velocity: float, apoapsis_velocity: float, gravitational_parameter: float = gravitational_parameter()) -> Self:
        # semimajor_axis = gravitational_parameter / (apoapsis_velocity * periapsis_velocity)
        periapsis = 2 * gravitational_parameter / ((apoapsis_velocity * periapsis_velocity) * (1 + (periapsis_velocity / apoapsis_velocity)))
        apoapsis = periapsis * periapsis_velocity / apoapsis_velocity
        return cls.from_apsis(periapsis, apoapsis, gravitational_parameter=gravitational_parameter)

    @classmethod
    def from_periapsis(cls, periapsis_distance: float, periapsis_velocity: float, gravitational_parameter: float = gravitational_parameter()) -> Self:
        semimajor_axis = -1 * gravitational_parameter * periapsis_distance / (periapsis_velocity ** 2 * periapsis_distance - 2 * gravitational_parameter)
        apoapsis = 2 * semimajor_axis - periapsis_distance
        return cls.from_apsis(periapsis_distance, apoapsis, gravitational_parameter=gravitational_parameter)

    @classmethod
    def from_apoapsis(cls, apoapsis_distance: float, apoapsis_velocity: float, gravitational_parameter: float = gravitational_parameter()) -> Self:
        semimajor_axis = -1 * gravitational_parameter * apoapsis_distance / (apoapsis_velocity ** 2 * apoapsis_distance - 2 * gravitational_parameter)
        periapsis = 2 * semimajor_axis - apoapsis_distance
        return cls.from_apsis(periapsis, apoapsis_distance, gravitational_parameter=gravitational_parameter)

    def periapsis_burn(self, delta_v: float) -> Self:
        return self.from_periapsis(self.periapsis, self.periapsis_velocity + delta_v)

    def apoapsis_burn(self, delta_v: float) -> Self:
        return self.from_apoapsis(self.apoapsis, self.apoapsis_velocity + delta_v)

    @property
    def semiminor_axis(self) -> float:
        return self.semimajor_axis * math.sqrt(1 - self.eccentricity ** 2)

    @property
    def specific_angular_momentum(self) -> float:
        return (self.gravitational_parameter * self.semimajor_axis * (1 - self.eccentricity**2)) ** 0.5

    @property
    def periapsis(self) -> float:
        return self.semimajor_axis * (1 - self.eccentricity)

    @property
    def apoapsis(self) -> float:
        return self.semimajor_axis * (1 + self.eccentricity)

    @property
    def periapsis_velocity(self) -> float:
        return ((self.gravitational_parameter * self.semiminor_axis ** 2) / (self.semimajor_axis * self.periapsis ** 2)) ** 0.5
        # return (1 - self.eccentricity) * self.gravitational_parameter / self.specific_angular_momentum

    @property
    def apoapsis_velocity(self) -> float:
        return ((self.gravitational_parameter * self.semiminor_axis ** 2) / (self.semimajor_axis * self.apoapsis ** 2)) ** 0.5
        # return (1 + self.eccentricity) * self.gravitational_parameter / self.specific_angular_momentum

    @property
    def period(self) -> float:
        return ((self.semimajor_axis**3 * 4 * math.pi**2) / self.gravitational_parameter) ** 0.5
    
    @property
    def angular_velocity(self) -> float:
        return 2 * math.pi / self.period

    def apoapsis_altitude(self, body_radius: float = constants.EARTH_MEAN_RADIUS) -> float:
        return self.apoapsis - body_radius
    
    def periapsis_altitude(self, body_radius: float = constants.EARTH_MEAN_RADIUS) -> float:
        return self.periapsis - body_radius

    def __eq__(self, other: Self) -> bool:
        return self.semimajor_axis - other.semimajor_axis < 1e-6 and self.eccentricity - other.eccentricity < 1e-6


def period(
        semimajor_axis: float, # m
        gravitational_parameter: int = gravitational_parameter() # m^3 s^-2
        ) -> float: # s
    """
    Period of an object in elliptical orbit, based on Kepler's third law

    https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Third_law
    """
    return ((semimajor_axis**3 * 4 * math.pi**2) / gravitational_parameter) ** 0.5


def angular_velocity(
        period: float, # s
        ) -> float: # rad/s
    return 2 * math.pi / period


def nodal_precession(
        inclination: float, # deg
        orbit: Orbit2d,
        *,
        equitorial_radius: int = constants.EARTH_EQUITORIAL_RADIUS, # m
        j2: float = constants.EARTH_J2, # unitless
        ) -> float: # rad/s
    """
    Returns the approximate nodal precession rate of a satellite

    https://en.wikipedia.org/wiki/Nodal_precession#Rate_of_precession
    """
    return ((-3 * equitorial_radius**2 * j2 * orbit.angular_velocity * math.cos(math.radians(inclination)))
     / (2 * (orbit.semimajor_axis * (1 - orbit.eccentricity**2))**2))


def inclination_change_dv(orbit: Orbit2d, inclination_delta: float) -> float:
    return 2 * orbit.apoapsis_velocity * math.sin(math.radians(inclination_delta / 2))


def tangential_orbit_change_dv(initial: Orbit2d, final: Orbit2d) -> float:
    if abs(initial.periapsis - final.periapsis) <= 1e-6:
        return abs(initial.periapsis_velocity - final.periapsis_velocity)
    if abs(initial.periapsis - final.apoapsis) <= 1e-6:
        return abs(initial.periapsis_velocity - final.apoapsis_velocity)
    if abs(initial.apoapsis - final.periapsis) <= 1e-6:
        return abs(initial.apoapsis_velocity - final.periapsis_velocity)
    if abs(initial.apoapsis - final.apoapsis) <= 1e-6:
        return abs(initial.apoapsis_velocity - final.apoapsis_velocity)

    intermediate_a = Orbit2d.from_apsis(initial.periapsis, final.apoapsis)
    intermediate_b = Orbit2d.from_apsis(final.periapsis, initial.apoapsis)

    return min(
        tangential_orbit_change_dv(initial, intermediate_a) + tangential_orbit_change_dv(intermediate_a, final),
        tangential_orbit_change_dv(initial, intermediate_b) + tangential_orbit_change_dv(intermediate_b, final))


def arg_of_periapsis_change_dv(orbit: Orbit2d, arg_of_periapsis_delta: float) -> float:
    return 2 * orbit.eccentricity * (orbit.gravitational_parameter / (orbit.semimajor_axis * (1 - orbit.eccentricity**2)))**0.5 * math.sin(math.radians(arg_of_periapsis_delta)/2)
