import math
import constants

def gravitational_parameter(mass_1: float = constants.EARTH_MASS, mass_2: float = 0, /) -> float:
    return constants.GRAVITATIONAL_CONSTANT * (mass_1 + mass_2)

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
        semimajor_axis: float, # m
        eccentricity: float = 0, # unitless
        /,
        equitorial_radius: int = constants.EARTH_EQUITORIAL_RADIUS, # m
        gravitational_parameter: int = gravitational_parameter(), # kg
        j2: float = constants.EARTH_J2, # unitless
        ) -> float: # rad/s
    """
    Returns the approximate nodal precession rate of a satellite

    https://en.wikipedia.org/wiki/Nodal_precession#Rate_of_precession
    """
    return ((-3 * equitorial_radius**2 * j2 * angular_velocity(period(semimajor_axis, gravitational_parameter)) * math.cos(math.radians(inclination)))
     / (2 * (semimajor_axis * (1 - eccentricity**2))**2))

def semimajor_axis(*, mean_motion: float, gravitational_parameter: float = gravitational_parameter()) -> float:
    return (gravitational_parameter / (2 * math.pi * mean_motion) ** 2) ** (1 / 3)

def apoapsis(semimajor_axis: float, eccentricity: float = 0) -> float:
    return semimajor_axis * (1 + eccentricity)

def periapsis(semimajor_axis: float, eccentricity: float = 0) -> float:
    return semimajor_axis * (1 - eccentricity)

def keplerian(periapsis: float, apoapsis: float) -> dict[str, float]:
    return {'semimajor_axis': (periapsis + apoapsis) / 2, 'eccentricity': (apoapsis - periapsis) / (apoapsis + periapsis)}

def angular_momentum_per_kg(*, semimajor_axis: float, eccentricity: float = 0, gravitational_parameter: float = gravitational_parameter()) -> float:
    return (gravitational_parameter * semimajor_axis * (1 - eccentricity**2)) ** 0.5

def max_velocity(eccentricity: float, angular_momentum_per_kg: float, gravitational_parameter: float = gravitational_parameter()) -> float:
    return (1 + eccentricity) * gravitational_parameter / angular_momentum_per_kg

def min_velocity(eccentricity: float, angular_momentum_per_kg: float, gravitational_parameter: float = gravitational_parameter()) -> float:
    return (1 - eccentricity) * gravitational_parameter / angular_momentum_per_kg

def inc_change_dv(*, semimajor_axis: float, eccentricity: float, inclination_delta: float):
    v_min = min_velocity(eccentricity, angular_momentum_per_kg(semimajor_axis=semimajor_axis, eccentricity=eccentricity))
    return 2 * v_min * math.sin(math.radians(inclination_delta / 2))

def planar_orbit_change_dv(initial_semimajor_axis: float, final_semimajor_axis: float, initial_eccentricity: float = 0, final_eccentricity: float = 0) -> float:
    intermediate_keplerian_a = keplerian(periapsis(initial_semimajor_axis, initial_eccentricity), apoapsis(final_semimajor_axis, final_eccentricity))
    intermediate_keplerian_b = keplerian(periapsis(final_semimajor_axis, final_eccentricity), apoapsis(initial_semimajor_axis, initial_eccentricity))

    dv_a = abs(max_velocity(initial_eccentricity, angular_momentum_per_kg(semimajor_axis=initial_semimajor_axis, eccentricity=initial_eccentricity)) - max_velocity(intermediate_keplerian_a['eccentricity'], angular_momentum_per_kg(semimajor_axis=intermediate_keplerian_a['semimajor_axis'], eccentricity=intermediate_keplerian_a['eccentricity']))) + abs(min_velocity(final_eccentricity, angular_momentum_per_kg(semimajor_axis=final_semimajor_axis, eccentricity=final_eccentricity)) - min_velocity(intermediate_keplerian_a['eccentricity'], angular_momentum_per_kg(semimajor_axis=intermediate_keplerian_a['semimajor_axis'], eccentricity=intermediate_keplerian_a['eccentricity'])))
    dv_b = abs(min_velocity(initial_eccentricity, angular_momentum_per_kg(semimajor_axis=initial_semimajor_axis, eccentricity=initial_eccentricity)) - min_velocity(intermediate_keplerian_b['eccentricity'], angular_momentum_per_kg(semimajor_axis=intermediate_keplerian_b['semimajor_axis'], eccentricity=intermediate_keplerian_b['eccentricity']))) + abs(max_velocity(final_eccentricity, angular_momentum_per_kg(semimajor_axis=final_semimajor_axis, eccentricity=final_eccentricity)) - max_velocity(intermediate_keplerian_b['eccentricity'], angular_momentum_per_kg(semimajor_axis=intermediate_keplerian_b['semimajor_axis'], eccentricity=intermediate_keplerian_b['eccentricity'])))

    return min(dv_a, dv_b)
