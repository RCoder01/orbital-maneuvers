from . import lib
from . import inclination_change
import math

def resources(object_1: dict, object_2: dict) -> tuple[float, float]:
    orbit_1 = {'periapsis': float(object_1['periapsis']), 'apoapsis': float(object_1['apoapsis'])}
    orbit_2 = {'periapsis': float(object_2['periapsis']), 'apoapsis': float(object_2['apoapsis'])}
    # intermediate orbit needed for arg of periapsis change
    radius_i = None
    # not perfect, but good enough
    if orbit_1['apoapsis'] <= orbit_2['periapsis']:
        radius_i = orbit_2['periapsis']
    elif orbit_1['apoapsis'] >= orbit_2['apoapsis']:
        radius_i = orbit_2['apoapsis']
    else:
        radius_i = orbit_1['apoapsis']

    inc_1 = float(object_1['INCLINATION'])
    inc_2 = float(object_2['INCLINATION'])
    inc_delta = abs(inc_1 - inc_2)
    inc_change_dv = min(
        inclination_change.inc_change_dv(semimajor_axis=float(object_1['SEMIMAJOR_AXIS']), eccentricity=float(object_1['ECCENTRICITY']), inclination_delta=inc_delta),
        inclination_change.inc_change_dv(semimajor_axis=float(object_2['SEMIMAJOR_AXIS']), eccentricity=float(object_2['ECCENTRICITY']), inclination_delta=inc_delta),
        inclination_change.inc_change_dv(semimajor_axis=radius_i, eccentricity=0, inclination_delta=inc_delta))

    RAAN_1 = float(object_1['RA_OF_ASC_NODE'])
    RAAN_2 = float(object_2['RA_OF_ASC_NODE'])
    RAAN_delta = RAAN_1 - RAAN_2
    if RAAN_delta < 0:
        RAAN_delta += 360

    precession_1 = math.degrees(lib.nodal_precession(inc_1, float(object_1['SEMIMAJOR_AXIS'])))
    precession_i = math.degrees(lib.nodal_precession(inc_delta, radius_i))
    precession_2 = math.degrees(lib.nodal_precession(inc_2, float(object_2['SEMIMAJOR_AXIS'])))
    precession_delta_1 = precession_1 - precession_2
    precession_delta_i = precession_i - precession_2

    time_to_precess = min(
        ((RAAN_delta - 360) if precession_delta_1 < 0 else RAAN_delta) / precession_delta_1,
        ((RAAN_delta - 360) if precession_delta_i < 0 else RAAN_delta) / precession_delta_i)

    orbit_1_keplerian = lib.keplerian(**orbit_1)
    orbit_2_keplerian = lib.keplerian(**orbit_2)
    orbit_1_to_orbit_i = lib.planar_orbit_change_dv(orbit_1_keplerian['semimajor_axis'], radius_i, orbit_1_keplerian['eccentricity'], 0)
    orbit_i_to_orbit_2 = lib.planar_orbit_change_dv(radius_i, orbit_2_keplerian['semimajor_axis'], 0, orbit_2_keplerian['eccentricity'])

