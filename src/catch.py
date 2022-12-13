from functools import partial
import lib
import math
import constants


class Writer:
    def __init__(self):
        self.string = ''
    
    def write(self, s):
        self.string += s



def resources_to_transfer(
        object_1: dict,
        object_2: dict,
        time_offset: float = 0,
        extra_budget: float = 0,
        debug_log = None
        ) -> tuple[float, float]:
    if debug_log is None:
        debug_log = Writer()

    # set up orbits
    orbit_1 = lib.Orbit2d.from_dict(object_1)
    orbit_2 = lib.Orbit2d.from_dict(object_2)
    orbit_i = find_intermediate_orbit(orbit_1, orbit_2)
    orbit_j = find_extra_budget_orbit(extra_budget, orbit_1, orbit_2)
    print(f'Changing from {orbit_1.semimajor_axis/1000:.0f} km to {orbit_2.semimajor_axis/1000:.0f} km', file=debug_log)

    inc_1 = float(object_1['INCLINATION'])
    inc_2 = float(object_2['INCLINATION'])

    nodal_precession_1 = lib.nodal_precession(inc_1, orbit_1)
    nodal_precession_2 = lib.nodal_precession(inc_2, orbit_2)

    # Calculate time to match RAAN
    RAAN_1 = (float(object_1['RA_OF_ASC_NODE']) + nodal_precession_1 * time_offset) % 360
    RAAN_2 = (float(object_2['RA_OF_ASC_NODE']) + nodal_precession_2 * time_offset) % 360
    RAAN_delta = RAAN_1 - RAAN_2
    if RAAN_delta < 0:
        RAAN_delta += 360
    print(f'initial RAAN: {RAAN_1}, target RAAN: {RAAN_2}, delta: {RAAN_delta}', file=debug_log)

    time_to_precess_to_2 = partial(match_RAAN, target_precession_rate=nodal_precession_2, RAAN_delta=RAAN_delta)

    j_precession_time = time_to_precess_to_2(orbit_j, inc_2)

    min_time_to_precess, period = min(
        (time_to_precess_to_2(orbit_1, inc_1), orbit_1.period),
        (time_to_precess_to_2(orbit_i, inc_2), orbit_i.period),
        (j_precession_time, orbit_j.period),
        key=lambda x: x[0])

    # == comaprison is ok because min returns the same object
    if period == orbit_1.period:
        print(f'match RAAN at 1, {min_time_to_precess:.0f} s', file=debug_log)
    if period == orbit_i.period:
        print(f'match RAAN at i, {min_time_to_precess:.0f} s', file=debug_log)
    if period == orbit_j.period:
        print(f'match RAAN at j, {min_time_to_precess:.0f} s', file=debug_log)

    orbits_to_precess = min_time_to_precess / period

    j_used = min_time_to_precess is j_precession_time
    j_change_dv_total = extra_budget if j_used else 0

    # Calculate cheapest place to match inclination
    inc_delta = abs(inc_1 - inc_2)
    necessary_inc_change_dv = partial(lib.inclination_change_dv, inclination_delta=inc_delta)
    inc_change_opportunities = [
        necessary_inc_change_dv(orbit_1),
        necessary_inc_change_dv(orbit_2),
        necessary_inc_change_dv(orbit_i)]
    if j_used:
        inc_change_opportunities.append(necessary_inc_change_dv(orbit_j))
    inc_change_dv = min(inc_change_opportunities)

    print(f'inc change from {inc_1} to {inc_2} (delta {inc_delta}) requires {inc_change_dv} m/s', file=debug_log)
    print(f'change inc at {["1", "2", "i", "j"][inc_change_opportunities.index(inc_change_dv)]}', file=debug_log)

    # Calculate time to match mean anomaly
    orbits_to_match_mean_anomaly = match_mean_anomaly(orbit_2, period, orbits_to_precess)
    print(f'time to match mean anomaly: {orbits_to_match_mean_anomaly * period} s', file=debug_log)

    # Calculate hohmann delta v
    orbit_1_to_i = lib.coaxial_elliptic_orbit_change_dv(orbit_1, orbit_i)
    orbit_i_to_2 = lib.coaxial_elliptic_orbit_change_dv(orbit_i, orbit_2)
    print(f'orbit_1: {orbit_1.semimajor_axis/1000:.0f} km', file=debug_log)
    if j_used:
        print(f'orbit_j: {orbit_j.semimajor_axis/1000:.0f} km (1 to j to 1: {extra_budget:.2f} m/s)', file=debug_log)
    print(f'orbit_i: {orbit_i.semimajor_axis/1000:.0f} km (1 to i: {orbit_1_to_i:.2f} m/s)', file=debug_log)
    print(f'orbit_2: {orbit_2.semimajor_axis/1000:.0f} km (i to 2: {orbit_i_to_2:.2f} m/s)', file=debug_log)


    total_dv = orbit_1_to_i + orbit_i_to_2 + inc_change_dv + j_change_dv_total
    total_time = orbits_to_precess * period + orbits_to_match_mean_anomaly * orbit_2.period
    return total_dv, total_time, debug_log


def match_RAAN(
        initial_orbit: lib.Orbit2d,
        inclination: float,
        target_precession_rate: float,
        RAAN_delta: float
        ) -> float:
    """
    Calcualtes orbits needed to precess from one RAAN to another.
    """
    initial_precession_rate = lib.nodal_precession(inclination, initial_orbit)
    precession_rate_delta = target_precession_rate - initial_precession_rate
    if precession_rate_delta < 0:
        RAAN_delta -= 360
    if precession_rate_delta:
        return RAAN_delta / precession_rate_delta
    return float('inf')


def find_extra_budget_orbit(extra_budget: float, orbit_1: lib.Orbit2d, orbit_2: lib.Orbit2d) -> lib.Orbit2d:
    orbit_j: lib.Orbit2d
    if orbit_1.period < orbit_2.period:
        orbit_j = orbit_1.apoapsis_burn(-extra_budget / 2)
    else:
        orbit_j = orbit_2.periapsis_burn(extra_budget / 2)
    return orbit_j


def find_intermediate_orbit(orbit_1: lib.Orbit2d, orbit_2: lib.Orbit2d) -> lib.Orbit2d:
    """
    Intermediate orbit needed for arg of periapsis change

    Simpler to just use a circular intermediate orbit than periapse change
    because most orbits are roughly circular

    Not perfect, but good enough
    """
    if orbit_1.apoapsis <= orbit_2.periapsis:
        radius = orbit_2.periapsis
    elif orbit_1.apoapsis >= orbit_2.apoapsis:
        radius = orbit_2.apoapsis
    else:
        radius = orbit_1.apoapsis
    return lib.Orbit2d(radius)


def match_mean_anomaly(orbit: lib.Orbit2d, period: float, orbits_to_precess: float) -> float:
    """
    Calculates number of orbits to match mean anomaly,
    assuming both orbits are 360 degrees separated
    """
    precession = (period - orbit.period) / orbit.period
    return math.fabs(1 / precession)


def get_objects():
    """
    Returns a list of objects from leo_debris.json
    """
    import json

    try:
        with open('../data/leo_debris.json', mode='r', encoding='UTF-8') as f:
            objects = json.load(f)
    except FileNotFoundError:
        print('leo_debris.json not found, ensure you have run get.py')
    return objects


def collect(
        objects: list[dict],
        start: int,
        per_catch_fuel_budget: float,
        per_catch_time_target: float,
        total_fuel_budget: float
        ) -> tuple[list[dict], float, float, list[tuple[float, float, int]]]:
    """
    Collects objects continuously until fuel budget is exhausted
    """
    v = t = 0
    caught = [objects[start]]
    metadata = []
    while v < total_fuel_budget:
        try:
            catch, dv, dt, index, debug_info = collect_one(objects, caught, start, per_catch_fuel_budget, per_catch_time_target, t)
        except ValueError:
            break
        v += dv
        t += dt
        caught.append(catch)
        metadata.append((dv, dt, index, debug_info))
    return caught, v, t, metadata


def collect_one(
        objects: list[dict],
        caught: list[dict],
        start: int,
        per_catch_fuel_budget: float,
        per_catch_time_target: float,
        current_time: float):
    """
    Collects one object

    Makes sure dt is within target time
    """
    index = start
    while index < len(objects) - 1:
        index += 1
        if objects[index] in caught:
            continue
        try:
            dv, dt, debug_info = resources_to_transfer(caught[-1], objects[index], current_time, per_catch_fuel_budget)
        except ZeroDivisionError:
            continue
        if dt < per_catch_time_target:
            return objects[index], dv, dt, index, debug_info
    raise ValueError('No object found')


def deorbit_dv(orbit: lib.Orbit2d) -> float:
    """
    Calculates delta-v needed to deorbit an object to the ground
    """
    return lib.coaxial_elliptic_orbit_change_dv(
        orbit,
        lib.Orbit2d.from_apsides(orbit.apoapsis, constants.EARTH_MEAN_RADIUS))


if __name__ == '__main__':
    objects = get_objects()
    # filter out eccentric orbits otherwise cost to circularize is high
    s_objects = sorted(
        list(obj for obj in objects if float(obj['ECCENTRICITY']) < 0.007 and 400 < float(obj['APOAPSIS']) < 600),
        key=lambda obj: (float(obj['INCLINATION']), float(obj['RA_OF_ASC_NODE'])))

    import sys

    ### CHANGE THESE ###
    PER_CATCH_FUEL_BUDGET = 100 # m/s
    PER_CATCH_TIME_BUDGET = 10**7.7 # s
    TOTAL_FUEL_BUDGET = 1100 # m/s
    ####################
    try:
        START = int(sys.argv[1])
    except IndexError:
        print('No start index given')
        print('Please enter a number between 0 and', len(s_objects) - 1)
        print('Example: py -3.11', sys.argv[0], '0')
        exit(1)

    caught, v, t, meta = collect(s_objects, START, PER_CATCH_FUEL_BUDGET, PER_CATCH_TIME_BUDGET, TOTAL_FUEL_BUDGET)

    for dv, dt, i, debug_info in meta:
        print(f'{i:3}: {dv:3.0f} m/s, 10^{math.log10(dt):4.2f} s')
    print(f'cumulative ({len(caught) - 1}): {v:3.0f} m/s, {t:.0f} s (10^{math.log10(t):4.2f} s, {t/60/60/24/365:5.2f} years)')

    dv_to_deorbit = deorbit_dv(lib.Orbit2d.from_dict(caught[-1]))
    print(f'deorbit dv: {dv_to_deorbit:.2f} m/s')
    print(f'total dv: {v + dv_to_deorbit:.2f} m/s')
