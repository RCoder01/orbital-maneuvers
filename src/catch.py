from functools import partial
import lib
import math
import constants

def resources_to_transfer(
        object_1: dict,
        object_2: dict,
        time_offset: float = 0,
        extra_budget: float = 0
        ) -> tuple[float, float]:
    orbit_1 = lib.Orbit2d.from_dict(object_1)
    orbit_2 = lib.Orbit2d.from_dict(object_2)
    orbit_i = find_intermediate_orbit(orbit_1, orbit_2)
    orbit_j = find_extra_budget_orbit(extra_budget, orbit_1, orbit_2)


    inc_1 = float(object_1['INCLINATION'])
    inc_2 = float(object_2['INCLINATION'])

    nodal_precession_1 = lib.nodal_precession(inc_1, orbit_1)
    nodal_precession_2 = lib.nodal_precession(inc_2, orbit_2)

    RAAN_1 = (float(object_1['RA_OF_ASC_NODE']) + nodal_precession_1 * time_offset) % 360
    RAAN_2 = (float(object_2['RA_OF_ASC_NODE']) + nodal_precession_2 * time_offset) % 360
    RAAN_delta = RAAN_1 - RAAN_2
    if RAAN_delta < 0:
        RAAN_delta += 360

    time_to_precess_to_2 = partial(match_RAAN, target_precession_rate=nodal_precession_2, RAAN_delta=RAAN_delta)

    j_precession_time = time_to_precess_to_2(orbit_j, inc_2)

    min_time_to_precess, period = min(
        (time_to_precess_to_2(orbit_1, inc_1), orbit_1.period),
        (time_to_precess_to_2(orbit_i, inc_2), orbit_i.period),
        (j_precession_time, orbit_j.period),
        key=lambda x: x[0])

    orbits_to_precess = min_time_to_precess / period

    j_used = min_time_to_precess is j_precession_time
    j_change_dv_total = extra_budget if j_used else 0

    inc_delta = abs(inc_1 - inc_2)
    necessary_inc_change_dv = partial(lib.inclination_change_dv, inclination_delta=inc_delta)
    inc_change_opportunities = [
        necessary_inc_change_dv(orbit_1),
        necessary_inc_change_dv(orbit_2),
        necessary_inc_change_dv(orbit_i)]
    if j_used:
        inc_change_opportunities.append(necessary_inc_change_dv(orbit_j))
    inc_change_dv = min(inc_change_opportunities)

    orbits_to_match_mean_anomaly = match_mean_anomaly(orbit_2, period, orbits_to_precess)

    orbit_1_to_i = lib.tangential_orbit_change_dv(orbit_1, orbit_i)
    orbit_i_to_2 = lib.tangential_orbit_change_dv(orbit_i, orbit_2)

    total_dv = orbit_1_to_i + orbit_i_to_2 + inc_change_dv + j_change_dv_total
    total_time = orbits_to_precess * period + orbits_to_match_mean_anomaly * orbit_2.period
    return total_dv, total_time


def match_RAAN(
        initial_orbit: lib.Orbit2d,
        inclination: float,
        target_precession_rate: float,
        RAAN_delta: float
        ) -> float:
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
    '''
    Intermediate orbit needed for arg of periapsis change

    Simpler to just use a circular intermediate orbit than periapse change because most orbits are roughly circular

    Not perfect, but good enough
    '''
    if orbit_1.apoapsis <= orbit_2.periapsis:
        radius = orbit_2.periapsis
    elif orbit_1.apoapsis >= orbit_2.apoapsis:
        radius = orbit_2.apoapsis
    else:
        radius = orbit_1.apoapsis
    return lib.Orbit2d(radius)


def match_mean_anomaly(orbit: lib.Orbit2d, period: float, orbits_to_precess: float) -> float:
    precession = (period - orbit.period) / orbit.period
    return math.fabs(1 / precession)
    precession_delta = precession * orbits_to_precess % 360
    total_delta = 360 # worst case scenario
    final_delta = (total_delta + precession_delta) % 360
    if final_delta < 0:
        final_delta += 360 # normalize to 0-360
    if precession < 0:
        final_delta -= 360 # ensure precession has the same sign as delta
    return final_delta / precession


def get_objects():
    import json

    with open('../data/leo_debris.json', mode='r', encoding='UTF-8') as f:
        objects = json.load(f)
    return objects


def collect(
        objects: list[dict],
        start: int,
        per_catch_fuel_budget: float,
        per_catch_time_target: float,
        total_fuel_budget: float
        ) -> tuple[list[dict], float, float, list[tuple[float, float, int]]]:
    v = t = 0
    caught = [objects[start]]
    metadata = []
    while v < total_fuel_budget:
        try:
            catch, dv, dt, index = collect_one(objects, caught, start, per_catch_fuel_budget, per_catch_time_target, t)
        except ValueError:
            break
        v += dv
        t += dt
        caught.append(catch)
        metadata.append((dv, dt, index))
    return caught, v, t, metadata


def collect_one(
        objects: list[dict],
        caught: list[dict],
        start: int,
        per_catch_fuel_budget: float,
        per_catch_time_target: float,
        current_time: float):
    index = start
    while index < len(objects) - 1:
        index += 1
        if objects[index] in caught:
            continue
        dv, dt = resources_to_transfer(caught[-1], objects[index], current_time, per_catch_fuel_budget)
        if dt < per_catch_time_target:
            return objects[index], dv, dt, index
    raise ValueError('No object found')


def deorbit_dv(orbit: lib.Orbit2d) -> float:
    return lib.tangential_orbit_change_dv(
        orbit,
        lib.Orbit2d.from_apsides(orbit.apoapsis, constants.EARTH_MEAN_RADIUS + 200_000))


if __name__ == '__main__':
    objects = get_objects()
    s_objects = sorted(
        list(obj for obj in objects if float(obj['ECCENTRICITY']) < 0.007 and 400 < float(obj['APOAPSIS']) < 600),
        key=lambda obj: (float(obj['INCLINATION']), float(obj['RA_OF_ASC_NODE'])))

    import sys

    PER_CATCH_FUEL_BUDGET = 300 # m/s
    PER_CATCH_TIME_BUDGET = 10**7.25 # s
    TOTAL_FUEL_BUDGET = 5000 # m/s
    START = int(sys.argv[1])

    caught, v, t, meta = collect(s_objects, START, PER_CATCH_FUEL_BUDGET, PER_CATCH_TIME_BUDGET, TOTAL_FUEL_BUDGET)

    for dv, dt, i in meta:
        print(f'{i:3}: {dv:3.0f} m/s, 10^{math.log10(dt):4.2f} s')
    print(f'cumulative ({len(caught) - 1}): {v:3.0f} m/s, {t:.0f} s (10^{math.log10(t):4.2f} s, {t/60/60/24/365:5.2f} years)')

    print(f'deorbit dv: {deorbit_dv(lib.Orbit2d.from_dict(caught[-1])):.2f} m/s')

    for _, _, i in [(0, 0, START)] + meta:
        print(f'{i}: {s_objects[i]["OBJECT_ID"]:13} {s_objects[i]["OBJECT_NAME"]}')
