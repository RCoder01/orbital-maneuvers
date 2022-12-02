import lib
import math
import constants

def resources(object_1: dict, object_2: dict, time_offset: float = 0, extra_budget: float = 0) -> tuple[float, float]:
    orbit_1 = lib.Orbit2d.from_dict(object_1)
    orbit_2 = lib.Orbit2d.from_dict(object_2)
    # intermediate orbit needed for arg of periapsis change
    # simpler to just use a circular intermediate orbit than periapse change because most orbits are roughly circular
    # not perfect, but good enough
    radius_i: float
    if orbit_1.apoapsis <= orbit_2.periapsis:
        radius_i = orbit_2.periapsis
    elif orbit_1.apoapsis >= orbit_2.apoapsis:
        radius_i = orbit_2.apoapsis
    else:
        radius_i = orbit_1.apoapsis
    orbit_i = lib.Orbit2d(radius_i)


    orbit_j: lib.Orbit2d
    if orbit_1.period < orbit_2.period:
        orbit_j = orbit_1.apoapsis_burn(-extra_budget / 2)
    else:
        orbit_j = orbit_2.periapsis_burn(extra_budget / 2)
    j_change_dv_total = extra_budget


    inc_1 = float(object_1['INCLINATION'])
    inc_2 = float(object_2['INCLINATION'])
    inc_delta = abs(inc_1 - inc_2)

    nodal_precession_2 = math.degrees(lib.nodal_precession(inc_2, orbit_2))

    nodal_precession_1 = math.degrees(lib.nodal_precession(inc_1, orbit_1))
    nodal_precession_i = math.degrees(lib.nodal_precession(inc_2, orbit_i))
    nodal_precession_j = math.degrees(lib.nodal_precession(inc_2, orbit_j))
    nodal_precession_delta_1 = nodal_precession_2 - nodal_precession_1
    nodal_precession_delta_i = nodal_precession_2 - nodal_precession_i
    nodal_precession_delta_j = nodal_precession_2 - nodal_precession_j

    RAAN_1 = (float(object_1['RA_OF_ASC_NODE']) + nodal_precession_1 * time_offset) % 360
    RAAN_2 = (float(object_2['RA_OF_ASC_NODE']) + nodal_precession_2 * time_offset) % 360
    RAAN_delta = RAAN_1 - RAAN_2
    if RAAN_delta < 0:
        RAAN_delta += 360

    time_to_precess_1 = ((RAAN_delta - 360) if nodal_precession_delta_1 < 0 else RAAN_delta) / nodal_precession_delta_1 if nodal_precession_delta_1 else float('inf')
    time_to_precess_i = ((RAAN_delta - 360) if nodal_precession_delta_i < 0 else RAAN_delta) / nodal_precession_delta_i if nodal_precession_delta_i else float('inf')
    time_to_precess_j = ((RAAN_delta - 360) if nodal_precession_delta_j < 0 else RAAN_delta) / nodal_precession_delta_j if nodal_precession_delta_j else float('inf')

    period: float
    orbits_to_precess: float
    j_used = False
    min_time = min(time_to_precess_1, time_to_precess_i, time_to_precess_j)
    if min_time is time_to_precess_1:
        period = orbit_1.period
        orbits_to_precess = time_to_precess_1 / period
    elif min_time is time_to_precess_i:
        period = orbit_i.period
        orbits_to_precess = time_to_precess_i / period
    else:
        period = orbit_j.period
        orbits_to_precess = time_to_precess_j / period
        j_used = True
        # print(f'j used {math.log10(time_to_precess_1):4.2f}, {math.log10(time_to_precess_i):4.2f}, {math.log10(time_to_precess_j):4.2f}')


    if not j_used:
        j_change_dv_total = 0

    inc_change_opportunities = [
        lib.inclination_change_dv(orbit_1, inc_delta),
        lib.inclination_change_dv(orbit_2, inc_delta),
        lib.inclination_change_dv(orbit_i, inc_delta)]
    if j_used:
        inc_change_opportunities.append(lib.inclination_change_dv(orbit_j, inc_delta))
    inc_change_dv = min(inc_change_opportunities)


    mean_anomaly_precession = (period - orbit_2.period) / orbit_2.period
    mean_anomaly_precession_delta = mean_anomaly_precession * orbits_to_precess % 360
    mean_anomaly_delta = 360 # worst case scenario #float(object_2['MEAN_ANOMALY']) - float(object_1['MEAN_ANOMALY'])
    final_mean_anomaly_delta = (mean_anomaly_delta + mean_anomaly_precession_delta) % 360
    if final_mean_anomaly_delta < 0:
        final_mean_anomaly_delta += 360 # normalize to 0-360
    if mean_anomaly_precession < 0:
        final_mean_anomaly_delta -= 360 # ensure precession has the same sign as delta
    orbits_to_match_mean_anomaly = final_mean_anomaly_delta / mean_anomaly_precession

    orbit_1_to_i = lib.tangential_orbit_change_dv(orbit_1, orbit_i)
    orbit_i_to_2 = lib.tangential_orbit_change_dv(orbit_i, orbit_2)

    # print(f'{orbit_1_to_i=}, {orbit_i_to_2=}, {inc_change_dv=}, {lib.tangential_orbit_change_dv(orbit_1, orbit_2)}')
    total_dv = orbit_1_to_i + orbit_i_to_2 + inc_change_dv + j_change_dv_total
    total_time = orbits_to_precess * period + orbits_to_match_mean_anomaly * orbit_2.period
    return total_dv, total_time

def get_objects():
    import json

    with open('../data/leo_debris.json', mode='r', encoding='UTF-8') as f:
        objects = json.load(f)
    return objects

def collect(objects, start, per_catch_budget):
    v = 0
    t = 0

    index = start
    caught = [objects[start]]
    metadata = []
    num = 10
    while v < 5000 and index < len(objects) - 1:
        # print(i, s_objects[i])
        index += 1
        if objects[index] in caught:
            continue
        dv, dt = resources(caught[-1], objects[index], t, per_catch_budget)
        if dt < 10**7.25:
            v += dv
            t += dt
            caught.append(objects[index])
            metadata.append((v, t, index))
            # print(f'{index}: {dv:3.0f} m/s, 10^{math.log10(dt):4.2f} s')
            index = start
    # print(f'cumulative ({len(caught) - 1}): {v:3.0f} m/s, {t:.0f} s (10^{math.log10(t):4.2f} s, {t/60/60/24/365:5.2f} years)')
    return caught, v, t, metadata


if __name__ == '__main__':
    objects = get_objects()
    s_objects = sorted(
        list(obj for obj in objects if float(obj['ECCENTRICITY']) < 0.007 and 400 < float(obj['APOAPSIS']) < 600),
        key=lambda obj: (float(obj['INCLINATION']), float(obj['RA_OF_ASC_NODE'])))

    PER_CATCH_BUDGET = 300 # m/s
    START = 370
    caught, v, t, meta = collect(s_objects, START, PER_CATCH_BUDGET)
    for dv, dt, i in meta:
        print(f'{i:3}: {dv:3.0f} m/s, 10^{math.log10(dt):4.2f} s')
    print(f'cumulative ({len(caught) - 1}): {v:3.0f} m/s, {t:.0f} s (10^{math.log10(t):4.2f} s, {t/60/60/24/365:5.2f} years)')
