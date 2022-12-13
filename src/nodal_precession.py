import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

from lib import nodal_precession, Orbit2d
import constants


def main():
    altitudes = range(400_000, 1600_000, 1_000)
    alt_delta = 100_000
    inclinations = range(0, 180, 1)
    inc_delta = 15


    omegas = np.zeros((len(inclinations), len(altitudes)))
    for x, inc in enumerate(inclinations):
        for y, alt in enumerate(altitudes):
            omegas[x, y] = math.radians(nodal_precession(inc, Orbit2d(alt + constants.EARTH_MEAN_RADIUS)))


    alt_ticks = range(0, len(altitudes) + 1, alt_delta // 1000)
    inc_ticks = range(0, len(inclinations) + 1, inc_delta)

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    im = plt.imshow(
        omegas,
        aspect='auto',
        extent=(alt_ticks.start, alt_ticks.stop - 1, inc_ticks.stop - 1, inc_ticks.start),
        cmap='seismic')
    plt.gca().invert_yaxis()

    plt.title('Nodal Precession in Circular Orbits as a Function of Inclination and Altitude', fontsize=12, pad=15)
    plt.xticks(alt_ticks, map(str, map(lambda alt:  (alt * altitudes.step + altitudes.start) // 1000, alt_ticks)))
    plt.yticks(inc_ticks, map(str, map(lambda inc: inc * inclinations.step + inclinations.start, inc_ticks)))
    plt.xlabel('Altitude (km)')
    plt.ylabel('Inclination (deg)')

    bar = plt.colorbar(label='Nodal precession rate ($rad \ s^{-1}$)', orientation='horizontal')
    bar.formatter.set_useMathText(True)

    clines = plt.contour(omegas, [tick for tick in bar.get_ticks() if omegas.min() <= tick <= omegas.max()], cmap=ListedColormap([[0, 0, 0]]))

    plt.clabel(clines, inline=True, manual=[])

    trans = plt.gca().transData
    cpoint = trans.transform(
        (alt_ticks[round(len(alt_ticks) * 0.15)],
        inc_ticks[round(len(inc_ticks) * 0.5)]))

    for index in range(len(clines.collections)):
        x, y = clines.find_nearest_contour(*cpoint, indices=(index,))[3:5]
        clines.add_label_near(x, y, inline=True, transform=False)

    plt.gcf().set_size_inches(6, 6)
    plt.savefig('graphs/nodal_precession.png', bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
