import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
# from matplotlib.font_manager import FontProperties
# from matplotlib.ticker import 
import numpy as np

from lib import nodal_precession
import constants

altitudes = range(400_000, 1600_000, 1_000)
alt_delta = 100_000
inclinations = range(0, 180, 1)
inc_delta = 15


omegas = np.zeros((len(inclinations), len(altitudes)))
for x, inc in enumerate(inclinations):
    for y, alt in enumerate(altitudes):
        omegas[x, y] = nodal_precession(inc, alt + constants.EARTH_MEAN_RADIUS)


alt_ticks = range(0, len(altitudes) + 1, alt_delta // 1000)
inc_ticks = range(0, len(inclinations) + 1, inc_delta)

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# aspect = len(altitudes)/len(inclinations)
im = plt.imshow(
    omegas,
    aspect='auto',
    extent=(alt_ticks.start, alt_ticks.stop - 1, inc_ticks.stop - 1, inc_ticks.start),
    cmap='seismic')
plt.gca().invert_yaxis()
# fig = plt.gcf()
# bbox = im.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
# print(bbox.get_points()[0][0] * fig.dpi, bbox.get_points()[0][1] * fig.dpi, sep='\t')

plt.title('Nodal Precession in Circular Orbits as a Function of Inclination and Altitude', fontsize=12, pad=15)
plt.xticks(alt_ticks, map(str, map(lambda alt:  (alt * altitudes.step + altitudes.start) // 1000, alt_ticks)))
plt.yticks(inc_ticks, map(str, map(lambda inc: inc * inclinations.step + inclinations.start, inc_ticks)))
plt.xlabel('Altitude (km)')
plt.ylabel('Inclination (deg)')

bar = plt.colorbar(label='Nodal precession rate ($rad \ s^{-1}$)', orientation='horizontal')
bar.formatter.set_useMathText(True)
# bar.ax.yaxis.label.set_font_properties(FontProperties(family='Computer Modern Roman'))
clines = plt.contour(omegas, [tick for tick in bar.get_ticks() if omegas.min() <= tick <= omegas.max()], cmap=ListedColormap([[0, 0, 0]]))

# bbox = im.clipbox #im.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
# print(bbox.get_points()[0][0] * fig.dpi, bbox.get_points()[0][1] * fig.dpi, sep='\t')
plt.clabel(clines, inline=True, manual=[])
# print(bbox.x0, bbox.x1, bbox.y0, bbox.y1)
trans = plt.gca().transData
# print()
cpoint = trans.transform((alt_ticks[round(len(alt_ticks) * 0.1)], inc_ticks[round(len(inc_ticks) * 0.5)]))#(bbox.x0 + (bbox.x1 - bbox.x0) * 0.5, bbox.y0 + (bbox.y1 - bbox.y0) * 0.5)
# print(cpoint)
# clines.add_label_near(*cpoint, inline=True)
# pts = []
for index in range(len(clines.collections)):
    x, y = clines.find_nearest_contour(*cpoint, indices=(index,))[3:5]
    # print(trans.inverted().transform((x, y)))
    # print(f'{x:0.2f} {y:0.2f}')
    clines.add_label_near(x, y, inline=True, transform=False)
    # pts.append((x, y))
# print(cpoint, (bbox.bounds[2] * fig.dpi, bbox.bounds[3] * fig.dpi / aspect), sep='\n')
# for label in plt.clabel(clines, inline=True, manual=[cpoint]):
#     label.set_rotation(0)
# print(bbox.x0, bbox.x1, bbox.y0, bbox.y1)

# for pt in pts:
#     p = plt.gca().add_patch(Rectangle((0, 0), 10, 10))

plt.gcf().set_size_inches(6, 6)
plt.savefig('nodal_precession.png', bbox_inches='tight')
# plt.show()
