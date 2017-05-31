import scipy
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

# geometryFilename = 'OCS_222_7fs_LT_geometries.txt'
# geometryFilename = 'OCS_222_30fs_LT_geometries.txt'
# geometryFilename = 'OCS_222_60fs_LT_geometries.txt'
geometryFilename = 'OCS_222_100fs_LT_geometries.txt'
# geometryFilename = 'OCS_222_200fs_LT_geometries.txt'

print('geometryFilename = {:s}'.format(geometryFilename))

geometries = np.loadtxt(geometryFilename)
nGeometriesAll = geometries.shape[0]

# Deleting geometries (rows) that are all zero
geometries = geometries[~np.all(geometries == 0.0, axis=1)]

nGeometries = geometries.shape[0]
print('g={:d}, g_zero={:d}, g_good={:d}'.format(nGeometriesAll, nGeometriesAll-nGeometries, nGeometries))

r12 = geometries[:, 0]
r23 = geometries[:, 1]
theta = geometries[:, 2]

amu = 1.66053886e-27  # [kg], 1 atomic mass unit
m1, m2, m3 = amu*np.array([15.9994, 12.0107, 32.065])
M = m1 + m2 + m3

x1 = r12 * np.cos(np.deg2rad(90 + theta/2))
y1 = r12 * np.sin(np.deg2rad(90 + theta/2))

x2 = np.zeros(nGeometries)
y2 = np.zeros(nGeometries)

x3 = r23 * np.cos(np.deg2rad(90 - theta/2))
y3 = r23 * np.sin(np.deg2rad(90 - theta/2))

# Calculate COM
xCOM = (m1*x1 + m2*x2 + m3*x3) / M
yCOM = (m1*y1 + m2*y2 + m3*y3) / M

# Shift all coordinates so that the origin coincides with the COM. And convert to angstroms.
x1 = 1e10*(x1 - xCOM)
x2 = 1e10*(x2 - xCOM)
x3 = 1e10*(x3 - xCOM)
y1 = 1e10*(y1 - yCOM)
y2 = 1e10*(y2 - yCOM)
y3 = 1e10*(y3 - yCOM)

# Calculate mean atomic positions.
xO_mean = np.mean(x1)
yO_mean = np.mean(y1)
xC_mean = np.mean(x2)
yC_mean = np.mean(y2)
xS_mean = np.mean(x3)
yS_mean = np.mean(y3)

x = np.concatenate((x1, x2, x3))
y = np.concatenate((y1, y2, y3))

d = {'xO': x1, 'yO': y1, 'xC': x2, 'yC': y2, 'xS': x3, 'yS': y3}
df = pd.DataFrame(data=d)

d2 = {'xAll': x, 'yAll': y}
df2 = pd.DataFrame(data=d2)

# g = sns.JointGrid(x="xAll", y="yAll", data=df2, space=0)

# g = g.plot_joint(sns.kdeplot, cmap="Blues_d")
# g = g.plot_marginals(sns.kdeplot, shade=True)

# g = sns.jointplot("xAll", "yAll", data=df2, kind="kde", space=0, color="g")

# g = g.plot_joint(plt.scatter, color="b", edgecolor="white")
# g = g.plot_marginals(sns.distplot, kde=False, color="g")

# g = sns.jointplot(x="xAll", y="yAll", data=df2, kind="kde", color="m")
# g.plot_joint(plt.scatter, c="w", s=30, linewidth=1, marker="+")
# g.ax_joint.collections[0].set_alpha(0)
# g.set_axis_labels("$X$", "$Y$")

sns.set(style="ticks", color_codes=True)

# Plotting geometries.
p = sns.JointGrid(x=df['xO'], y=df['yO'])
p.plot_joint(sns.kdeplot, cmap="Reds", shade=True, shade_lowest=False, bw='scott')
p.plot_joint(plt.scatter, s=2, color='#7D0112', alpha=0.2)

p.x = df['xC']
p.y = df['yC']
p.plot_joint(sns.kdeplot, cmap="Greys", shade=True, shade_lowest=False, bw='scott')
p.plot_joint(plt.scatter, s=2, color='#4B4B4B', alpha=0.2)

p.x = df['xS']
p.y = df['yS']
p.plot_joint(sns.kdeplot, cmap="YlOrBr", shade=True, shade_lowest=False, bw='scott')
p.plot_joint(plt.scatter, s=2, color='#E99A2C', alpha=0.2)

# Marginal distributions = kernel density estimate plots
sns.kdeplot(df['xO'], ax=p.ax_marg_x, vertical=False, color='#7D0112', shade=True)
sns.kdeplot(df['yO'], ax=p.ax_marg_y, vertical=True, color='#7D0112', shade=True)
sns.kdeplot(df['xC'], ax=p.ax_marg_x, vertical=False, color='#4B4B4B', shade=True)
sns.kdeplot(df['yC'], ax=p.ax_marg_y, vertical=True, color='#4B4B4B', shade=True)
sns.kdeplot(df['xS'], ax=p.ax_marg_x, vertical=False, color='#E99A2C', shade=True)
sns.kdeplot(df['yS'], ax=p.ax_marg_y, vertical=True, color='#E99A2C', shade=True)

# Calculates medium from kde but I wanted center of 2D KDE and this works for 1D. Doing x,y doesn't give the same
# center as the 2D KDE. Will just find them by hand.
# x, y = p.ax_marg_x.get_lines()[0].get_data()
# cdf = scipy.integrate.cumtrapz(y, x, initial=0)
# nearest_05 = np.abs(cdf-0.5).argmin()
# xO_median = x[nearest_05]

# Marginal distributions = histograms
# p.ax_marg_x.hist(df['xO'], alpha=0.5, color='#7D0112')
# p.ax_marg_y.hist(df['yO'], orientation='horizontal', color='#7D0112', alpha=0.5)
# p.ax_marg_x.hist(df['xC'], color='#4B4B4B', alpha=0.5, range=(np.min(df['xC']), np.max(df['xC'])))
# p.ax_marg_y.hist(df['yC'], orientation='horizontal', color='#4B4B4B', alpha=0.5,
#                  range=(np.min(df['yC']), np.max(df['yC'])))
# p.ax_marg_x.hist(df['xS'], color='#E9B62D', alpha=0.5, range=(np.min(df['xS']), np.max(df['xS'])))
# p.ax_marg_y.hist(df['yS'], orientation='horizontal', color='#E9B62D', alpha=0.5,
#                  range=(np.min(df['yS']), np.max(df['yS'])))

o_patch = mpatches.Patch(color='#7D0112', label='oxygen')
c_patch = mpatches.Patch(color='#4B4B4B', label='carbon')
s_patch = mpatches.Patch(color='#E99A2C', label='sulfur')
plt.legend(handles=[o_patch, c_patch, s_patch])

# Plot "average geometry" using centroid positions.
# p.x = [xO_mean, xC_mean, xS_mean]
# p.y = [yO_mean, yC_mean, yS_mean]
# p.plot_joint(plt.scatter, marker='x', color='black', s=50)
# plt.plot([xO_mean, xC_mean], [yO_mean, yC_mean], linewidth=1, color='black')
# plt.plot([xC_mean, xS_mean], [yC_mean, yS_mean], linewidth=1, color='black')

# 7 fs
# xO_modal, yO_modal = -2.121, 0.019
# xC_modal, yC_modal = -0.384, -0.091
# xS_modal, yS_modal = 1.188, 0.010

# 30 fs
# xO_modal, yO_modal = -2.594, 0.022
# xC_modal, yC_modal = -0.2727, -0.127
# xS_modal, yS_modal = 1.459, 0.0113

# 60 fs
# xO_modal, yO_modal = -2.814, 0.0453
# xC_modal, yC_modal = -0.3958, -0.1304
# xS_modal, yS_modal = 1.545, 0.0196

# 100 fs
xO_modal, yO_modal = -3.007, 0.0349
xC_modal, yC_modal = -0.499, -0.1671
xS_modal, yS_modal = 1.651, 0.0120

# 200 fs
# xO_modal, yO_modal = -3.479, 0.0635
# xC_modal, yC_modal = -0.3239, -0.1618
# xS_modal, yS_modal = 1.908, 0.0145

# Plot "modal geometry" or "most likely geometry" at the peaks of the bivariate KDE.
# Coordinates were taken using the built-in plot viewer. Calculating them would have been too much work (see above).
p.x = [xO_modal, xC_modal, xS_modal]
p.y = [yO_modal, yC_modal, yS_modal]
# p.plot_joint(plt.scatter, marker='o', color='black', s=50)
plt.plot([xO_modal, xC_modal], [yO_modal, yC_modal], linewidth=1, color='black', alpha=0.5)
plt.plot([xC_modal, xS_modal], [yC_modal, yS_modal], linewidth=1, color='black', alpha=0.5)

plt.xlabel('$x$ (angstroms)')
plt.ylabel('$y$ (angstroms)')
p.ax_marg_x.legend_.remove()
p.ax_marg_y.legend_.remove()
# plt.xlim([-4.1, 2.2])
# plt.ylim([-2, 2])

# Calculate and print average and modal geometries.
O = np.array([xO_mean, yO_mean])
C = np.array([xC_mean, yC_mean])
S = np.array([xS_mean, yS_mean])

CO = O - C
CS = S - C
rCO = np.sqrt(np.linalg.norm(CO))
rCS = np.sqrt(np.linalg.norm(CS))

theta = np.rad2deg(np.arccos(np.dot(CO, CS)/rCO**2/rCS**2))

print('Average geometry:')
print('rCO = {:f} A, rCS = {:f} A, theta = {:f} deg'.format(rCO, rCS, theta))

O = np.array([xO_modal, yO_modal])
C = np.array([xC_modal, yC_modal])
S = np.array([xS_modal, yS_modal])

CO = O - C
CS = S - C
rCO = np.sqrt(np.linalg.norm(CO))
rCS = np.sqrt(np.linalg.norm(CS))

theta = np.rad2deg(np.arccos(np.dot(CO, CS)/rCO**2/rCS**2))

print('Modal or "most likely" geometry:')
print('rCO = {:f} A, rCS = {:f} A, theta = {:f} deg'.format(rCO, rCS, theta))

plt.show()
