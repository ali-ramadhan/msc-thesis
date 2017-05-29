import scipy
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

geometryFilename = 'OCS_222_7fs_LT_geometries.txt'  # str(sys.argv[1])

geometries = np.loadtxt(geometryFilename)
nGeometriesAll = geometries.shape[0]

# Deleting geometries (rows) that are all zero
geometries = geometries[~np.all(geometries == 0.0, axis=1)]

nGeometries = geometries.shape[0]
print('g={:d}, g_zero={:d}'.format(nGeometriesAll, nGeometriesAll-nGeometries))

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

p.x = [xO_mean, xC_mean, xS_mean]
p.y = [yO_mean, yC_mean, yS_mean]
p.plot_joint(plt.scatter, marker='x', color='black', s=50)
plt.plot([xO_mean, xC_mean], [yO_mean, yC_mean], linewidth=1, color='black')
plt.plot([xC_mean, xS_mean], [yC_mean, yS_mean], linewidth=1, color='black')

plt.xlabel('$x$ (angstroms)')
plt.ylabel('$y$ (angstroms)')
p.ax_marg_x.legend_.remove()
p.ax_marg_y.legend_.remove()
# plt.xlim([-4.1, 2.2])
# plt.ylim([-2, 2])

plt.show()
