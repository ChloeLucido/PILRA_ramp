# -*- coding: utf-8 -*-
"""Figure2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Takwim6IDOUFHBlr_4o_8hkKJkmA8kPQ
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read in the line in the file that is pertinent
pilra_speeds = pd.read_csv('pilra_speeds.txt', skiprows=[0,2,3,4,5,6,7],
                           dtype=np.float64, header=None)

pilra_mutant_speed = pd.read_csv('pilra_mutant_speeds.txt', skiprows=[0],
                                dtype=np.float64, header=None)

# Calculate sliding window means 

speeds = pilra_speeds.values[0]
N = 9 # Use a 9-codon window (the default setting for ExtRamp)
windowMeans = pd.Series(speeds).rolling(window=N).mean().iloc[N-1:].values

mutant_speeds = pilra_mutant_speed.values[0]
N = 9
mutant_windowMeans = pd.Series(mutant_speeds).rolling(window=N).mean().iloc[N-1:].values

positions = range(1, 295) # for x-axis

# Plot the data

fig, ax = plt.subplots() # create a new figure with a default 111 subplot

ax.plot(positions, mutant_windowMeans, 'r', label="rs2405442")
ax.plot(positions, windowMeans, 'b', label="NP_038467.2")
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 2.5, loc=8)
axins.plot(positions, mutant_windowMeans, 'r')
axins.plot(positions, windowMeans, 'b')
ax.set_xlim([0,301]) 
ax.set_ylim([0,1.1])
x1, x2, y1, y2 = 0, 30, 0.7, 0.95 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2)
plt.yticks(visible=False)
plt.xticks(visible=False)
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")

ax.set_ylabel('Relative Codon Adaptiveness') 
ax.set_xlabel('Codon Position')
ax.set_title('Relative Codon Adaptiveness for PILRA') 
ax.annotate('Ramp sequence region', xy=(15,.6),
             xytext=(7, .3), arrowprops=dict(facecolor='black', shrink=0.05))
ax.legend(loc="lower right")


plt.savefig('pilra_speeds.png',dpi=300)
