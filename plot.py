from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mplcolors



data = np.loadtxt('bands.out')

k = data[:,0]

plt.figure(figsize=[12.,8.])
plt.title(r'$E(k)$')
plt.grid(True)

cNorm  = mplcolors.Normalize(vmin=0, vmax=data.shape[1])
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.jet)


for ii in range(data.shape[1]-1):
    color = scalarMap.to_rgba(ii+1)
    plt.plot(k,data[:,ii+1],color=color,label='n={}'.format(ii+1))
plt.legend()
plt.show()