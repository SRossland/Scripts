import os, sys, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


with fits.open('grid6x6A.fits') as hdul:
    Agrid = hdul[0].data

with fits.open('grid6x6B.fits') as hdul:
    Bgrid = hdul[0].data

with fits.open('fullmaskA_final.fits') as hdul:
    Am = hdul[0].data

with fits.open('fullmaskB_final.fits') as hdul:
    Bm = hdul[0].data

im = plt.imshow(Agrid*Am, interpolation='none', origin='lower')
values = np.ravel(Agrid.ravel())

values = list(set(values))

values = [int(i) for i in values]

vals = values[1:]
plt.set_cmap('jet')
colors = [im.cmap(im.norm(value)) for value in vals]

patches = [mpatches.Patch(color=colors[i], label="Pos. {l}".format(l=vals[i])) for i in range(len(vals))]

plt.legend(handles=patches, bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
plt.grid(True)
plt.savefig('GridLegendA.png')

plt.close()
im = plt.imshow(Bgrid*Bm, interpolation='none', origin = 'lower')
plt.set_cmap('jet')
color = [im.cmap(im.norm(value)) for value in vals]
patches = [mpatches.Patch(color=colors[i], label='Pos. {l}'.format(l=vals[i])) for i in range(len(vals))]

plt.legend(handles=patches, bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
plt.grid(True)
plt.savefig('GridLegendB.png')








