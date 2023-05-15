#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: PCA.py image_file 


import numpy as np
import os, sys
import matplotlib.pyplot as plt
from astropy.io import fits


image_file = sys.argv[1]

with fits.open(image_file) as hdul:
    img = hdul[0].data

U,S,V = np.linalg.svd(img)

##### Varience explained for PC's of given image, used for verifying

#var_explained = np.round(S**2/np.sum(S**2), decimals=3)

#plt.bar(list(range(1,21)),var_explained[0:20], color='b')
#plt.xlabel('Singular Vector')
#plt.ylabel('Varience Explained')
#plt.show()

###### Image reconstruction

reconst_image_5 = np.matrix(U[:,:5])*np.diag(S[:5])*np.matrix(V[:5,:])
reconst_image_10 = np.matrix(U[:,:10])*np.diag(S[:10])*np.matrix(V[:10,:]) 
reconst_image_30 = np.matrix(U[:,:30])*np.diag(S[:30])*np.matrix(V[:30,:])
reconst_image_50 = np.matrix(U[:,:50])*np.diag(S[:50])*np.matrix(V[:50,:])

fig,axs = plt.subplots(2,2, figsize=(10,10))
axs[0,0].imshow(reconst_image_5)
axs[0,0].set_title("Reconstructed Image: 5 SV's")
axs[0,1].imshow(reconst_image_10)
axs[0,1].set_title("Reconstructed Image: 10 SV's")
axs[1,0].imshow(reconst_image_30)
axs[1,0].set_title("Reconstructed Image: 30 SV's")
axs[1,1].imshow(reconst_image_50)
axs[1,1].set_title("Reconstructed Image: 50 SV's")
plt.tight_layout()
plt.show()




