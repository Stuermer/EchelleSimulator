import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import tables
from astropy.io import fits

h5file = tables.open_file('../MaroonX.hdf')
# h5file2 = tables.open_file('../FIES2.hdf')
data = h5file.root.image
# data2 = h5file2.root.image
# data = np.array(data)+np.array(data2)

plt.figure()
plt.imshow(data, interpolation='None', vmin=0.)

plt.figure()
plt.plot(data[:,1024])
plt.show()

#
# hdu = fits.PrimaryHDU()
# hdu.data = data
# hdu.writeto('FIES_format.fits')

# d = np.genfromtxt('../transform.csv', delimiter=";")
#
# f, axarr = plt.subplots(3, 2, sharex=True)
#
# # order = np.atleast_1d(d[:,0])
# # if np.ptp(order)>1:
# #     scaled_z = np.array((order - np.min(order)), dtype=float) / np.ptp(order)
# wl = d[:,1]
# # colors = plt.cm.rainbow(scaled_z)
# # for o,c in zip(order, colors):
# c='b'
# axarr[0, 0].scatter(wl, d[:,6], c=c)
# axarr[0, 0].set_title('translation x')
# axarr[1, 0].scatter(wl, d[:,7], c=c)
# axarr[1, 0].set_title('translation y')
# axarr[0, 1].scatter(wl, d[:,2], c=c)
# axarr[0, 1].set_title('scale x')
# axarr[1, 1].scatter(wl, d[:,3], c=c)
# axarr[1, 1].set_title('scale y')
# axarr[2, 0].scatter(wl, d[:,5], c=c)
# axarr[2, 0].set_title('rotation')
# axarr[2, 1].scatter(wl, d[:,4], c=c)
# axarr[2, 1].set_title('shear')
# plt.show()