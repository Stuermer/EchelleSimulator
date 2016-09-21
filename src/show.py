import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import tables

# data = np.genfromtxt('../o90.dat', delimiter=',')
h5file = tables.open_file('../image2.hdf')
data = h5file.root.image
# data = np.genfromtxt('../image2.dat', delimiter=',')
# data2 = np.genfromtxt('test.dat', delimiter=',')

# data = scipy.misc.imread('test.png')
plt.figure()
# plt.plot(data)
plt.imshow(data, interpolation='None')

plt.show()

