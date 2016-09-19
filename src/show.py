import numpy as np
import matplotlib.pyplot as plt
import scipy.misc

data = np.genfromtxt('../image2.dat', delimiter=',', skip_footer=3000)
# data = np.genfromtxt('image2.dat', delimiter=',')
# data2 = np.genfromtxt('test.dat', delimiter=',')

# data = scipy.misc.imread('test.png')
plt.figure()
# plt.plot(data2)
plt.imshow(data, interpolation='None')

plt.show()

