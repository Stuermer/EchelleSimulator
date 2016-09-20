import numpy as np
import matplotlib.pyplot as plt
import scipy.misc

# data = np.genfromtxt('../o90.dat', delimiter=',')
data = np.genfromtxt('../image2.dat', delimiter=',')
# data2 = np.genfromtxt('test.dat', delimiter=',')

# data = scipy.misc.imread('test.png')
plt.figure()
# plt.plot(data)
plt.imshow(data, interpolation='None')

plt.show()

