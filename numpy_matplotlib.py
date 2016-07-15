from numpy.distutils.system_info import numarray_info

import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

a = np.array([i for i in range(9)])
b = np.ones((3, 3))

print a
a = a.reshape((3, 3))

print a

a = a.T
print a


print a[0, :]
print b

b = 3 * b.reshape(9) + 1

m = np.matrix(a)
print m

# print solve(m, np.ones((3,1)))

x = np.linspace(0, 1, 30)
y = np.sin(x)


plt.figure()
plt.plot(x, y, 'g.', markersize=20)
plt.plot(y, x)
plt.show()
