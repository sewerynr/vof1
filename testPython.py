from fvMatrix import fvMatrix
import numpy as np

# tab = list(range(10))
#
#
# tab = [2*t for t in tab]

# print tab
#
# import numpy as np
#
#
# data = np.array([i for i in range(10)])
#
# indices = [[1,2],[0, 1]]
#
# print data[indices]
#
# macierz = fvMatrix.diagonal(4)
#
#
# macierz[1, 1] = 3
#
# macierz[2, 2] = 4
#
# macierz[1, 3] += 3
#
# macierz[0, 1] += 3
# macierz[0, 1] += 3
# macierz[0, 0] += 3
# macierz[3, 3] -= 1
#
# U = [1,3,5,2]
#
# print macierz.offdiagmul(U)
# print macierz.sparse.todense()

import numpy as np

index = []
for a in range(10):
    index.append(a)


print index