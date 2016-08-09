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

# index = []
# for a in range(10):
#     index.append(a)
#
#
# print index

m = 4
n = 1

boundaries_BC = list()
mat = list()
for i in range (m):
    mat.append([0, 0])

boundaries_BC.append(mat)

mat = list()
for i in range (n):
    mat.append([0, 0])

boundaries_BC.append(mat)

mat = list()
for i in range (m):
    mat.append([0, 0])

boundaries_BC.append(mat)

mat = list()
for i in range (n):
    mat.append([0, 0])

boundaries_BC.append(mat)

# boundaries_BC.append([[0, 0]] * n)
# boundaries_BC.append([[0, 0]] * m)
# boundaries_BC.append([[0, 0]] * n)

np.array(boundaries_BC)

print boundaries_BC
print boundaries_BC[0]
print boundaries_BC[0][0]
boundaries_BC[0][0][0] = 1
print boundaries_BC[0][0][0]

boundaries_BC[2][2][1] = 100
print boundaries_BC

# print  boundaries_BC
