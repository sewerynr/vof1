from field import *
from interpolacja import *
import numpy as np

DlPrzX = 1.
DlPrzY = 1.
n = 50
dx = DlPrzX/n
dy = DlPrzY/n
x0, y0, dl = (0, 0, 0)
diffusivity = 0.001

node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)
mesh = Mesh(node_c, cells, bound)


T = SurfField(mesh, Dirichlet)
T.setBoundaryCondition(Dirichlet(mesh, 2, 1))

np.set_printoptions(precision=3)

einterp = EdgeField.interp

Ux, Uy = generate_u(mesh, quadratic_velocity)
edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
phi = edgeU.dot(mesh.normals)
# phi = EdgeField(mesh)
# phi.data = np.multiply(generate_phi_r(mesh, quadratic_velocity), mesh.eLengths)

Md, Fd = laplace(diffusivity, T)
Mc, Fc = div(phi, T)

from scipy.sparse.linalg.isolve.iterative import bicgstab

M = Mc - Md
F = Fc - Fd
print M.sparse.shape

res, info = bicgstab(A=M.sparse, b=F, x0=T.data, maxiter=250e3)
T.setValues(res)

print ">>>> info: ", info

# animate_contour_plot([T.data.reshape(n, n)], skip=10, repeat=False, interval=75)
animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1,n+1))], skip=10, repeat=False, interval=75)

# magU = np.sqrt(Ux.data**2 + Uy.data**2)
# print magU.shape, n*n, Ux.data.shape
# animate_contour_plot([magU.reshape(n, n)], skip=10, repeat=False, interval=75)

plt.show()
