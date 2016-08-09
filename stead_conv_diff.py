from field import *
from interpolacja import *
import numpy as np
einterp = EdgeField.interp

adj = 0
DlPrzX = 1.
DlPrzY = 1.
n = 10
m = n
dx = DlPrzX/n
dy = DlPrzY/n
x0, y0, dl = (0, 0, 0)
diffusivity = 0.1

node_c, cells, bound = siatka_regularna_prost(n, m, dx, dy, x0, y0)
mesh = Mesh(node_c, cells, bound)

T = SurfField(mesh, Dirichlet)

T.setBoundaryCondition(Dirichlet(mesh, 2, 1))
# T.setBoundaryCondition(Dirichlet(mesh, 0, -100))
# T.setBoundaryCondition(Neuman(mesh, 1, 0))
# T.setBoundaryCondition(Neuman(mesh, 3, 0))


Ux, Uy = generate_u(mesh, constant_velocity)    # generuje 2 pola pr w sr komorek

Ux.setBoundaryCondition(Dirichlet(mesh, 0, 1))
Ux.setBoundaryCondition(Dirichlet(mesh, 1, 1))
Ux.setBoundaryCondition(Dirichlet(mesh, 2, 1))
Ux.setBoundaryCondition(Dirichlet(mesh, 3, 1))

edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
# print "einterpUx",einterp(Ux).data
# print mesh.normals
phi = edgeU.dot(mesh.normals)


# phi = EdgeField(mesh)
# phi.data = generate_phi_r(mesh, constant_velocity) #np.multiply(, mesh.eLengths)

adj = adjustPhi(phi)

Md, Fd = laplace(diffusivity, T)
Mc, Fc = div(phi, T)

from scipy.sparse.linalg.isolve.iterative import bicgstab

# pkt1 = n/2 + n*n/2                           # pkt srodek
pkt2 = n/2 + n*10                             # srodek 4 wiersze od spodu
# pkt3 = n/2 + n*(n-5)                         # srodek 4 wiersze od gory
# Fd[pkt1] += -300
Fd[pkt2] += 0.00001
# Fd[pkt3] += -200

M = Mc - Md  #
F = Fd - Fc  #

# print "Mc diag", Mc.diag
# print "Mc data",Mc.data
#
# print "Md data",Md.data
# print "Md diag",Md.diag

# print "Fc",Fc
# print "Fd",Fd
# print "F",F

# print "M data",M.data
# print "M diag",M.diag
# print M.sparse


print Mc.sparse
print F


res, info = bicgstab(A=M.sparse, b=F, x0=T.data, tol=1e-8, maxiter=250e3)

T.setValues(res)

print ">>>> info: ", info

# animate_contour_plot([T.data.reshape(n, m)], skip=10, repeat=False, interval=75)
animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1, m+1))], skip=10, nLevels=16, repeat=False, interval=75, nN=n, diff=diffusivity, adj=adj)

magU = np.sqrt(Ux.data**2 + Uy.data**2)
print magU.shape, n*m, Ux.data.shape
animate_contour_plot([magU.reshape(n, m)], skip=10, repeat=False, interval=75)

plt.show()
