from field import *
from interpolacja import *
import numpy as np
einterp = EdgeField.interp

adj = 0
DlPrzX = 1.
DlPrzY = 1.

n = 1           # wiersze
m = 20          # kolumny

dx = DlPrzX/m
dy = DlPrzY/n
x0, y0, dl = (0, 0, 0)
diffusivity = 1

node_c, cells, bound = siatka_regularna_prost(n, m, dx, dy, x0, y0)
mesh = Mesh(node_c, cells, bound)

# RYSUJ MESH
# draw_edges(mesh.xy, mesh.list_kr)
# plt.show()

T = SurfField(mesh, Dirichlet)

# T.setBoundaryCondition(Dirichlet(mesh, 2, 1))
# T.setBoundaryCondition(Dirichlet(mesh, 0, -100))
# T.setBoundaryCondition(Neuman(mesh, 1, 0))
# T.setBoundaryCondition(Neuman(mesh, 3, 0))

T.setBoundaryCondition(Neuman(mesh, 0, 0))
T.setBoundaryCondition(Neuman(mesh, 2, 0))


Ux, Uy = generate_u(mesh, constant_velocity)    # generuje 2 pola pr w sr komorek

Ux.setBoundaryCondition(Dirichlet(mesh, 0, 10))
Ux.setBoundaryCondition(Dirichlet(mesh, 1, 10))
Ux.setBoundaryCondition(Dirichlet(mesh, 2, 10))
Ux.setBoundaryCondition(Dirichlet(mesh, 3, 10))

edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
# print "einterpUx",einterp(Ux).data
# print mesh.normals
phi = edgeU.dot(mesh.normals)


# phi = EdgeField(mesh)
# phi.data = generate_phi_r(mesh, constant_velocity) #np.multiply(, mesh.eLengths)

# adj = adjustPhi(phi, mesh)

Md, Fd = laplace(diffusivity, T)
Mc, Fc = div(phi, T)

from scipy.sparse.linalg.isolve.iterative import bicgstab

# pkt1 = n/2 + n*n/2                           # pkt srodek
# pkt2 = n/2 + n*10                             # srodek 4 wiersze od spodu
# pkt3 = n/2 + n*(n-5)                         # srodek 4 wiersze od gory
# Fd[pkt1] += -300
# Fd[pkt2] += 0.00001
# Fd[pkt3] += -200

M = Mc - Md  #
F = Fd - Fc  #

from math import *

def Ffun(x, y, ux, uy):
    return ux * cos(x)*cos(y) - uy * sin(x)*sin(y) + 2*sin(x)*cos(y)

def Ffun1(x, y, ux, uy):
    return ux * cos(x) + sin(x)

def Ffun2(x, y, ux, uy):
    return ux*2*x - 2

for i, celCXY in enumerate(mesh.cell_centers):
    F[i] += (Ffun(celCXY[0], celCXY[1], 1, 0)) * mesh.cells_areas[i]

res, info = bicgstab(A=M.sparse, b=F, x0=T.data, tol=1e-8, maxiter=250e3)

T.setValues(res)
# print "T.dat", T.data
print ">>>> info: ", info

animate_contour_plot([T.data.reshape(n, m)], skip=10, repeat=False, interval=75)   # nie dziala
animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1, m+1))], skip=10, nLevels=16, repeat=False, interval=75, nN=n, diff=diffusivity, adj=adj)

# magU = np.sqrt(Ux.data**2 + Uy.data**2)
# print magU.shape, n*n, Ux.data.shape
# animate_contour_plot([magU.reshape(n, m)], skip=10, repeat=False, interval=75)


# draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
# draw_edges(mesh.xy, mesh.list_kr)


plt.show()
