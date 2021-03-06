from field import *
from interpolacja import *
import numpy as np
einterp = EdgeField.interp

adj = 0
DlPrzX = 1.
DlPrzY = 1.
n = 30
m = n
dx = DlPrzX/n
dy = DlPrzY/n
x0, y0, dl = (0, 0, 0)

diffusivity = 0.0001

node_c, cells, bound = siatka_regularna_prost(n, m, dx, dy, x0, y0)
mesh = Mesh(node_c, cells, bound)

coeffFieldX = SurfField(mesh, Neuman)
coeffFieldX.data[:] = 1
coeffFieldX.updateBoundaryValues()

coeffFieldY = SurfField(mesh, Neuman)
coeffFieldY.data[:] = 1
coeffFieldY.updateBoundaryValues()

coeff_diff_Edge = EdgeField.vector(EdgeField.interp(coeffFieldX), EdgeField.interp(coeffFieldY))

T = SurfField(mesh, Dirichlet)

T.setBoundaryCondition(Dirichlet(mesh, 2, 1))
# T.setBoundaryCondition(Dirichlet(mesh, 0, -100))
# T.setBoundaryCondition(Neuman(mesh, 1, 0))
# T.setBoundaryCondition(Neuman(mesh, 3, 0))

Ux, Uy = generate_u(mesh, constant_velocity)             # generuje 2 pola pr w sr komorek

Ux.setBoundaryCondition(Dirichlet(mesh, 0, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 1, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 2, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 3, 0))

edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
phi = edgeU.dot(mesh.normals)

# phi = EdgeField(mesh)
# phi.data = generate_phi_r(mesh, constant_velocity) #np.multiply(, mesh.eLengths)

adj = adjustPhi(phi)
# Md, Fd = laplace(diffusivity, T)

Md, Fd = laplace(coeff_diff_Edge, T)
Mc, Fc = div(phi, T)

# pkt1 = n/2 + n*n/2                           # pkt srodek
pkt2 = n/2 + n*5                              # srodek 4 wiersze od spodu
# pkt3 = n/2 + n*(n-5)                         # srodek 4 wiersze od gory
# Fd[pkt1] += -300
# Fd[pkt2] += 0.01
# Fd[pkt3] += -200

# M = Mc - Md  #
# F = Fd - Fc  #

M = - Md  #
F = Fd   #


from scipy.sparse.linalg.isolve.iterative import bicgstab
res, info = bicgstab(A=M.sparse, b=F, x0=T.data, tol=1e-8, maxiter=250e3)

T.setValues(res)

print ">>>> info: ", info

# animate_contour_plot([T.data.reshape(n, m)], skip=10, repeat=False, interval=75)
animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1, m+1))], skip=10, nLevels=16, repeat=False, interval=75, nN=n, diff=diffusivity, adj=adj)

# magU = np.sqrt(Ux.data**2 + Uy.data**2)
# animate_contour_plot([magU.reshape(n, m)], nLevels=16, skip=10, repeat=False, interval=75, dataRange=[0, 1], nN=n, diff=diffusivity, adj=adj)


plt.show()
