from field import *
from interpolacja import *
from fvMatrix import fvMatrix
from fvm1 import *
einterp = EdgeField.interp


adj = 0

DlPrzX = 1.
DlPrzY = 1.

n = 20
m = n
dx = DlPrzX/m
dy = DlPrzY/n

dt = 0.004                                                      # CFL u*dt/dx <= 1   => dt = dx/u
tp = 0
tk = 20

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)
diffusivity = 0.001

node_c, cells, bound = siatka_regularna_prost(n, m,  dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                                # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

# Tdir = 1
# TdirWB = 0.

T = SurfField(mesh, Dirichlet)                                    # temp w srodkach komorek


T.setBoundaryCondition(Neuman(mesh, 0, 0))                       # zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)
# T.setBoundaryCondition(Dirichlet(mesh, 0, TdirWB))

T.setBoundaryCondition(Neuman(mesh, 1, 0))
# T.setBoundaryCondition(Dirichlet(mesh, 1, TdirWB))

T.setBoundaryCondition(Neuman(mesh, 2, 10))
# T.setBoundaryCondition(Dirichlet(mesh, 2, 1))

T.setBoundaryCondition(Neuman(mesh, 3, 0))
# T.setBoundaryCondition(Dirichlet(mesh, 3, 1.))

T.data[:] = 0.

np.set_printoptions(precision=6)

Ux, Uy = generate_u(mesh, constant_velocity)                   # stworz w srodkach komorek

Ux.setBoundaryCondition(Dirichlet(mesh, 3, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 1, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 0, 0))
Ux.setBoundaryCondition(Dirichlet(mesh, 2, 0))


edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))              # przenies za pomoca .interp do krawedzi
phi = edgeU.dot(mesh.normals)                                   # oblicz UnaKrawedzi * normalnaDoKrawedzi = (skalar)


# phi = EdgeField(mesh)
# phi.data = np.multiply(generate_phi_1(mesh), mesh.eLengths)

 #!!!!!!!!!!!!!!!!!    blad pola predkosci - zrodlowosc      !!!!!!!!!!!!!!!!!
# animate_contour_plot([inter(mesh.xy, mesh.cells, eInt(phi)).reshape((n+1, n+1))], skip=10, repeat=False, interval=75)

adj = adjustPhi(phi)               #pobiera i pracuje i zwraca phi zmodyfikowne

# animate_contour_plot([inter(mesh.xy, mesh.cells, eInt(phi)).reshape((n+1, n+1))], skip=10, repeat=False, interval=75)
# plt.show()


Md, Fd = laplace(diffusivity, T)
Mc, Fc = div(phi, T)

# pkt1 = n/2 + n*n/2                       # pkt srodek
# pkt2 = (n-1)/2 + n*6                     # srodek 4 wiersze od spodu
# pkt3 = n/2 + n*(n-5)                     # srodek 4 wiersze od gory
# Fd[pkt1] += 0.1
# Fd[pkt2] += 0.02
# Fd[pkt3] += -200

Mass = fvMatrix.diagonal(mesh, mesh.cells_areas / dt)

M = Mass + Mc - Md
F = Fd - Fc

# wypelnij polowe temperatura
# for i, point in enumerate(T.data):
#     if i < (n**2)/2:
#         T.data[i] = 1

Results = list()
Tn = T.data.reshape((n, n))
Results.append(Tn)

from scipy.sparse.linalg.isolve.iterative import bicgstab

licznik = 0
step = 10


for iter in range(int(nt)):
    licznik = licznik + 1

    solution = bicgstab(A=M.sparse, b=F + source(mesh, T.data/dt), x0=T.data, tol=1e-8)[0]

    # solution = np.linalg.solve(a=M.dens, b=F + source(mesh, T.data / dt))

    T.setValues(solution)
    if licznik == step:
        Results.append(T.data.reshape((n, n)))
        licznik = 0
        step += 10 + 1
    print "iter left: ", int(nt - iter)

anim = animate_contour_plot(Results, skip=2, repeat=False, interval=75, nLevels=20, nN=n, dataRange=[0., 2], diff=diffusivity, dt=dt, adj=adj)


# from interpolacja import inter
# from matplotlib.pyplot import quiver
# q = quiver(mesh.cell_centers[:, 0], mesh.cell_centers[:, 1], Ux[:], Uy[:])

# animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1, n+1))], skip=10, repeat=False, interval=75, dataRange=[0, 10])

# magU = np.sqrt(Ux.data**2 + Uy.data**2)
# print magU.shape, n*n, Ux.data.shape
# animate_contour_plot([magU.reshape(n, n)], skip=10, repeat=False, interval=75)
plt.show()