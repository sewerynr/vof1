from field import *
from interpolacja import *
from fvMatrix import fvMatrix
from fvm1 import *

DlPrzX = 1.
DlPrzY = 1.

n = 50
dx = DlPrzX/n
dy = DlPrzY/n

dt = 0.001                                                 # CFL u*dt/dx <= 1   => dt = dx/u
tp = 0
tk = 20

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)


node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                         # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

diffusivity = 0.001

Tdir = 1
TdirWB = 0

T = SurfField(mesh, Dirichlet)                  # temp w srodkach komorek


#T.setBoundaryCondition(Neuman(mesh, 0, 0))               # zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)
#T.setBoundaryCondition(Dirichlet(mesh, 0, TdirWB))

#T.setBoundaryCondition(Neuman(mesh, 1, 0))
#T.setBoundaryCondition(Dirichlet(mesh, 1, TdirWB))

#T.setBoundaryCondition(Neuman(mesh, 2, 0))
T.setBoundaryCondition(Dirichlet(mesh, 2, 10))

np.set_printoptions(precision=3)

einterp = EdgeField.interp

Ux, Uy = generate_u(mesh, quadratic_velocity)                   # stworz w srodkach komorek
edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))              # przenies za pomoca .interp do krawedzi
phi = edgeU.dot(mesh.normals)                                   # oblicz UnaKrawedzi * normalnaDoKrawedzi = (skalar)


Md, Fd = laplace(diffusivity, T)
Mc, Fc = div(phi, T)

# pkt1 = n/2 + n*n/2                       # pkt srodek
# pkt2 = n/2 + n*5                         # srodek 4 wiersze od spodu
# pkt3 = n/2 + n*(n-5)                     # srodek 4 wiersze od gory
# Fd[pkt1] += -300
# Fd[pkt2] += -1
# Fd[pkt3] += -200

Mass = fvMatrix.diagonal(mesh, mesh.cells_areas / dt)

M = Mass + Mc - Md
F = Fc - Fd

# wypelnij polowe temperatura
# for i, point in enumerate(T.data):
#     if i < (n**2)/2:
#         T.data[i] = 10

Results = list()
Tn = T.data.reshape((n, n))
Results.append(Tn)

from scipy.sparse.linalg.isolve.iterative import bicgstab

licznik = 0
#zliczacz = 0

for iter in range(int(nt)):
    licznik = licznik + 1
    #print 'time iteration:', iter
    solution = bicgstab(A=M.sparse, b=F + source(mesh, T.data/dt), x0=T.data)[0]
    T.setValues(solution)
    if licznik == 100:
        Results.append(T.data.reshape((n, n)))
        #zliczacz = zliczacz + 1
        licznik = 0
    print "pozostalo: ", int(nt - iter)
#print zliczacz

# Animate results:
anim = animate_contour_plot(Results, skip=10, repeat=False, interval=75, nLevels=20, dataRange=[0., 10])
plt.show()


# from interpolacja import inter
# from matplotlib.pyplot import quiver
#
# animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1, n+1))], skip=10, repeat=False, interval=75, dataRange=[0, 10])
#
# q = quiver(mesh.cell_centers[:, 0], mesh.cell_centers[:, 1], Ux[:], Uy[:])
#
# plt.show()


# magU = np.sqrt(Ux.data**2 + Uy.data**2)
# print magU.shape, n*n, Ux.data.shape
# animate_contour_plot([magU.reshape(n, n)], skip=10, repeat=False, interval=75)
# plt.show()

