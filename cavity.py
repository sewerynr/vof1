from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Zmienne w czasie zapis macierzy zadkich jako "wektory" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DlPrzX = 1.; DlPrzY = 1.

n = 50                                                    # ilosc podzialow

dx = DlPrzX/n
dy = DlPrzY/n

dt = 0.001                                                 # CFL u*dt/dx <= 1
tp = 0
tk = 0.1

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)


import time

# wsp_wezl, cells, bounduary = siatka_regularna_prost(n, dx, dy, x0, y0)  czyli metoda siatka_regularna_prost zwroci to co Mesh potrzebuje czyli nodes, cells, boundary

node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                         # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...


Ux = SurfField(mesh)                                       # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))
Uy = SurfField(mesh)
p = SurfField(mesh)

#print mesh.cell_centers


#T.setBoundaryCondition(Neuman(mesh, 0, 0))               # zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)
Ux.setBoundaryCondition(Dirichlet(mesh, 0, 0))
Uy.setBoundaryCondition(Dirichlet(mesh, 0, 0))

Ux.setBoundaryCondition(Dirichlet(mesh, 1, 0))
Uy.setBoundaryCondition(Dirichlet(mesh, 1, 0))

Ux.setBoundaryCondition(Dirichlet(mesh, 2, 1))
Uy.setBoundaryCondition(Dirichlet(mesh, 2, 0))

Ux.setBoundaryCondition(Dirichlet(mesh, 3, 0))
Uy.setBoundaryCondition(Dirichlet(mesh, 3, 0))

p.setBoundaryCondition(Neuman(mesh, 0, 0))
p.setBoundaryCondition(Neuman(mesh, 1, 0))
p.setBoundaryCondition(Neuman(mesh, 2, 0))
p.setBoundaryCondition(Neuman(mesh, 3, 0))

np.set_printoptions(precision=3)


from fvMatrix import fvMatrix
einterp = EdgeField.interp

Mxd, Fxd = laplace(Ux, fvMatrix)
Myd, Fyd = laplace(Uy, fvMatrix)

Mpd, Fpd = laplace(p, fvMatrix)

edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
print edgeU.data
#print mesh.Se

phi = edgeU.dot(mesh.Se)

Mxc, Fxc = div(phi, Ux, fvMatrix)                      # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div
Myc, Fyc = div(phi, Uy, fvMatrix)


viscosity = 0.1

momX_M = Mxc - Mxd * viscosity
momY_M = Myc - Myd * viscosity

gradP = grad(p)

momX_F = Fxd - Fxd * viscosity - gradP[:, 0]*mesh.cell_area
momY_F = Fyd - Fyd * viscosity - gradP[:, 1]*mesh.cell_area


from scipy.sparse.linalg.isolve.iterative import bicgstab


Ux.data = bicgstab(A=momX_M.sparse, b=momX_F, x0=Ux.data)[0]
Uy.data = bicgstab(A=momY_M.sparse, b=momY_F, x0=Uy.data)[0]

# Results = list()
# Tn = T.data.reshape((n, n))
# Results.append(Tn)
#
#

#
#
# for iter in range(int(nt)):
#     print 'time iteration:',iter
#
#     F = Fconst + T.data
#     # T.data = np.array(np.linalg.solve(M, F))
#     T.data = bicgstab(A=M.sparse, b=F, x0=T.data)[0]
#
#     T.data = T.data.reshape((len(F), 1))
#
#     Tn = T.data.reshape((n, n))
#     Results.append(Tn)
#
# # Animate results:
# #animate_contour_plot(Results, skip=10, repeat=False, interval=75, dataRange=[0, 10])
#
# gradT = grad(T)
#
# from interpolacja import inter
# from matplotlib.pyplot import quiver
# from matplotlib.pyplot import contourf
# from numpy import meshgrid
#
#
# animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1,n+1))], skip=10, repeat=False, interval=75, dataRange=[0, 10])
#
# q = quiver(mesh.cell_centers[:,0], mesh.cell_centers[:,1], gradT[:, 0], gradT[:, 1])
#
# plt.show()



#draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
#draw_edges(mesh.xy, mesh.list_kr)


# start = time.clock()
# print ">>>> Solved in ", time.clock()-start
