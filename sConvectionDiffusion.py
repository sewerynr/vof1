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


T = SurfField(mesh)                                       # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))

#print mesh.cell_centers


Tdir = 1
TdirWB = 0

#T.setBoundaryCondition(Neuman(mesh, 0, 0))               # zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)
T.setBoundaryCondition(Dirichlet(mesh, 0, TdirWB))

#T.setBoundaryCondition(Neuman(mesh, 1, 0))
T.setBoundaryCondition(Dirichlet(mesh, 1, TdirWB))

#T.setBoundaryCondition(Neuman(mesh, 2, 0))
T.setBoundaryCondition(Dirichlet(mesh, 2, TdirWB))

#T.setBoundaryCondition(Neuman(mesh, 3, 0))              # symetria na krawedzi 3 (4)
T.setBoundaryCondition(Dirichlet(mesh, 3, TdirWB))

# d = Dirichlet(mesh, 3, TdirWB)
# d[:] = 1                                               # domyslnie zainicjalizowana zerami (dziedziczy po EdgeField) tu zapisuje te krawedz wartoscia = 1
# T.setBoundaryCondition(d)

np.set_printoptions(precision=3)

#print Fc

pkt1 = n/2 + n*n/2                       # pkt srodek
pkt2 = n/2 + n*5                         # srodek 4 wiersze od spodu
pkt3 = n/2 + n*(n-5)                     # srodek 4 wiersze od gory


from fvMatrix import fvMatrix

Md, Fd = laplace(T, fvMatrix)#sLaplace(T)                                        # ukladanie macierzy i wektora prawych stron laplace

Mc, Fc = div(generate_phi_r(mesh, quadratic_velocity), T, fvMatrix)                      # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div

M = fvMatrix.diagonal(mesh) - Md*(dt*0.2) + Mc*dt

Fconst = (-Fd*0.2 + Fc)*dt

Fconst[pkt2] += 20000*dt
Fconst[pkt3] += 20000*dt


T.data[:] = 0                          # War. Pocz. # [:] do listy przypisze wartosc 0, samo = przypisze inny obiekt przypisuje wszedzie wartosc 0

# do polowy przypisuje wartosc T = 10
# for i, point in enumerate(T.data):
#     if i < (n**2)/2:
#         T.data[i] = 10


Results = list()
Tn = T.data.reshape((n, n))
Results.append(Tn)


from scipy.sparse.linalg.isolve.iterative import bicgstab


for iter in range(int(nt)):
    print 'time iteration:',iter

    F = Fconst + T.data
    # T.data = np.array(np.linalg.solve(M, F))
    T.data = bicgstab(A=M.sparse, b=F, x0=T.data)[0]

    T.data = T.data.reshape((len(F), 1))

    Tn = T.data.reshape((n, n))
    Results.append(Tn)

# Animate results:
#animate_contour_plot(Results, skip=10, repeat=False, interval=75, dataRange=[0, 10])

gradT = grad(T)

from interpolacja import inter
from matplotlib.pyplot import quiver
from matplotlib.pyplot import contourf
from numpy import meshgrid


animate_contour_plot([inter(mesh.xy, mesh.cells, T.data).reshape((n+1,n+1))], skip=10, repeat=False, interval=75, dataRange=[0, 10])

q = quiver(mesh.cell_centers[:,0], mesh.cell_centers[:,1], gradT[:, 0], gradT[:, 1])

plt.show()



#draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
#draw_edges(mesh.xy, mesh.list_kr)


# start = time.clock()
# print ">>>> Solved in ", time.clock()-start
