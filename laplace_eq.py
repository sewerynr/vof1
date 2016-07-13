from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

DlPrzX = 1.; DlPrzY = 1.

n = 26                                                     # ilosc podzialow

dx = DlPrzX/n
dy = DlPrzY/n

dt = 0.01                                                 # CFL u*dt/dx <= 1
tp = 0
tk = 1

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)


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


M, F = laplace(T, dt)                                        # ukladanie macierzy i wektora prawych stron laplace
Mc, Fc = div(generate_phi(mesh), T, dt)                      # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div

np.set_printoptions(precision=3)

#print M
#M = M*(1.000001) + Mc*(0)

#F = F*(1.000001) + Fc*(0)
np.set_printoptions(precision=3)
#print Fc

pkt1 = n/2 + n*n/2                       # pkt srodek
pkt2 = n/2 + n*8                         # srodek 4 wiersze od spodu
pkt3 = n/2 + n*(n-8)                     # srodek 4 wiersze od gory

F[pkt1] += -100
#F[pkt2] += -100
#F[pkt3] += -100
#print F
# !!!!!!!!!!!!!!!!!!!!!!!!!!!! steady solution !!!!!!!!!!!!!!!!!!!!!!
# T.data[:] = 0                          # War. Pocz. # [:] do listy przypisze wartosc 0, samo = przypisze inny obiekt przypisuje wszedzie wartosc 0
#
# F = F + T.data
#
# A = solve(M, F)                        # rozw uklad rownan
#
# T.setValues(A)                         # pole temp do aktualizacji wartosci temp w WB Neumana

#print T.data
# print A

#draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
# # draw_edges(mesh.xy, mesh.list_kr)



I = np.matrix(np.identity(n*n))
M = I - M

T.data[:] = 0                          # War. Pocz. # [:] do listy przypisze wartosc 0, samo = przypisze inny obiekt przypisuje wszedzie wartosc 0
Fconst = F
F = - Fconst + T.data

Results = list()
Tn = T.data.reshape((n, n))
Results.append(Tn)
for iter in range(int(nt)):
#   print 'time iteration:',iter

    T.data = np.array(np.linalg.solve(M, F))
    F = -Fconst + T.data
    Tn = T.data.reshape((n, n))
    Results.append(Tn)
    # print T

# Animate results:
animate_contour_plot(Results, skip=1, repeat=False, interval=200, dataRange=[0, 10])

#draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
#draw_edges(mesh.xy, mesh.list_kr)