from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

DlPrzX = 1.; DlPrzY = 1.

n = 40                                                     # ilosc podzialow

dx = DlPrzX/n
dy = DlPrzY/n

x0, y0, dl = (0, 0, 0)


# wsp_wezl, cells, bounduary = siatka_regularna_prost(n, dx, dy, x0, y0)  czyli metoda siatka_regularna_prost zwroci to co Mesh potrzebuje czyli nodes, cells, boundary

node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                         # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

T = SurfField(mesh)                                       # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))

T.data = 1                                                # [:] do listy przypisze wartosc 0, samo = przypisze inny obiekt   przypisuje wszedzie wartosc 0

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


M, F = steady_laplace(T)                                        # ukladanie macierzy i wektora prawych stron laplace
Mc, Fc = steady_div(generate_phi_r(mesh, quadratic_velocity), T)                         # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div
#Mc, Fc = steady_div(generate_phi_1(mesh), T)                      # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div

np.set_printoptions(precision=3)

#print Mc
M = M*(0) + Mc*(1)

F = F*(0) + Fc*(1)
np.set_printoptions(precision=3)
#print Fc

pkt1 = n/2 + n*n/2                       # pkt srodek
pkt2 = n/2 + n*4                        # srodek 4 wiersze od spodu
pkt3 = n/2 + n*(n-4)                     # srodek 4 wiersze od gory

#F[pkt1] = -10
F[pkt2] = -10
F[pkt3] = -11

A = solve(M, F)                        # rozw uklad rownan

T.setValues(A)                         # pole temp do aktualizacji wartosci temp w WB Neumana

print T.data
# print A

draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
# draw_edges(mesh.xy, mesh.list_kr)

