from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

DlPrzX = 1.; DlPrzY = 1.

n = 10       # ilosc podzialow

dx = DlPrzX/n
dy = DlPrzY/n

x0 = 0; y0 = 0; d1 = 0

#   wsp_wezl, cells, bounduary = siatka_regularna_prost(n, dx, dy, x0, y0)  czyli metoda siatka_regularna_prost zwroci to co Mesh potrzebuje czyli nodes, cells, boundary

node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)          # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

T = SurfField(mesh)    # tworzy obiekt klasy SurfField , pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))


T.data = 0 # [:] do listy przypisze wartosc 5, samo = przypisze inny obiekt   przypisuje wszedzie wartosc 5

Tdir = 1

# neuman    def __init__(self, mesh, bId, derivativeValue )
T.setBoundaryCondition(Neuman(mesh, 0, -10))     #  zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)

#T.setBoundaryCondition(Neuman(mesh, 1, -10))
T.setBoundaryCondition(Dirichlet(mesh, 1, Tdir))

T.setBoundaryCondition(Dirichlet(mesh, 2, Tdir))
T.setBoundaryCondition(Dirichlet(mesh, 3, Tdir))

# d = Dirichlet(mesh, 3, Tdir)
# d[:] = 1            # domyslnie zainicjalizowana zerami (dziedziczy po EdgeField) tu zapisuje te krawedz wartoscia = 1
# T.setBoundaryCondition(d)

T.mesh.cells


M, F = laplace(T)    # po prostu ukladanie macierzy i wektora prawych stron


T.apply_bc(M, F)    # dodawanie warunkow brzegowych do macierzy i wektora PS funkcja pochodzi z klasy T czyli SurfField

# print M, F

A = solve(M, F)

T.setValues(A)        #  rozw uklad rownan

# print T.data
# print A

draw_values_edges(mesh.xy, mesh.cells, T, n, DlPrzX, DlPrzY, Tdir)

# index_draw_cell( mesh.xy, mesh.cells)

#draw_edges(mesh.xy, mesh.list_kr)


draw_values_centers(dx, dy, DlPrzX, DlPrzY, n, A)