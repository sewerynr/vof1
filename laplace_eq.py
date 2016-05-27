from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

DlPrzX = 1.; DlPrzY = 1.

n = 20          # ilosc podzialow

dx = DlPrzX/n
dy = DlPrzY/n

x0, y0, dl = (0, 0, 0)


# wsp_wezl, cells, bounduary = siatka_regularna_prost(n, dx, dy, x0, y0)  czyli metoda siatka_regularna_prost zwroci to co Mesh potrzebuje czyli nodes, cells, boundary

node_c, cells, bound = siatka_regularna_prost(n, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)          # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

T = SurfField(mesh)    # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))


T.data = 0          # [:] do listy przypisze wartosc 0, samo = przypisze inny obiekt   przypisuje wszedzie wartosc 0

Tdir = 1

# neuman    def __init__(self, mesh, bId, derivativeValue )
#T.setBoundaryCondition(Neuman(mesh, 0, -5))         # zero odpowiada zerowej krawedzi pobiera obiekt klasy Dirichlet (wywoluje go i tworzy)
T.setBoundaryCondition(Dirichlet(mesh, 0, 1))

T.setBoundaryCondition(Neuman(mesh, 1, 0))

#T.setBoundaryCondition(Dirichlet(mesh, 1, Tdir))
#T.setBoundaryCondition(Neuman(mesh, 2, -10))

T.setBoundaryCondition(Dirichlet(mesh, 2, Tdir))
T.setBoundaryCondition(Neuman(mesh, 3, 0))

#T.setBoundaryCondition(Dirichlet(mesh, 3, Tdir))

# d = Dirichlet(mesh, 3, Tdir)
# d[:] = 1            # domyslnie zainicjalizowana zerami (dziedziczy po EdgeField) tu zapisuje te krawedz wartoscia = 1
# T.setBoundaryCondition(d)


M, F = laplace(T)    # po prostu ukladanie macierzy i wektora prawych stron
Mc, Fc = div(generate_phi(mesh), T)    # po prostu ukladanie macierzy i wektora prawych stron

Nt = 20
dt = 1./(Nt-1)

Mass = np.zeros((n*n, n*n))
Mass[range(n), range(n)] = 1./dt

M = M + Mc
F = F + Fc



# def anim(t):
#     Fc = np.copy(F) + T.data/dt
#     X = solve(M, F)
#     T.setValues(X)
#     draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
#
# from matplotlib.animation import FuncAnimation
# fig = plt.figure()
# animation = FuncAnimation(fig, anim, frames=range(Nt))
#
# plt.show()


A = solve(M, F)         #  rozw uklad rownan


T.setValues(A)          # pole temp do aktualizacji wartosci temp w WB Neumana

# print T.data
# print A

draw_values_edges(mesh.xy, mesh.cells, mesh.list_kr, T, n, DlPrzX, DlPrzY, Tdir)
# draw_edges(mesh.xy, mesh.list_kr)

# print "poch", (T.boundaries[0].data[n/2] - T[n/2])/(dy/2)
#
# plt.figure()
# plt.plot(T.boundaries[0].data,color="blue")
# plt.plot(T[:n], color="green")
# plt.plot(T[n:2*n], color="red")
# plt.show()


#index_draw_cell( mesh.xy, mesh.cells)



#draw_values_centers(dx, dy, DlPrzX, DlPrzY, n, A)