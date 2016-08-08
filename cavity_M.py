from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Zmienne w czasie zapis macierzy zadkich jako "wektory" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DlPrzX = 1.; DlPrzY = 1.

n = 10                                                    # ilosc podzialow
m = n
dx = DlPrzX/m
dy = DlPrzY/n

dt = 0.001                                                 # CFL u*dt/dx <= 1
tp = 0
tk = 0.1

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)

viscosity = 1

node_c, cells, bound = siatka_regularna_prost(n, m, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                         # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

Ux = SurfField(mesh, Dirichlet)                                       # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))
Uy = SurfField(mesh, Dirichlet)
p = SurfField(mesh, Neuman)

Ux.setBoundaryCondition(Dirichlet(mesh, 2, 1.))

np.set_printoptions(precision=3)

from fvMatrix import fvMatrix
einterp = EdgeField.interp

Mxd, Fxd = laplace(viscosity, Ux)
Myd, Fyd = laplace(viscosity, Uy)

from scipy.sparse.linalg.isolve.iterative import bicgstab

edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))  # pole wektorowe predkosci [Ux, Uy] wyinterpolowanych wartosci na krawedzie (ze sr komurek  EdgeField.interp)
phi = edgeU.dot(mesh.normals)  # phi = v n A  gdzie An tu rowne jest dl_krawedzi obruconej o 90 stopni

gradP = grad(p)                         # gradient cisnienia z Green - Gaussa  ( ze sr na krawedzie z krawedzi komorki do jej srodka)

for i in range(400):
    print "iter", i

    Mxc, Fxc = div(phi, Ux)             # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div
    Myc, Fyc = div(phi, Uy)

    # ukladanie rownan transportu dla pedu
    momX_M = Mxc - Mxd * viscosity
    momY_M = Myc - Myd * viscosity
    # RHS:
    momX_F = -(Fxc - Fxd * viscosity) - gradP[:, 0] * mesh.cells_areas
    momY_F = -(Fyc - Fyd * viscosity) - gradP[:, 1] * mesh.cells_areas

    momX_M.relax(0.7)
    momY_M.relax(0.7)

    # rozwiazuje rownania pedu:
    xSol = bicgstab(A=momX_M.sparse, b=momX_F, x0=Ux.data)[0]           # A * x0 = b   =>   momX_M * Ux.data = momX_F
    ySol = bicgstab(A=momY_M.sparse, b=momY_F, x0=Uy.data)[0]

    # uaktualniam pole predkosci (setValues()):
    Ux.setValues(xSol)
    Uy.setValues(ySol)
    print Ux.data
    edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))
    phi = edgeU.dot(mesh.normals)

    # czlon dyfuzyjny dla rown popr cisnina
    Mpd, Fpd = laplace(1, p)

    # do RHS rownania na popr cisnienia dodaj obliczone strumienie (edgeDiv(predkosci) = RC)
    Fpd = -Fpd + edgeDiv(phi)

    # Solve the pressure equation and apply under-relaxation.
    pres = bicgstab(A=Mpd.sparse, b=Fpd, x0=p.data)[0]


    # uaktualnij pole cisnienia
    p.setValues(p.data * 0.7 + 0.3 * pres)

    gradP = grad(p)


animate_contour_plot([Ux.data.reshape((n,m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("Ux")

animate_contour_plot([Uy.data.reshape((n,m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("Uy")

animate_contour_plot([p.data.reshape((n,m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("p")

Umag = np.sqrt(np.multiply(Ux.data, Ux.data) + np.multiply(Uy.data, Uy.data))
animate_contour_plot([inter(mesh.xy, mesh.cells, Umag).reshape((n+1, m+1))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
#
from matplotlib.pyplot import quiver
q = quiver(mesh.cell_centers[:, 0], mesh.cell_centers[:, 1], Ux[:], Uy[:])
plt.title("magU")
plt.show()



