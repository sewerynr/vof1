from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Zmienne w czasie zapis macierzy zadkich jako "wektory" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DlPrzX = 1.; DlPrzY = 1.

n = 50                                                   # ilosc podzialow
m = n
dx = DlPrzX/m
dy = DlPrzY/n

dt = 0.001                                                 # CFL u*dt/dx <= 1
tp = 0
tk = 0.1

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)

viscosity = 0.01

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


for i in range(300):
    print "iter", i
    edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))          # interpoluje krawedziowe pole predkosci
    phi = edgeU.dot(mesh.normals)                               # do RC 
    gradP = grad(p)

    Mxc, Fxc = div(phi, Ux)             # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div
    Myc, Fyc = div(phi, Uy)

    # ukladanie rownan transportu dla pedu
    momX_M = Mxc - Mxd * viscosity
    momY_M = Myc - Myd * viscosity
    # RHS:
    momX_F = -(Fxc - Fxd * viscosity) - gradP[:, 0] * mesh.cells_areas
    momY_F = -(Fyc - Fyd * viscosity) - gradP[:, 1] * mesh.cells_areas

    print "Initial continuity residual:", np.linalg.norm(edgeDiv(phi))      # spr spelnienie RC w obszarze

    # print "Initial Ux residual:", np.linalg.norm(momX_M.dot(Ux.data) - momX_F)
    # print "Initial Uy residual:", np.linalg.norm(momY_M.dot(Uy.data) - momY_F)

    momX_M.relax(0.7)
    momY_M.relax(0.7)

    # 2.  rozwiazuje rownania pedu:
    # oblicz i uaktualniam pole predkosci (setValues()):

    Ux.setValues(bicgstab(A=momX_M.sparse, b=momX_F, x0=Ux.data, tol=1e-8)[0])
    Uy.setValues(bicgstab(A=momY_M.sparse, b=momY_F, x0=Uy.data, tol=1e-8)[0])

    A = np.array([np.array(momX_M.diag), np.array(momY_M.diag)]).T

    coeffFieldX = SurfField(mesh, Neuman)
    coeffFieldX.setValues(mesh.cells_areas / A[:, 0])

    coeffFieldY = SurfField(mesh, Neuman)
    coeffFieldY.setValues(mesh.cells_areas / A[:, 1])

    coeffEdge = EdgeField.vector(EdgeField.interp(coeffFieldX), EdgeField.interp(coeffFieldY))          # interpoluje wsp ze srodkow komorek na fejsy

    edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))      # pole predkosci na krawedziach my mamy w centrach (tak ukladam i rozwiazuje uklad rownan)
    phi = edgeU.dot(mesh.normals)

    # ukladam rownania poprawki cisnienia
    Mpd, Fpd = laplace1(coeffEdge, p)
    Fpd = -Fpd + edgeDiv(phi)

    pressP = SurfField(mesh, Neuman)
    # rozwiazuje rownanie poporawki cisinenia
    pressP.setValues(bicgstab(A=Mpd.sparse, b=Fpd, x0=p.data, tol=1e-8)[0])

    # uaktualniam pole cisnienia
    p.setValues(p.data + (1. - 0.7) * pressP.data)

    # licze gradienty w srodkach komorek z obliczonych cisnieni w srodkach kom.
    gradP = grad(pressP)

    # uaktualniam pole predkosci (*) korzystajac z uproszczonego rownania pedu
    Ux.setValues(Ux.data - gradP[:, 0] * mesh.cells_areas / A[:, 0])
    Uy.setValues(Uy.data - gradP[:, 1] * mesh.cells_areas / A[:, 1])

animate_contour_plot([Ux.data.reshape((n, m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("Ux")

animate_contour_plot([Uy.data.reshape((n, m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("Uy")

animate_contour_plot([p.data.reshape((n, m))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
plt.title("p")

Umag = np.sqrt(np.multiply(Ux.data, Ux.data) + np.multiply(Uy.data, Uy.data))
animate_contour_plot([inter(mesh.xy, mesh.cells, Umag).reshape((n+1, m+1))], skip=1, nLevels=20, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
#
from matplotlib.pyplot import quiver
q = quiver(mesh.cell_centers[:, 0], mesh.cell_centers[:, 1], Ux[:], Uy[:])
plt.title("magU")
plt.show()



