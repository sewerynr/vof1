DlPrzX = 1.
DlPrzY = 1.

n = 34                                                   # ilosc podzialow
m = n
dx = DlPrzX/m
dy = DlPrzY/n

dt = 0.001                                              # CFL u*dt/dx <= 1
tp = 0
tk = 0.1

nt = (tk - tp)/dt

x0, y0, dl = (0, 0, 0)

viscosity = 0.01
mom_relax = 0.7

from mesh import Mesh
from numpy.linalg import solve
from fvm1 import *
from field import *
from interpolacja import *
from fvMatrix import fvMatrix
from scipy.sparse.linalg.isolve.iterative import bicgstab

einterp = EdgeField.interp

node_c, cells, bound = siatka_regularna_prost(n, m, dx, dy, x0, y0)

mesh = Mesh(node_c, cells, bound)                         # 1. tworzy obiekt mesh klasy Mesh, 2. generujemy siatke dla tego obiektu funkcja siatka_reg...

Ux = SurfField(mesh, Dirichlet)                                       # tworzy obiekt klasy SurfField, pobierajacy obirkt mesh klasy Mesh ( na tej siatce ma tworzyc i przechowywac rozwiazanie (wartosci))
Uy = SurfField(mesh, Dirichlet)
p = SurfField(mesh, Neuman)

Ux.setBoundaryCondition(Dirichlet(mesh, 0, 0.))
Ux.setBoundaryCondition(Dirichlet(mesh, 1, 1.))
Ux.setBoundaryCondition(Dirichlet(mesh, 2, 0.))
Ux.setBoundaryCondition(Dirichlet(mesh, 3, 1.))

# p.setBoundaryCondition(Dirichlet(mesh, 3, 10))
# p.setBoundaryCondition(Dirichlet(mesh, 1, -10))

np.set_printoptions(precision=3)

Mxd, Fxd = laplace(viscosity, Ux)
Myd, Fyd = laplace(viscosity, Uy)

for i in range(80):
    print "iter", i
    edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))          # interpoluje krawedziowe pole predkosci z nowych wartosci w kazdym kroku
    phi = edgeU.dot(mesh.normals)                               # do RC phi = v * n
    gradP = grad(p)                                             # obl gradienty cisnienia w srodkach komorek z nowych wart cisnienia w kazdym kroku

    Mxc, Fxc = div(phi, Ux)             # ukladanie macierzy i wektora prawych stron, dostaje D i Rhs z div
    Myc, Fyc = div(phi, Uy)

    # ukladanie rownan transportu dla pedu
    momX_M = Mxc - Mxd
    momY_M = Myc - Myd
    # RHS:
    momX_F = -(Fxc - Fxd) - gradP[:, 0] * mesh.cells_areas
    momY_F = -(Fyc - Fyd) - gradP[:, 1] * mesh.cells_areas

    print "Initial continuity residual:", np.linalg.norm(edgeDiv(phi))      # spr spelnienie RC w obszarze
    print "Initial Ux residual:", np.linalg.norm(momX_M.dot(Ux.data) - momX_F)
    print "Initial Uy residual:", np.linalg.norm(momY_M.dot(Uy.data) - momY_F)

    momX_M.relaxM1(momX_F, Ux, mom_relax)
    momY_M.relaxM1(momY_F, Uy, mom_relax)

    # 2.  rozwiazuje rownania pedu:
    # oblicz i uaktualniam pole predkosci (setValues()):

    Ux.setValues(bicgstab(A=momX_M.sparse, b=momX_F, x0=Ux.data, tol=1e-8)[0])
    Uy.setValues(bicgstab(A=momY_M.sparse, b=momY_F, x0=Uy.data, tol=1e-8)[0])

    edgeU = EdgeField.vector(einterp(Ux), einterp(Uy))       # pole predkosci na krawedziach my mamy w centrach (tak ukladam i rozwiazuje uklad rownan)
    phi = edgeU.dot(mesh.normals)

    # R-H

    A = np.array([np.array(momX_M.diag), np.array(momY_M.diag)]).T

    coeffFieldX = SurfField(mesh, Neuman)
    coeffFieldX.setValues(mesh.cells_areas / A[:, 0])

    coeffFieldY = SurfField(mesh, Neuman)
    coeffFieldY.setValues(mesh.cells_areas / A[:, 1])

    coeffEdge = EdgeField.vector(EdgeField.interp(coeffFieldX), EdgeField.interp(coeffFieldY))          # interpoluje wsp ze srodkow komorek na fejsy

    # ukladam rownania poprawki cisnienia
    Mpd, Fpd = laplace(coeffEdge, p)           # dla nowych wartosci predkosci juz policzonych
    Fpd = -Fpd + edgeDiv(phi)                  # plus RC

    pressP = SurfField(mesh, Neuman)
    # rozwiazuje rownanie poporawki cisinenia
    pressP.setValues(bicgstab(A=Mpd.sparse, b=Fpd, x0=p.data, tol=1e-8)[0])

    # uaktualniam pole cisnienia
    p.setValues(p.data + (1 - mom_relax) * pressP.data)

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
animate_contour_plot([inter(mesh.xy, mesh.cells, Umag).reshape((n+1, m+1))], skip=1, dataRange= [0,2], nLevels=25, repeat=False, interval=75, diff=viscosity, adj=0, nN=n)
#
from matplotlib.pyplot import quiver
q = quiver(mesh.cell_centers[:, 0], mesh.cell_centers[:, 1], Ux[:], Uy[:])
plt.title("magU")
plt.show()



