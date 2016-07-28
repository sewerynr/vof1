import numpy as np
from mesh import *
from fvm1 import *



class BoundaryField:          # pole dla krawedzi zawiera inf co na krawedziach w konstr pobiera siatke mesh oraz numer krawedzi boundaryId
    def __init__(self, mesh, boundaryId):
        self.data = np.array([0.] * len(mesh.boundaries[boundaryId]))  # inicjuje zerami pobierajac dlugosc z mesh.boundaries[boundaryId]
        self.id = boundaryId       # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh
        #print self.mesh.boundaries
#!!!!!!!!!!!!!!!!!! te funkcje maja byc w warunkach brzeg ale tu sa nie koniecznie potrzebne (pokazuja strukture WB)!!!!!!!!!!!!!!!!!!!

    def upadate(self, solution_field):   # ta metoda musi tu byc odwolujemy sie do tej ktora nic nie robi dla Dirichleta, a nadpisujem metoda lokalna dla Neumanna
        '''
        Placeholder for future updates of values at boundary. It will be called after field solution change
        :return: None
        '''
        pass

    def insertDiffusiveFlux(self, edgeFieldCoeff, EqMat, Rhs):
        '''
        Function which sets specific boundary condition into provided equation system. Core function of this class
        :param EqMat: Equation system matrix to be modified by this boundary condition
        :param Rhs: Rhs vector of equation system to be modified by this boundary condition
        :return: None
        '''
        pass

    def sInsertDiffusiveFlux(self, data, indices, rhs):
        pass


    def insertConvectiveFlux(self, phi, EqMat, Rhs):            # inne dla Neuman inne dla Dirichet
        pass


    def setUniformValue(self, val):
        '''
        Sets uniform value along whole boundary in respective field
        :param val: value which should be set to all edges
        :return: None
        '''
        self.data = val                                        # przypisze do tego po czym dziedziczy czyli po np.array

    def setField(self, surfaceFieldObject):
        '''
        Joins boundary conditions with its filed.
        :param surfaceFieldObject: Surface Field type object (storing cell center solutions)
        :return: None
        '''
        self.field = surfaceFieldObject

    def __getitem__(self, item):
        return self.data[item]



class Dirichlet(BoundaryField):                                         # clasa dla kazdej krawedzi o warunku dirich
    def __init__(self, mesh, bId, fi=0):
        BoundaryField.__init__(self, mesh, bId)
        self.data[:] = fi
        self.id = bId                                               # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh

    def __setitem__(self, key, value):
        return self.data.__setitem__(key, value)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def insertDiffusiveFlux(self, edgeFieldCoeff, EqMat, Rhs):                      # dopisuje do macierzy_K WB
        for i, _ in enumerate(self.mesh.boundaries[self.id]):
            id_edge = self.mesh.boundaries[self.id, i]              # indeks krawedzi w WB
            field = self.field
            c = self.mesh.list_kr[id_edge, 2]                       # pobierz wlascicela tej krawedzi

            pC = sum(field.mesh.xy[field.mesh.cells[c], :]) / len(field.mesh.cells[c])              # srodek komorki
            nr_kr_1, nr_kr_2 = field.mesh.list_kr[id_edge, 0], field.mesh.list_kr[id_edge, 1]
            pK = (field.mesh.xy[nr_kr_1, :] + field.mesh.xy[nr_kr_2, :]) / 2                        # srodek scianki to samo co ck1 i ck2 ale od razu dla 2 wsp srodka sciany [x,y]

            CK = pK - pC
            Snorm = field.mesh.Se[id_edge]
            a = edgeFieldCoeff.data[id_edge] * np.dot(CK, Snorm) / np.dot(CK, CK)

            EqMat[c, c] += - a
            Rhs[c] -= self.data[i] * a                 # wartosc na brzegu przemnozona przez wspolczynnik. -= bo przerzucamy do wekt. prawych stron


    def insertConvectiveFlux(self, EqMat, Rhs, phi):
        for i, _ in enumerate(self.mesh.boundaries[self.id]):
            id_edge = self.mesh.boundaries[self.id, i]    # indeks krawedzi w WB
            #field = self.field
            c = self.mesh.list_kr[id_edge, 2]             # pobierz wlascicela tej krawedzi
            edgeLen = self.mesh.eLengths[id_edge]
            phiEdge = phi[id_edge]
            if phiEdge > 0:                               # gdy wylatuje
                EqMat[c, c] -= phiEdge * edgeLen
            else:                                         # gdy wlatuje
                Tbrzeg = self.data[i]
                Rhs[c] += Tbrzeg * phiEdge * edgeLen


class Neuman(BoundaryField):                                        # klasa dla kazdej krawedzi o warunku neuman
    def __init__(self, mesh, bId, derivativeValue = 0 ):            # derivativeValue po prostu zadana wart poch na krawedzi WB
        BoundaryField.__init__(self, mesh, bId)

        self.deriv = derivativeValue
        self.id = bId                                           # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh

    # mamy pochodna (wartosci szukanej) na krawedzi ale nie wiemy jaka sama wartosci wiec do np wizualizacji przyda nam sie wartosc rozwiazaznia (calka z poch)
    def upadate(self, rozw_ukl_row):                            # gdy sie rozwiaze to ma uaktualnic sama siebie ta metoda te klase
        #print self.mesh.boundaries_points
        self.extrapolate(rozw_ukl_row)

        # print (self.mesh.boundaries[0, :].size)               # zwraca rozmiar pod macierzy (rozmiar wiersza)
        # dane numery pkt w tablicy mesh.xy[12] = x12, y12 mozemy odczytac wsp pkt z ktorych jest krawedz i policzyc sr krawedzi
        '''
        TODO - calc. edge values from cells canter and derivative value on edge
        :return: None
        '''
        pass
    # do upadate
    def extrapolate(self, sol):                                      # i to numer krawedzi w WB
        for i, id_edge in enumerate(self.mesh.boundaries[self.id, :]):
            c = self.mesh.list_kr[id_edge, 2]                                            # pobierz wlascicela tej krawedzi
            cc1 = sum(self.mesh.xy[self.mesh.cells[c], :]) / len(self.mesh.cells[c])  # srodek komorki  [x,y]

            n = self.mesh.edge_normal(id_edge)                   # normalny do wektora (krawedzi WB) ale o jego dlugosci (trzeba pobrac numer pod krawedzi )
            n = n / np.linalg.norm(n)                            # wersor normalny (podzielony przez swoja dlugosc)
            sr_pod_krawBC = sum([self.mesh.xy[nid] for nid in self.mesh.list_kr[id_edge, :2]])/2
            e = sr_pod_krawBC - cc1
            ds = np.dot(n, e)                                    # wektor od srodka komorki do srodka sciany
            elenSq = np.dot(e, e)
            flux = self.deriv
            Tsr = sol[c]
            Tbrzeg = Tsr - flux*elenSq/ds
            self.data[i] = Tbrzeg
            #print self.data[i], self.deriv, ds


    def insertDiffusiveFlux(self, edgeFieldCoeff, EqMat, Rhs):  # pobiera macierz K i wektor pr stron
        #print len(self.mesh.boundaries[self.id])
        for i in range(len(self.mesh.boundaries[self.id])):
            id_edge = self.mesh.boundaries[self.id, i]                          # indeks krawedzi w WB
            c = self.field.mesh.list_kr[id_edge, 2]                             # indeks wlasciciela do niego dopicac w rhs
            Rhs[c] += edgeFieldCoeff.data[id_edge] * self.deriv * self.field.mesh.eLengths[id_edge]      #  dodac razy dlugosc

    def insertConvectiveFlux(self, EqMat, Rhs, phi):
        for i, _ in enumerate(self.mesh.boundaries[self.id]):
            id_edge = self.mesh.boundaries[self.id, i]                          # indeks krawedzi w WB
            field = self.field
            c = self.mesh.list_kr[id_edge, 2]                                   # pobierz wlascicela tej krawedzi
            edgeLen = self.mesh.eLengths[id_edge]

            phiEdge = phi[id_edge]

            if phiEdge > 0:                                                     # wylatuje
                EqMat[c, c] -= phiEdge * edgeLen
            else:                                                               # gdy wlatuje
                Tbrzeg = self.data[i]                                           # powinno byc z zewnetrznej krawedzi nie istniejacej bo wlatuje z zewnatrz
                Rhs[c] -= Tbrzeg * phiEdge * edgeLen
                # nalezy poprawic to powyzej, tak aby z pochodnej i wartosc w srodku sasiedniej
                # komorki pisac rownanie na wartosc wielkosci wlatujacej do srodka
                # to powinno doprowadzic do wstawienia czegos do macierzy i czegos do wekt. pr. stron




class symmetry(Neuman):                                                                # klasa dla kazdej krawedzi o warunku neuman o wartosci poch w kier normalym = 0
    def __init__(self, mesh, bId):
        Neuman.__init__(self, mesh, bId, np.array([0.] * len(mesh.boundaries[bId])))    # konstruktor pierwotnej klasy po ktorej dziedziczy, derrevativeValue przyjmuje teraz liste o wartosicach 0 i rozmiarze dlugosci mesh.boundaries
        self.id = bId                                                                   # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh



class SurfField:
    def __init__(self, mesh, bcGenerator=BoundaryField):                        # pobiera mesh a z nim jego rozmiar i boundaries czyli WB
        self.data = np.zeros(mesh.n)
        self.boundaries = []

        for i, _ in enumerate(mesh.boundaries):
            bc = bcGenerator(mesh, i)                              # wywoluje BonduaryField z meshem i z numerem krawedzi
            bc.setField(self)
            self.boundaries.append(bc)

        self.mesh = mesh

    def __setitem__(self, key, value):
        return self.data.__setitem__(key, value)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def setBoundaryCondition(self, bcObject):           # to zada warunki brzegowe bcObject( np class/object nemuam) na surface field
        bcObject.setField(self)                         # na przyklad  neuman.setFirld(self)
        self.boundaries[bcObject.id] = bcObject         # przypisuje do tej listy (***) (tylko z numerami krawedzi WB) zadane przez urzytkownika WB

    def updateBoundaryValues(self):
        for b in self.boundaries:            # b kolejny warunek brzegowy tu 4 warunki
            b.upadate(self.data)             # przekaz pole temp do funkcji upadate

    def setValues(self, lista_wart):           # pobiera rozwiazanie ukladu rownan (pole temp) - uaktualnia wartosci pola temp tak aby na krawedziach WB neumana byly temp nie strumienie ciepla (do wizualizacji)
        self.data = lista_wart                 # rozwiazanie tu tablica "A"
        self.updateBoundaryValues()

    def apply_bc_diffusiveFlux(self, edgeFieldCoeff, EqMat, Rhs):
        for b in self.boundaries:
            b.insertDiffusiveFlux(edgeFieldCoeff, EqMat, Rhs)


    def apply_bc_convectiveFlux(self, EqMat, Rhs, phi):          # ma wrzucic ten warunek na macierz M i rhs
        for b in self.boundaries:                                # petla po 4 war brzeg (te brzegi juz zapisalem uzywajac setBoundaryCondition
            b.insertConvectiveFlux(EqMat, Rhs, phi)              # jesli neuman to wywola z klasy neuman jesil dirichlet to z dirichlet zalezy czym jest b



class EdgeField:
    def __init__(self, mesh, dim=1):                             #   wymiar dla T - dim = 1 dla v = 1 do 3
        if dim > 1:
            self.data = np.zeros( (len(mesh.list_kr), dim) )
        else:
            self.data = np.zeros(len(mesh.list_kr))

        self.edgesdef = mesh.list_kr
        self.mesh = mesh

    @staticmethod
    def interp(surfField):                                              # def interp(self, surfField):   nie statyczna musi byc z self

        mesh = surfField.mesh
        efield = EdgeField(mesh)

        for i, kraw in enumerate(mesh.list_kr):                          # interoplacja ze srodkow na krawedz
            if kraw[3] > -1:                                             # jesli ma sasiada
                # przypadek szczegolny odl do srodkow krawedzi nie przeciec wektorow
                vwl = surfField[kraw[2]]
                vsas = surfField[kraw[3]]

                # wsp sr kom wl: mesh.cell_centers[kraw[2]]
                dl = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[3]], mesh.cell_centers[kraw[2]])
                dl = mesh.dl_wekt(dl[0], dl[1])                              # dlugosc miedzy sr komorek

                dlc = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[2]], (mesh.xy[kraw[0]] + mesh.xy[kraw[1]]) / 2)
                dlc = mesh.dl_wekt(dlc[0], dlc[1])

                dlf = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[3]], (mesh.xy[kraw[0]] + mesh.xy[kraw[1]]) / 2)
                dlf = mesh.dl_wekt(dlf[0], dlf[1])

                v = vwl * (dlf / dl) + vsas * (dlc / dl)                    # przypadek szczegolny
                efield.data[i] = v
            else:                                                           # gdy krawedz brzegowa
                for bId, bEdges in enumerate(mesh.boundaries):
                    for eLocalId, e in enumerate(bEdges):
                        if e == i:                                                        # jesli krawedz z kraw jest w danym WB to:
                            efield.data[i] = surfField.boundaries[bId][eLocalId]          # przyjmij to co na WB [bId] i krawedzi [eLocalId]

        return efield

    @staticmethod
    def vector(scalarX, scalarY):
        efield = EdgeField(scalarX.mesh, dim=2)
        efield.data[:, 0] = scalarX.data                # wszystkie wiersze z 0 kolumny
        efield.data[:, 1] = scalarY.data
        return efield

    def dot(self, other):                                # iloczyn skalarny self.data(to co wywolujw) z podanym wektorem (lista)
        ret = EdgeField(self.mesh)

        if isinstance(other, EdgeField):                # jesli numpy ?
            data = other.data
        else:
            data = other

        ret.data = np.sum(self.data[:]*data[:], axis=1, dtype=float)
        return ret


def quadratic_velocity(pc, tanPc, r):
    U = 10.
    if r <= 0.5:
        return tanPc*U*(-(4*r-1.)**2+1.)
    else:
        return np.array([0., 0.])

#  odrazu daje normalne skl pr na sciankach
def generate_phi_r(mesh, velcity_function):           #po krawedziach
    vals = np.zeros(len(mesh.list_kr), dtype=float)
    for i,kraw in enumerate(mesh.list_kr):
        p1, p2 = mesh.xy[kraw[:2]]
        pc = (p1 + p2)/2. - [0.5, 0.5]
        r = np.sqrt(pc.dot(pc))
        tan = np.array([-pc[1], pc[0]])
        tan = tan / np.sqrt(tan.dot(tan))
        U = velcity_function(pc, tan, r)
        vals[i] = U.dot(mesh.normals[i])

    return vals

def generate_phi_1(mesh):
    vals = np.zeros(len(mesh.list_kr))
    for i, kraw in enumerate(mesh.list_kr):
        p1, p2 = mesh.xy[kraw[:2]]              # zczytaj punkty krawedzi  (dwie pierwsze liczby z list_kr) i pobierz ich wsp x , y
        pc = (p1 + p2)/2.                       # srodek krawedzi
        U = [0, 0]
        U = np.array([U[0], U[1]])
        vals[i] = U.dot(mesh.edge_normal(i))    # rzut na normalna do konkretnej krawedzi - phi
    return vals


def constant_velocity(pc, tanPc, r):
    return np.array([10, 0.])

def generate_u(mesh, velcity_function):
    Ux = SurfField(mesh, Dirichlet)
    Uy = SurfField(mesh, Dirichlet)
    # Ux = SurfField(mesh, Neuman)
    # Uy = SurfField(mesh, Neuman)


    vals = np.zeros((len(mesh.cell_centers), 2), dtype=float)

    for c, pc in enumerate(mesh.cell_centers):
        pc = pc - [0.5, 0.5]
        r = np.sqrt(pc.dot(pc))
        tan = np.array([-pc[1], pc[0]])                 # normalna do pc[x, y] = pcn[-y, x]
        tan = tan / np.sqrt(tan.dot(tan))
        vals[c, :] = velcity_function(pc, tan, r)

    Ux.setValues(vals[:, 0])
    Uy.setValues(vals[:, 1])

    return Ux, Uy


def generate_phi_r_2(mesh):
    vals_return = np.zeros(len(mesh.list_kr))
    vals = np.array([[0.] * 2] * len(mesh.cells))     # lista wartosci 0 i dl odpowiadajacej ilosci krawedzi
    # generuje predkosci w srodkach komorek
    for i, cell in enumerate(mesh.cells):             # dla kazdej krawedzi
        pc = mesh.cell_centers[i] - [0.5, 0.5]
        r = np.square(pc.dot(pc))
        tan = np.array([-pc[1], pc[0]])               # normalna do pc[x, y] = pcn[-y, x]
        tan = tan / np.sqrt(tan.dot(tan))             # kierunek normalej
        U = tan * r * 3200                             # razy wsp.
        #print len(mesh.list_kr)
        vals[i] = U                                   # wektor predkosci w srodkach komorek

    # trzeba wyinterpolowac predkosci na scianki

    # trzeba warunek na krawedzie brzegowe
    for i, kraw in enumerate(mesh.list_kr):
        if kraw[3] > 0:    # jesli ma sasiada
            # przypadek szczegolny odl do srodkow krawedzi nie przeciec wektorow
            # zczytuje z listy_kr wlasciciela i sasiada i pobieram predkosci w ich srodkach komorek jako wektor [x, y]
            vwl = vals[kraw[2]]
            vsas = vals[kraw[3]]
            #wsp sr kom wl: mesh.cell_centers[kraw[2]]
            dl = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[3]], mesh.cell_centers[kraw[2]])
            dl = mesh.dl_wekt(dl[0], dl[1])                      # dlugosc miedzy sr komurek

            dlc = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[2]], (mesh.xy[kraw[0]] + mesh.xy[kraw[1]])/2)
            dlc = mesh.dl_wekt(dlc[0], dlc[1])

            dlf = mesh.wsp_wekt_z_wsp(mesh.cell_centers[kraw[3]], (mesh.xy[kraw[0]] + mesh.xy[kraw[1]]) / 2)
            dlf = mesh.dl_wekt(dlf[0], dlf[1])

            v = vwl*(dlf/dl) + vsas*(dlc/dl)                      # przypadek szczegolny
            vals_return[i] = v.dot(mesh.edge_normal(i))           # rzut na normalna do konkretnej krawedzi
            print vals_return[i]
        else:                                                     # gdy krawedz brzegowa
            #vals_return[i] = 0.
            pass                                                  # vals_return[i] = 0. nie trzeba bo zainicjalizowana zerami

    #print vals_return
    return vals_return


def grad(surfField):                     # Green - Gauss

    import numpy as np

    mesh = surfField.mesh

    efield = EdgeField.interp(surfField)

    cellGrad = np.zeros((len(mesh.cells), 2), dtype=float)

    for e, (defE, valE) in enumerate(zip(mesh.list_kr, efield.data)):
        P1 = mesh.xy[defE[0]]
        P2 = mesh.xy[defE[1]]
        dP = P2 - P1
        cellGrad[defE[2], :] += [-valE * dP[1], valE * dP[0]]

        if defE[3] >= 0:
            cellGrad[defE[3], :] += [valE * dP[1], -valE * dP[0]]

    cellGrad[:, 0] = cellGrad[:, 0] / mesh.cells_areas
    cellGrad[:, 1] = cellGrad[:, 1] / mesh.cells_areas

    return cellGrad


def edgeDiv(edgeField):
    mesh = edgeField.mesh
    ret = np.zeros(len(mesh.cells), dtype=float)

    for e, eDef in enumerate(mesh.list_kr):
        ev = mesh.Se[e]
        # dS = np.sqrt(ev.dot(ev))
        val = edgeField.data[e]

        ret[eDef[2]] += val#*dS

        if eDef[3] >= 0:
            ret[eDef[3]] -= val#*dS

    return ret

#
#
# def values(obiekt):
#     return obiekt._sol


# class OpenFoamWriter:
#     def __init__(self, mesh, boundaryId, field):
#         self.id = boundaryId  # zapisuje numer krawedzi pod oznaczeniem id
#         self.mesh = mesh
#         self.field = field
#
#     def write_mesh(self):
#
#         plik = open("mesh.txt", "w")
#
#         plik.write(np.array_str(self.mesh.cells))
#
#         plik.write(np.array_str(self.mesh.xy))
#
#         if not plik.closed:
#             print "closing file ..."
#             plik.close()
#
#         if not plik.closed:
#             print "closing file ..."
#             plik.close()
#
#
#     def write_solution(self):
#         plik = open("solution.txt", "w")
#         if not plik.closed:
#             print "closing file ..."
#             plik.close()






