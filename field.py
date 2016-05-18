import numpy as np
from mesh import *
from fvm1 import *



class EdgeField:          # pole dla krawedzi zawiera inf co na krawedziach w konstr pobiera siatke mesh oraz numer krawedzi boundaryId
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

    def manipulateMatrixVector(self, EqMat, Rhs):
        '''
        Function which sets specific boundary condition into provided equation system. Core function of this class
        :param EqMat: Equation system matrix to be modified by this boundary condition
        :param Rhs: Rhs vector of equation system to be modified by this boundary condition
        :return: None
        '''
        pass

    def setUniformValue(self, val):
        '''
        Sets uniform value along whole boundary in respective field
        :param val: value which should be set to all edges
        :return: None
        '''
        self.data = val      # przypisze do tego po czym dziedziczy czyli po np.array

    def setField(self, surfaceFieldObject):             #
        '''
        Joins boundary conditions with its filed.
        :param surfaceFieldObject: Surface Field type object (storing cell center solutions)
        :return: None
        '''
        self.field = surfaceFieldObject

    def __getitem__(self, item):
        return self.data[item]





class Dirichlet(EdgeField):         # clasa dla kazdej krawedzi o warunku dirich
    def __init__(self, mesh, bId, fi=0):
        EdgeField.__init__(self, mesh, bId)
        self.data[:] = fi
        self.id = bId  # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh

    def __setitem__(self, key, value):
        return self.data.__setitem__(key, value)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def manipulateMatrixVector(self, EqMat, Rhs):           # dopisuje do macierzy_K WB
        for i, _ in enumerate(self.mesh.boundaries[self.id]):
            id_edge = self.mesh.boundaries[self.id, i]  # indeks krawedzi w WB
            field = self.field
            c = self.mesh.list_kr[id_edge, 2]      #pobierz wlascicela tej krawedzi
            cc1 = sum(field.mesh.xy[field.mesh.cells[c], :]) / len(field.mesh.cells[c])  # srodek komorki
            nr_kr_1, nr_kr_2 = field.mesh.list_kr[id_edge, 0], field.mesh.list_kr[id_edge, 1]
            ck3 = (field.mesh.xy[nr_kr_1, :] + field.mesh.xy[nr_kr_2, :]) / 2  # srodek scianki to samo co ck1 i ck2 ale od razu dla 2 wsp srodka sciany [x,y]
            ad = wspolczynnik_d(ck3, cc1, field.mesh.xy[nr_kr_1, :], field.mesh.xy[nr_kr_2, :])    #pobiera wsp punkt z ktorych sklada sie krawedz
            EqMat[c, c] -= ad           #   do kom wlasciciela przypisuje wsp
            Rhs[c] += - self.data[i] * ad





class Neuman(EdgeField):            # klasa dla kazdej krawedzi o warunku neuman
    def __init__(self, mesh, bId, derivativeValue ):            # derivativeValue po prostu zadana wart poch na krawedzi WB
        EdgeField.__init__(self, mesh, bId)

        self.deriv = derivativeValue
        self.id = bId            # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh

    # mamy pochodna (wartosci szukanej) na krawedzi ale nie wiemy jaka sama wartosci wiec do np wizualizacji przyda nam sie wartosc rozwiazaznia (calka z poch)
    def upadate(self, rozw_ukl_row):           # gdy sie rozwiaze to ma uaktualnic sama siebie ta metoda te klase
        self.sol = rozw_ukl_row
        print self.mesh.boundaries_points
        self.extrapolate()

        # print (self.mesh.boundaries[0, :].size)        # zwraca rozmiar pod macierzy (rozmiar wiersza)
        # dane numery pkt w tablicy mesh.xy[12] = x12, y12 mozemy odczytac wsp pkt z ktorych jest krawedz i policzyc sr krawedzi
        '''
        TODO - calc. edge values from cells canter and derivative value on edge
        :return: None
        '''
        pass

    def manipulateMatrixVector(self, EqMat, Rhs):  # pobiera macierz K i wektor pr stron
        #print len(self.mesh.boundaries[self.id])
        for i in range(len(self.mesh.boundaries[self.id])):
            id_edge = self.mesh.boundaries[self.id, i]   # indeks krawedzi w WB
            c = self.field.mesh.list_kr[id_edge, 2]   # indeks wlasciciela do niego dopicac w rhs
            nr_kr_1, nr_kr_2 = self.field.mesh.list_kr[id_edge, 0], self.field.mesh.list_kr[id_edge, 1]
            wekt_ws = self.field.mesh.wsp_wekt_z_wsp(self.field.mesh.xy[nr_kr_1, :], self.field.mesh.xy[nr_kr_2, :])
            Rhs[c] += - self.deriv * self.field.mesh.dl_wekt(wekt_ws[0], wekt_ws[1])          #  dodac razy dlugosc

    def extrapolate(self):       # i to numer krawedzi w WB
        for i in range(self.mesh.boundaries[0, :].size):

            id_edge = self.mesh.boundaries[self.id, i]  # indeks krawedzi w WB
            field = self.field
            c = self.mesh.list_kr[id_edge, 2]  # pobierz wlascicela tej krawedzi
            cc1 = sum(field.mesh.xy[field.mesh.cells[c], :]) / len(field.mesh.cells[c])  # srodek komorki  [x,y]
            # print "cc1:", cc1
            #print self.mesh.boundaries[self.id, i]
            #print i

            n = self.mesh.edge_normal(self.mesh.boundaries[self.id, i])        # normalny do wektora (krawedzi WB) ale o jego dlugosci (trzeba pobrac numer pod krawedzi )
            n = n / np.linalg.norm(n)           # wersor normalny (podzielony przez swoja dlugosc)

            sr_pod_krawBC = ( (self.field.mesh.xy[self.mesh.boundaries_points[self.id, i, 0]] + self.field.mesh.xy[self.mesh.boundaries_points[self.id, i, 1]]) )/2  # wsp srodka podkrawedzi WB [x,y]
            # print "sr_pod_kr:", sr_pod_krawBC

            e = self.mesh.wsp_wekt_z_wsp(cc1, sr_pod_krawBC)

            #print "n: ", n
            #print "e: ", e
            ds = np.dot(-n, e)    # wektor od srodka komurki do srodka sciany
            #print "ds: ", ds
            #print " "

            flux = self.deriv
            Tsr = self.sol[c]
            print "Tsr:", Tsr
            #Tbrzeg = Tsr + flux*ds
            self.data[i] = Tsr + flux*ds
            print self.data[i]






class symmetry(Neuman):             #  # klasa dla kazdej krawedzi o warunku neuman o wartosci poch w kier normalym = 0
    def __init__(self, mesh, bId):
        Neuman.__init__(self, mesh, bId, np.array([0.] * len(mesh.boundaries[bId])))    # konstruktor pierwotnej klasy po ktorej dziedziczy, derrevativeValue przyjmuje teraz liste o wartosicach 0 i rozmiarze dlugosci mesh.boundaries
        self.id = bId        # zapisuje numer krawedzi pod oznaczeniem id
        self.mesh = mesh






class SurfField:              # to jest po prostu field z wartosciami rozwiazania dziedziczy po np.array
    def __init__(self, mesh):           # pobiera mesh a z nim jego rozmiar i boundaries czyli WB
        self.data = np.array([0.]*mesh.n)        # [0.]*mesh.n lista zer o rozmiarze mesh.n
        self.boundaries = [EdgeField(mesh, i) for i, _ in enumerate(mesh.boundaries)]          # (***) lista z pustymi warunkami brzegowymi dodatkowa do przechowania _ zmienna ktorej nikt nie uzyje interesuje nas tylko ilosc nie wartosc
        self.mesh = mesh
    # @property
    # def values(self):
    #     return self._sol

    def __setitem__(self, key, value):
        return self.data.__setitem__(key, value)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def setBoundaryCondition(self, bcObject):           # to zada warunki brzegowe bcObject( np class/object nemuam) na surface field
        bcObject.setField(self)                 # na przyklad  neuman.setFirld(self)
        self.boundaries[bcObject.id] = bcObject   # przypisuje do tej listy (***) (tylko z numerami krawedzi WB) zadane przez urzytkownika WB

    def setValues(self, lista_wart):       # pobiera rozwiazanie ukladu rownan (pole temp) - uaktualnia wartosci pola temp tak aby na krawedziach WB neumana byly temp nie strumienie ciepla (do wizualizacji)
        self.data = lista_wart          # rozwiazanie tu tablica "A"
        for b in self.boundaries:       # b kolejny warunek brzegowy tu 4 warunki
            b.upadate(self.data)        # przekaz pole temp do funkcji upadate


    def apply_bc(self, EqMat, Rhs):          # ma wrzucic ten warunek na macierz K i WPS
        for b in self.boundaries:   # petla po 4 war brzeg (te brzegi juz zapisalem uzywajac setBoundaryCondition
            b.manipulateMatrixVector(EqMat, Rhs)            # jesli neuman to wywola z klasy neuman jesil dirichlet to z dirichlet zalezy czym jest b




#
#
# def values(obiekt):
#     return obiekt._sol






