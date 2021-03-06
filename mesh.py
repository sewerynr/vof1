from duplicity.pexpect import which
import numpy as np
import time

def iterExtra(tab):     #  bierze wiersz z cells np. [0 1 11 12]  i dodaje do niego na koncu pierwszy element [0 1 11 12 0]
    for i in tab:
        yield i
    yield tab[0]


class Mesh:

    """obiekt_operujcy_na_danych_punktach_siatki"""
    def __init__(self, nodes, cells, boundaries):
        self.xy = nodes
        self.cells = cells
        self.n = cells.shape[0]                         # shape returns the dimension of array. shape[0] daje ilosc wierszy shape[1] kolumn
        self.cells_areas = self.__cell_area__()
        self.list_kr = self.__lista_krawedzi__()
        self.list_kr = self.__wlasciciel_sasiad__(self.list_kr)
        self.edgeCenters = self.__edgeCenters__()
        self.boundaries_points = boundaries
        self.boundaries = self.__bound_to_edge_bound__(boundaries)
        self.Se = self.__Se__()
        self.normals, self.eLengths = self.__normals_and_edge_lengths__()
        self.cell_centers = self.__cell_center__()
        self.center_to_center_edge = self.__center_to_center_edge__()

    #gdy chcemy wywolac prywatna funkcje ale dopiero gdy ktos o nia zapyta (dostanie wtedy macierz) (*****)
    # @property
    # def wsp_dl(self):
    #     if not self.__wsp_dl__:
    #         self.__wsp_dl__ = self.__wektor_wsp_idl__()
    #     return self.__wsp_dl__

    def __Se__(self):               #krawedz normalna do kr o dl krawedzi
        se = np.array([[0.] * 2] * len(self.list_kr))
        for i, edge in enumerate(self.list_kr):
            p1 = self.xy[edge[0], :]     # edge[0] mowi ktory wiesz czyli ktory pkt. : to wez x, y
            p2 = self.xy[edge[1], :]
            x, y = self.wsp_wekt_z_wsp(p2, p1)
            se[i] = [-y, x]
        return se

    def __normals_and_edge_lengths__(self):
        magSqrt = np.sqrt(np.sum(np.multiply(self.Se, self.Se), axis=1))
        return np.array([self.Se[:, 0] / magSqrt, self.Se[:, 1] / magSqrt]).T, magSqrt


    def __cell_center__(self):
        cc = np.array([[0.] * 2] * (self.n))  # cell centers _ srodki komorek
        for i, cell in enumerate(self.cells):
            c1, c2, c3, c4 = self.cells[i]
            dl = len(self.cells[i])
            cc[i][0] = (self.xy[c1][0] + self.xy[c2][0] + self.xy[c3][0] + self.xy[c4][0]) / dl  # srodki komomorek pobiera numery wezlow z cells i wczytuje wsp z wsp_wezl

            cc[i][1] = (self.xy[c1][1] + self.xy[c2][1] + self.xy[c3][1] + self.xy[c4][1]) / len(self.cells[i])  # srodki komomorek pobiera numery wezlow z cells i wczytuje wsp z wsp_wezl
        return cc


    def __cell_area__(self):
        area = np.zeros(self.n)
        for i in range(self.n):  # dla komorki,
            punkty = []
            for node in iterExtra(self.cells[i]):  # podaje kolejne punkty z czytac punkty
                pointk_xy = self.xy[node]  # joty punkt z komorki i
                punkty.append(pointk_xy[0])  # 0 2 4 6... wsp X
                punkty.append(pointk_xy[1])  # 1 3 5 7... wsp y
            area[i] = self.area_from_coordinates(punkty)
        return area


    def area_from_coordinates(self, coord_list):
        area = 0
        n = len(coord_list)
        for i in range(0, n - 2, 2):
            area += 0.5 * (coord_list[i + 2] + coord_list[i]) * (coord_list[i + 3] - coord_list[i + 1])
        return area


    # zapisuje 4 war brzegowe jako liste z numerami krawedzi
    def __bound_to_edge_bound__(self, boundaries):                      # boundaries zawiera 4 krawedzie z siatki reg prost
        bond_to_edge = np.array( [ [0]*len(b) for b in boundaries ] )
        for idBoundry, b in enumerate(boundaries):                      # petla po 4 War Brzeg bo tyle mamy w naszej siatce gdyby wiecej to po wiekszej liczbie elementow
            for idBEdge, edge_nodes in enumerate(b):                    # petla po konkretnym WB (jego krawedziach p1----p2)
                en0 = edge_nodes[0]
                en1 = edge_nodes[1]
                for idEdge, edge in enumerate(self.list_kr):            # przelec po wszystkich krawedziach (list_kr) i porownaj z tym co w BC
                    e0 = edge[0]
                    e1 = edge[1]
                    if (en0 == e0 and en1 == e1) or (en1 == e0 and en0 == e1):
                        bond_to_edge[idBoundry][idBEdge] = idEdge       # przypisz krawedzi jej numer jesli punky sobie odpowiadaja
        return bond_to_edge


    def edge_vector(self, edgeId):
        p1 = self.xy[self.list_kr[edgeId, 0], :]
        p2 = self.xy[self.list_kr[edgeId, 1], :]
        dp = p2 - p1
        return [dp[1], -dp[0]]

    def wsp_wekt(self, i, j):                    # zamiast wsp pobiera numery poczatkowego i koncowego punktu wektora ( z listy wsp pkt xy)
        xwekt = self.xy[j, 0] - self.xy[i, 0]
        ywekt = self.xy[j, 1] - self.xy[i, 1]
        return xwekt, ywekt

    def wsp_wekt_z_wsp(self, w1, w2):            # pobiera dwa pkt w1[x1,y1],  w2[x2,y2]
        xwekt = w2[0] - w1[0]                    # roznica wsp x z pktow w1 i w2
        ywekt = w2[1] - w1[1]
        return xwekt, ywekt

    def dl_wekt(self, xwekt, ywekt):
        dl_w = (xwekt**2 + ywekt**2)**0.5
        return dl_w

    def sr_komX(self, x1, x2, x3, x4):
        srkX = (x1+x2+x3+x4)/4
        return srkX

    def sr_komY(self, y1, y2, y3, y4):
        srkY = (y1+y2+y3+y4)/4
        return srkY


    def wektor_norm(self, x, y):
        xn = y
        yn = -x
        dlw = (x**2 + y**2)**0.5
        nx = xn / dlw
        ny = yn / dlw
        return np.array([nx, ny])

    def edge_normal(self, i):                            # i indeks krawedzi   liczy normalna do krawedzi zapisanej w liscie pod indeksem i
        def_edge = self.list_kr[i, :2]                   # wyciagam numery (indeksy) punktow z ktorch sklada sie dana krawedz  dwa pierwsze elementy wiersza
        wektor = self.xy[def_edge[1]] - self.xy[def_edge[0]]

        return self.wektor_norm(*wektor)                 # pierwszy z wektor to pierwsza zmienna itd

    def __lista_krawedzi__(self):
        s = 0
        for cell in self.cells:                          # [ 1 0 12 11]  = [kr1 kr2 kr3 kr4]
            s += len(cell)                                # zlicza ile w kolejnych komurkach pkt (krawedzi) (liczy ile jest w wierszu [1 2 12 11] czyli 4 krawedzie
        lista_kr = np.array([[0] * 4] * s)

        for i, cell in enumerate(self.cells):
            ilosc = lista_kr.shape[1]                   # ilosc kolumn ( w cells ilosc kolumn odp ilosci krawedzi )
            for licz in range(0, ilosc, 1):
                if licz < ilosc -1:
                    lista_kr[i * ilosc + licz] = [cell[licz], cell[licz + 1], i, -1]
                elif licz == ilosc -1:
                    lista_kr[i * ilosc + licz] = [cell[licz], cell[licz - (ilosc - 1)], i, -1]

        return lista_kr
        # [ 1 0 0 1]  = [pkt1 pkt2 wl sasiad]

    def __edgeCenters__(self):
        s = len(self.list_kr)
        edgeCenters = np.array([[0.] * 2] * s)

        for krid, kr in enumerate(self.list_kr):
            p1, p2 = kr[: 2]
            x1, y1 = self.xy[p1]
            x2, y2 = self.xy[p2]
            edgeCenters[krid] = (x1 + x2) / 2, (y1 + y2) / 2

        return  edgeCenters


    def __center_to_center_edge__(self):
        ecf = np.array([[0.] * 2] * len(self.list_kr))
        for i, edge in enumerate(self.list_kr):

            if edge[3] > -1:  # jesli sasiad to nie WB
                ccw = self.cell_centers[edge[2]]   # cell center wlasciciela
                ccs = self.cell_centers[edge[3]]
                ecf[i] = ccs - ccw

            else:   # pobierz sr kom i srodek scianki WB i policz
                ccw = self.cell_centers[edge[2]]  # cell center wlasciciela
                p1, p2 = edge[: 2]
                x1, y1 = self.xy[p1]
                x2, y2 = self.xy[p2]
                ccs = (x1 + x2) / 2, (y1 + y2) / 2
                ecf[i] = ccs - ccw
        return ecf


    def __wlasciciel_sasiad__(self, lista_krawedzi):
        k_ids = lista_krawedzi[:, :2]
        # tworze liste o nazwie pairs (do listy() mozna dopisywac elementy przez pairs.append() )
        pairs = list()

        for i1, k in enumerate(k_ids.tolist()):
            m1 = k[0] == k_ids[:, 1]
            m2 = k[1] == k_ids[:, 0]
            m = m1*m2
            id = np.argwhere(m).flatten()

            if id.shape[0] > 0:
                pairs.append((i1, int(id[0])))

        do_wyrzucenia = list()

        # zapisuje w lista_krawedzi sasiada (biorac wlascicela z powtarzajacej sie jako sasiada we wczesjniejszej) oraz wpisuje powtarzajace sie do_wyrzucenia
        for p in pairs:
            if p[0] not in do_wyrzucenia:
                # wlas podobnej kraw. staje sie sasiadem w krawedzi procesowanej
                lista_krawedzi[p[0], 3] = lista_krawedzi[p[1], 2]
                do_wyrzucenia.append(p[1])

        list_kr_unik = list()
        for i, krawed in enumerate(lista_krawedzi):
            if i not in do_wyrzucenia:
                list_kr_unik.append(krawed)

        lista_krawedzi = np.array(list_kr_unik)

        del list_kr_unik, do_wyrzucenia

        return lista_krawedzi

