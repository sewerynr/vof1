import numpy as np


class Mesh:
    """obiekt_operujcy_na_danych_punktach_siatki"""
    def __init__(self, nodes, cells, boundaries):
        self.xy = nodes
        self.cells = cells
        self.n = cells.shape[0]    # shape returns the dimension of array. shape[0] daje ilosc wierszy shape[1] kolumn

        # self.__wsp_dl__ = None  # tez do (*****)
        self.wsp_dl = self.__wektor_wsp_idl__()         # wywolaj __prywatna metode__ i oddaj tablice wsp.dl

        self.list_kr = self.__lista_krawedzi__()
        self.list_kr = self.__wlasciciel_sasiad__(self.list_kr)
        self.boundaries_points = boundaries
        self.boundaries = self.__bound_to_edge_bound__(boundaries)

        # self.boundaries = boundaries

    #gdy chcemy wywolac prywatna funkcje ale dopiero gdy ktos o nia zapyta (dostanie wtedy macierz) (*****)
    # @property
    # def wsp_dl(self):
    #     if not self.__wsp_dl__:
    #         self.__wsp_dl__ = self.__wektor_wsp_idl__()
    #     return self.__wsp_dl__


    def __bound_to_edge_bound__(self, boundaries):
        bond_to_edge = np.array( [ [0]*len(b) for b in boundaries ] )

        for idBoundry, b in enumerate(boundaries):        # petla po 4 War Brzeg
            for idBEdge, edge_nodes in enumerate(b):             # petla po konkretnym WB
                for idEdge, edge in enumerate(self.list_kr):
                    if (edge_nodes[0] == edge[0] and edge_nodes[1] == edge[1]) or (edge_nodes[1] == edge[0] and edge_nodes[0] == edge[1]):
                        bond_to_edge[idBoundry][idBEdge] = idEdge
        #print bond_to_edge
        return bond_to_edge


    def __wektor_wsp_idl__(self):

        wekt_wsp_dlw = np.array([[0.] * 14] * (self.n))
        for l in range(self.n):

            wekt_wsp_dlw[l, 0] = self.xy[self.cells[l, 1], 0] - self.xy[self.cells[l, 0], 0]  # po x
            wekt_wsp_dlw[l, 1] = self.xy[self.cells[l, 1], 1] - self.xy[self.cells[l, 0], 1]  # po y

            wekt_wsp_dlw[l, 2] = self.xy[self.cells[l, 2], 0] - self.xy[self.cells[l, 1], 0]
            wekt_wsp_dlw[l, 3] = self.xy[self.cells[l, 2], 1] - self.xy[self.cells[l, 1], 1]

            wekt_wsp_dlw[l, 4] = self.xy[self.cells[l, 3], 0] - self.xy[self.cells[l, 2], 0]
            wekt_wsp_dlw[l, 5] = self.xy[self.cells[l, 3], 1] - self.xy[self.cells[l, 2], 1]

            wekt_wsp_dlw[l, 6] = self.xy[self.cells[l, 0], 0] - self.xy[self.cells[l, 3], 0]
            wekt_wsp_dlw[l, 7] = self.xy[self.cells[l, 0], 1] - self.xy[self.cells[l, 3], 1]

            # wpisuje dl wektorow

            wekt_wsp_dlw[l, 8] = self.dl_wekt(wekt_wsp_dlw[l, 0], wekt_wsp_dlw[l, 1])
            wekt_wsp_dlw[l, 9] = self.dl_wekt(wekt_wsp_dlw[l, 2], wekt_wsp_dlw[l, 3])
            wekt_wsp_dlw[l, 10] = self.dl_wekt(wekt_wsp_dlw[l, 4], wekt_wsp_dlw[l, 5])
            wekt_wsp_dlw[l, 11] = self.dl_wekt(wekt_wsp_dlw[l, 6], wekt_wsp_dlw[l, 7])

            # wpisuje sr komorek pobieram x z wsp_wezlow wczytujac jakie wezly tworza komorke z cells

            wekt_wsp_dlw[l, 12] = self.sr_komX(self.xy[self.cells[l, 1], 0], self.xy[self.cells[l, 0], 0], self.xy[self.cells[l, 3], 0], self.xy[self.cells[l, 2], 0])
            wekt_wsp_dlw[l, 13] = self.sr_komX(self.xy[self.cells[l, 1], 1], self.xy[self.cells[l, 0], 1], self.xy[self.cells[l, 3], 1], self.xy[self.cells[l, 2], 1])

        return wekt_wsp_dlw

    def edge_vector(self, edgeId):
        p1 = self.xy[self.list_kr[edgeId, 0], :]
        p2 = self.xy[self.list_kr[edgeId, 1], :]
        dp = p2 - p1
        return [dp[1], -dp[0]]

    def wsp_wekt(self, i, j):   # zamiast wsp pobiera numery poczatkowego i koncowego punktu wektora ( z listy wsp pkt xy)
        xwekt = self.xy[j, 0] - self.xy[i, 0]
        ywekt = self.xy[j, 1] - self.xy[i, 1]
        return xwekt, ywekt

    def wsp_wekt_z_wsp(self, w1, w2):  # pobiera dwa pkt w1[x1,y1],  w2[x2,y2]
        xwekt = w1[0] - w2[0]       #roznica wsp x z pktow w1 i w2
        ywekt = w1[1] - w2[1]
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
        yn = - x
        dlw = (x**2 + y**2)**0.5
        nx = xn / dlw
        ny = yn / dlw
        return nx, ny

    def edge_normal(self, i):           # i indeks krawedzi
        def_edge = self.list_kr[i, :2]          # wyciagam numery (indeksy) punktow z ktorch sklada sie dana krawedz  dwa pierwsze elementy wiersza
        #print i
        #print self.xy[def_edge[1]], self.xy[def_edge[0]]
        wektor = self.xy[def_edge[1]] - self.xy[def_edge[0]]

        return self.wektor_norm(*wektor)           # pierwszy z wektor to pierwsza zmienna itd

    def __lista_krawedzi__(self):

        s = 0
        for cell in self.cells:   # [ 1 0 12 11]  = [kr1 kr2 kr3 kr4]
           s += len(cell)       #zlicza ile w kolejnych komurkach pkt (krawedzi) (liczy ile jest w wierszu [1 2 12 11] czyli 4 krawedzie

        lista_kr = np.array([[0] * 4] * s)

        for i, cell in enumerate(self.cells):

            ilosc = lista_kr.shape[1]    # ilosc kolumn ( w cells ilosc kolumn odp ilosci krawedzi )
            for licz in range(0, ilosc, 1):

                if licz < ilosc -1:
                    lista_kr[i * ilosc + licz] = [cell[licz], cell[licz + 1], i, -1]
                elif licz == ilosc -1:
                    lista_kr[i * ilosc + licz] = [cell[licz], cell[licz - (ilosc - 1)], i, -1]
        return lista_kr     # [ 1 0 0 1]  = [pkt1 pkt2 wl sasiad]

    def __wlasciciel_sasiad__(self, lista_krawedzi):
        # do k_ids przypisuje dwie kolumny ze wszystkich wierszy czyli numery wezlow tworzacych krawedzie
        k_ids = lista_krawedzi[:, :2]

        # tworze liste o nazwie pairs (do listy() mozna dopisywac elementy przez pairs.append() )
        pairs = list()

        for i1, k in enumerate(k_ids):
            for i2, j in enumerate(k_ids):
                if k[0] == j[1] and k[1] == j[0]:
                    pairs.append((i1, i2))

        do_wyrzucenia = list()
        count = 0

        # zapisuje w lista_krawedzi sasiada (biorac wlascicela z powtarzajacej sie jako sasiada we wczesjniejszej) oraz wpisuje powtarzajace sie do_wyrzucenia
        for p in pairs:
            if p[0] not in do_wyrzucenia:
                # wlas podobnej kraw. staje sie sasiadem w krawedzi procesowanej
                lista_krawedzi[p[0], 3] = lista_krawedzi[p[1], 2]
                do_wyrzucenia.append(p[1])

        # wybiera tak aby sie nie powtarzaly
        # sprawdza czy numer krawedzi z lista_kr jest w do_wyrzucenia jesli nie to dopisuje do uniklnych i nadpisuje stara tablice lista_kr nowa bez powtarzajacych sie

        list_kr_unik = list()
        for i, krawed in enumerate(lista_krawedzi):
            if i not in do_wyrzucenia:
                list_kr_unik.append(krawed)

        # tu nadpisuje nowa
        # lista_kr_all = lista_krawedzi
        lista_krawedzi = np.array(list_kr_unik)

        del list_kr_unik, do_wyrzucenia
        return lista_krawedzi