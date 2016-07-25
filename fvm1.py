import numpy as np
import matplotlib.pyplot as plt
from mesh import *
#from functions import *

def siatka_regularna_prost(n, dx, dy, x0, y0):
    node_coordinates = np.array([[0.] * 2] * (n + 1) ** 2)
    for i in range(0, n + 1, 1):
        for j in range(0, n + 1, 1):
            index = i * (n + 1) + j
            node_coordinates[index, :] = [x0 + j * dx, y0 + i * dy]

    cells = np.array([[0] * 4] * (n ** 2))          # kazda komurka ma 4 pkt a komorek jest n^2
    for i in range(n):
        for j in range(n):
            cid = i * n + j
            cells[cid, :] = [i * (n + 1) + j, i * (n + 1) + j + 1, (i + 1) * (n + 1) + 1 + j, (i + 1) * (n + 1) + j]       # [wpisuje w wiersz o nr cid, : czyli  w kazda kolumne] odp w miejsce 1,2,3,4



    boundaries_BC = np.array([[[0]*2]*n]*4)       # jedna kolumna i 4 wiersze bo cztery brzegi (w kolumnie nowa lista tyle wierszy ile krawedzi i tyle kolumn ile pkt tworzacych krawedz
    temp = np.array([[[0] * 1] * (n+1)])

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ??????????????????????????  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[1] == 0:
            temp[0, licznik] = id
            licznik = licznik + 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[0, i, 0] = temp[0, i]
        boundaries_BC[0, i, 1] = temp[0, i + 1]

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[0] == 1:
            temp[0, licznik] = id
            licznik = licznik + 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[1, i, 0] = temp[0, i]
        boundaries_BC[1, i, 1] = temp[0, i + 1]

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[1] == 1:
            temp[0, licznik] = id
            licznik = licznik + 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[2, i, 0] = temp[0, i]
        boundaries_BC[2, i, 1] = temp[0, i + 1]

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[0] == 0:
            temp[0, licznik] = id
            licznik = licznik + 1
    for i in range(0, licznik - 1, 1):
       boundaries_BC[3, i, 0] = temp[0, i]
       boundaries_BC[3, i, 1] = temp[0, i + 1]

    #print boundaries_BC   # krawedzie przez numery punktow [0 1] [1 2] itd...

    return node_coordinates, cells, boundaries_BC
    # Utworzenie komorek na podstawie wezlow (zbior nr wezlow budujacych kom.).
    # Definicja 1 komorki: [ nr wz1, nr wz2, nr wz3, nr wz4 ]

def npArrayMatrix(mesh):
    return np.array([[0.] * mesh.n] * mesh.n)

from fvMatrix import fvMatrix

# ta do ustalonej w czasie
def laplace(coeff, field, matrixGeneratorFunction = fvMatrix):    # matrixProvider wywoluje konstruktor metody fvMatrix  z polem T lambda funkcja anonimowa
    n, lista_kra = field.mesh.n, field.mesh.list_kr
    macierz_K_e = matrixGeneratorFunction(field.mesh)
    mesh = field.mesh

    from field import EdgeField, SurfField, Neuman
    coeffField = SurfField(mesh, bcGenerator=Neuman)

    if not hasattr(coeff, "__iter__"):                  # czy to liczba czy tablica
        coeff = np.ones(len(mesh.cells)) * coeff

    coeffField.setValues(np.array(coeff))

    edgeCoeff = EdgeField.interp(coeffField)


    for i, kraw in enumerate(lista_kra):
        if kraw[3] > -1:
            k1, k2, c, f = kraw
            cC = mesh.cell_centers[c]
            cF = mesh.cell_centers[f]

            CF = cF - cC
            Snorm = mesh.Se[i]
            a = edgeCoeff.data[i] * np.dot(CF, Snorm) / np.dot(CF, CF)

            # Wstawiamy wsp. do ukladu rownan
            macierz_K_e[c, c] += - a                         # to co odp wlascicielowi, diagonalny element
            macierz_K_e[c, f] += a                           # to co odp sasiadowi z przeciwnym znakiem, poza diagonala
            # kazda krawedz tylko raz ale strumien na niej jest bilansowany w dwoch sasiadujacych komorkach, dlatego
            # takie same wsp (nawet znak, bo kierunek normalnej oraz kierunek od srodka do drugiej sie razem obracaja)
            # wiec dokladnie te same wspol. ida do rownania o indeksie takim jak sasiad
            macierz_K_e[f, f] += - a # krawedz wplywa na rownanie sasiada, teraz sasiad jest w centrum komorki
            macierz_K_e[f, c] += a # a poza diagonala jest wlasciciel

    rhs = np.zeros(n)

    field.apply_bc_diffusiveFlux(edgeCoeff, macierz_K_e, rhs)

    return macierz_K_e, rhs



    #   div to dywergencja, czlon zachowawczy adwekcji
def div(phi, field, matrixGeneratorFunction = fvMatrix):                  # phi to pole predkosci na scianach skalarne bo przemnozone skalarnie razy wektor normalny (ro * wektor predkosci * wekt normalny) = phi
    n, lista_kra = field.mesh.n, field.mesh.list_kr                      # lista kr: [ 1 0 0 1]  = [pkt1 pkt2 wl sasiad]
    D = matrixGeneratorFunction(field.mesh)

    Rhs = np.zeros(n)
    mesh = field.mesh

    for i, k in enumerate(mesh.list_kr):
        w, s = k[2:]
        edgeLen = mesh.eLengths[i]
        phiEdge = phi.data[i]                                # pobiera wartosci predkosci z macierzy phi[dla elementu i]

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wiersz mowi ktora komorka kolumna co i skad wlata wylata  (strumien o jakiejs temp)   !!!!!!!!!!!!!!!!!!!!!!!!!!

        # Upiwnd
        if s > -1:
            if phiEdge > 0:                              # od wlasiciela do sasiada
                D[w, w] += phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [od wl , temp wl]
                D[s, w] -= phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [do sas, temp wl]
            else:                                        # phiedge < 0 mniejsze od sasi ada do wlasciciela
                D[s, s] -= phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [od sasiada , z temp sasiada]
                D[w, s] += phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [do wl , temp sasiada]

        #Central
        # if k[3] > -1:
        #         D[w, w] -= phiEdge * edgeLen / 2
        #         D[w, s] -= phiEdge * edgeLen / 2
        #         D[s, s] += phiEdge * edgeLen / 2
        #         D[s, w] += phiEdge * edgeLen / 2


    field.apply_bc_convectiveFlux(D, Rhs, phi.data)
    return D, Rhs


def eInt(edgeField):
    import numpy as np

    mesh = edgeField.mesh
    res = np.zeros(len(mesh.cell_centers))

    for i, (v, k) in enumerate(zip(edgeField.data, mesh.list_kr)):
        eLen = mesh.eLengths[i]
        res[k[2]] += v * eLen
        if k[3] > -1:
            res[k[3]] -= v * eLen

    return res


def eInt_implicit(mesh):
    pass


def adjustPhi(phiEdgeField):
    pass






def WB_dir_T1(lista_kr, Td, wsp_wezl, macierz_K_e, rhs):
    # warunki brzegowe 1) jako T w centrum komorki wiec po prostu w macierz K wpisuje w miejsce odp sr. kom 1 a w wektor pr stron = danej temp
    for kraw in lista_kr:
        if kraw[3] == -1:
            c = kraw[2]                       # indeks komorki wlascicela
            macierz_K_e[c, :] = 0             # wakazuje na wiersz czyli na komorke i zapisuje we wszystkich elementach wiersza 0
            macierz_K_e[c, c] = 1
            if wsp_wezl[kraw[0], 1] == 0 and wsp_wezl[kraw[1], 1] == 0:          # spr wspl obu wezlow krawedzi czy sa rowne zero po y
                rhs[c] = Td

    return rhs


# wydruk konturow temperatur
# Wyswietl rozw. w wezlach interpolujemy otrzymane wyniki na srodki komorek PO INTERPOLACJI

from interpolacja import inter

def draw_values_edges(wsp_wezl, cells, listak, T, n, Tdirich, DlPrzX, DlPrzY):
    '''
    :param  T jako pole temperatur z rozwiazania Tdirich WB dir
    '''
    Tinter = inter(wsp_wezl, cells, T.data)

    for bId, b in enumerate(T.mesh.boundaries):
        varB = T.boundaries[bId]
        for eLocal, eGlobal in enumerate(b):
            valOnEdge = varB[eLocal]
            node1 = T.mesh.list_kr[eGlobal, 0]
            node2 = T.mesh.list_kr[eGlobal, 1]
            Tinter[[node1, node2]] = 0


    for bId, b in enumerate(T.mesh.boundaries):
        varB = T.boundaries[bId]
        for eLocal, eGlobal in enumerate(b):
            valOnEdge = varB[eLocal]
            node1 = T.mesh.list_kr[eGlobal, 0]
            node2 = T.mesh.list_kr[eGlobal, 1]
            Tinter[ [node1, node2] ] += valOnEdge/2


    # numpy linespace dzieli przedzial na n podzialow o rownym odstepie
    X, Y = np.meshgrid(np.linspace(0, DlPrzX, n+1), np.linspace(0, DlPrzY, n+1))
    T_new = Tinter.reshape((n+1, n+1))
    plt.figure()
    cont = plt.contourf(X, Y, T_new, 14)
    plt.colorbar(cont)

    draw_edges(wsp_wezl, listak)

    plt.show()


# Wyswietl wart. w srodkach
# numpy linespace dzieli przedzial na n podzialow o rownym odstepie / T pole temp po rozwiazaniu
def draw_values_centers(dx, dy,  DlPrzX, DlPrzY, n, T):
    X, Y = np.meshgrid(np.linspace(dx/2, DlPrzX-dx/2, n), np.linspace(dy/2, DlPrzY-dy/2, n))
    T_new = T.reshape((n, n))
    plt.figure()
    cont = plt.contourf(X, Y, T_new)
    plt.colorbar(cont)
    plt.show()


#po indeksach rysuje komorke

def index_draw_cell(cells, wsp_wezl):
    plt.figure()                                    # inicjalizuje okno
    for id, cell in enumerate(cells):               #pobiera numery wierszy macierzy "cell" i leci po nich
        plt.plot(wsp_wezl[cell, 0], wsp_wezl[cell, 1], 'b-')        # rysuje linie z x i y
    plt.gca().set_xlim([-0.1, 1.1])                 # granice rysowania od -0.1 do 1.1
    plt.gca().set_ylim([-0.1, 1.1])
    plt.show()                                      # dopoki sie tego nie napisze dodaje nastepne rysownia do tego samego okna


# rysuje krawedzie siatki
def draw_edges(wsp_wezl, lista_kr):

    #plt.figure()
    for k in lista_kr:
        p1 = wsp_wezl[k[0]]
        p2 = wsp_wezl[k[1]]
        if k[3] == -1:
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'r-')
        else:
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color="black")


    # plt.gca().set_xlim([-0.1, 1.1])           # granice rysowania od -0.1 do 1.1
    # plt.gca().set_ylim([-0.1, 1.1])
    # plt.show()


def animate_contour_plot(framesDatas, sizeX=(0, 1), sizeY=(0, 1), dataRange=None, nLevels=10, skip=1, repeat=False, interval=5):
    """
    Function which make animation from set of 2D data on cartesian grid
    :param framesDatas: List of 2D numpy.arrays containing nodal values
    :param sizeX: tuple holding domain range in X dir
    :param sizeY: tuple holding domain range in Y dir
    :param skip: number of frames to be skipped
    :param nLevels: number of color levels
    :param dataRange: tuple holding min and max value for data range to be displayed in plot and colorbar
    :return: None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import animation

    if len(framesDatas) == 0:
        raise Exception("Data frames number should be at least one")

    Nx, Ny = framesDatas[0].shape

    X, Y = np.meshgrid(np.linspace(sizeX[0], sizeX[1], Nx), np.linspace(sizeY[0], sizeY[1], Ny))

    if not dataRange:
        minD = min(framesDatas[0].flatten())
        maxD = max(framesDatas[0].flatten())
    else:
        minD, maxD = dataRange

    fig = plt.figure()
    plt.axes().set_aspect('equal', 'datalim')
    ticks = np.linspace(minD, maxD, nLevels + 1)
    cs = plt.contourf(X, Y, framesDatas[0], ticks)
    cbar = fig.colorbar(cs, ticks=ticks)
    cbar.ax.set_yticklabels(map(str, ticks))

    if len(framesDatas) > 1:
        def animate(i):
            i = i * skip
            cs = plt.contourf(X, Y, framesDatas[i], ticks)
            # cbar.ax.set_yticklabels(map(str, ticks))
            # cbar.update_ticks()

            cs.zmin = minD
            cs.zmmax = maxD
            plt.title('Frame %d' % (i + 1))
            return cs

        return animation.FuncAnimation(fig, animate, frames=len(framesDatas) / skip, interval=interval, repeat=repeat)
# for id, cell in enumerate(k_ids):
#     print id, cell

# print np.allclose(np.dot(macierz_K_e, T), wektor_pr_str)
# print T[:n], len(T)
# print macierz_K_e[:n,:]
# print wektor_pr_str[:n]
# print wektor_pr_str

def source(mesh, cellData):
    import numpy as np
    return np.multiply(mesh.cells_areas, cellData)


