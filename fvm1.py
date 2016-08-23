import numpy as np
import matplotlib.pyplot as plt
from mesh import *
#from functions import *

def siatka_regularna_prost(n, m,  dx, dy, x0, y0):
    node_coordinates = np.array([[0.] * 2] * (n + 1) * (m + 1))
    for i in range(0, n + 1, 1):
        for j in range(0, m + 1, 1):
            index = i * (m + 1) + j
            node_coordinates[index, :] = [x0 + j * dx, y0 + i * dy]

    cells = np.array([[0] * 4] * n * m)

    for i in range(n):
        for j in range(m):
            cid = i * m + j
            cells[cid, :] = [i * (m + 1) + j, i * (m + 1) + j + 1, (i + 1) * (m + 1) + 1 + j, (i + 1) * (m + 1) + j]
            # [wpisuje w wiersz o nr cid, : czyli  w kazda kolumne] odp w miejsce 1,2,3,4
    # print "s", cells.shape
    boundaries_BC = list()
    mat = list()
    for i in range(m):
        mat.append([0, 0])

    boundaries_BC.append(mat)

    mat = list()
    for i in range(n):
        mat.append([0, 0])

    boundaries_BC.append(mat)

    mat = list()
    for i in range(m):
        mat.append([0, 0])

    boundaries_BC.append(mat)

    mat = list()
    for i in range(n):
        mat.append([0, 0])

    boundaries_BC.append(mat)
    # print boundaries_BC

    temp = np.array([[[0] * 1] * (max(m, n)+1)])

    # boundaries_BC[2][[1][0]] = [100, 78]
    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[1] == 0:
            temp[0, licznik] = id
            licznik += 1

    for i in range(0, licznik - 1, 1):
        boundaries_BC[0][i][0] = int(temp[0, i])
        boundaries_BC[0][i][1] = int(temp[0, i + 1])

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[0] == 1:
            temp[0, licznik] = id
            licznik += 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[1][i][0] = int(temp[0, i])
        boundaries_BC[1][i][1] = int(temp[0, i + 1])

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[1] == 1:
            temp[0, licznik] = id
            licznik += 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[2][i][0] = int(temp[0, i])
        boundaries_BC[2][i][1] = int(temp[0, i + 1])

    licznik = 0
    for id, wspw in enumerate(node_coordinates):
        if wspw[0] == 0:
            temp[0, licznik] = id
            licznik += 1
    for i in range(0, licznik - 1, 1):
        boundaries_BC[3][i][0] = int(temp[0, i])
        boundaries_BC[3][i][1] = int(temp[0, i + 1])

    # print boundaries_BC   # krawedzie przez numery punktow [0 1] [1 2] itd...

    return node_coordinates, cells, boundaries_BC

def npArrayMatrix(mesh):
    return np.array([[0.] * mesh.n] * mesh.n)

from fvMatrix import fvMatrix


def laplace(coeff, field, matrixGeneratorFunction=fvMatrix):         #matrixGeneratorFunction wywoluje konstruktor metody fvMatrix  z polem T
    n, lista_kra = field.mesh.n, field.mesh.list_kr
    macierz_K_e = matrixGeneratorFunction(field.mesh)
    mesh = field.mesh

    from field import EdgeField, SurfField, Neuman

    if not hasattr(coeff, "__iter__") and not isinstance(coeff, EdgeField):          # czy to liczba czy tablica o wym. [len(mesh.cells) x 2]
        coeff = np.ones((len(mesh.cells), 2), dtype=float) * coeff                   # liczba wiec stworz wektor i pomnoz coeff razy wektor jedynek [1,1]
    elif isinstance(coeff, np.ndarray) and (len(coeff.shape) == 1 or coeff.shape[1] == 1):          # tablica wiec tylko ja obroc
        coeff = np.array([coeff, coeff]).T

    # another but different if
    if isinstance(coeff, EdgeField):                # gdy wspolczynniki to pole edgefield ( na kazdej krawedzi inny wsp)
        edgeCoeff = coeff
    else:
        coeffFieldX = SurfField(mesh, bcGenerator=Neuman)
        coeffFieldY = SurfField(mesh, bcGenerator=Neuman)
        coeffFieldX.setValues(np.array(coeff[:, 0]))
        coeffFieldY.setValues(np.array(coeff[:, 1]))
        edgeCoeff = EdgeField.vector(EdgeField.interp(coeffFieldX), EdgeField.interp(coeffFieldY))

    for i, kraw in enumerate(lista_kra):
        if kraw[3] > -1:
            k1, k2, c, f = kraw
            CF = mesh.center_to_center_edge[i]
            Snorm = mesh.Se[i]
            coeff = edgeCoeff.data[i]
            a = (coeff * CF).dot(Snorm) / CF.dot(CF)
            macierz_K_e[c, c] += - a           # to co odp wlascicielowi, diagonalny element
            macierz_K_e[c, f] += a             # to co odp sasiadowi z przeciwnym znakiem, poza diagonala
            macierz_K_e[f, f] += - a           # krawedz wplywa na rownanie sasiada, teraz sasiad jest w centrum komorki
            macierz_K_e[f, c] += a             # a poza diagonala jest wlasciciel

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
        if s > -1:
            if phiEdge > 0:                            # od wlasiciela do sasiada, skad wylata tam dodaje do macierzy
                D[w, w] += phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [od wl , temp wl]
                D[s, w] -= phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [do sas, temp wl]
            else:                                      # phiedge < 0 mniejsze od sasi ada do wlasciciela
                D[s, s] -= phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [od sasiada , z temp sasiada]
                D[w, s] += phiEdge * edgeLen           # [ skad wylata/dokad , z jaka temp ] => [do wl , temp sasiada]

        # Pe = U*dx / a_diff
        # if k[3] > -1:
        #         D[w, w] += phiEdge * edgeLen / 2      # dodatnie phi to wylata od wlasciciela wiec do niego dodac
        #         D[w, s] += phiEdge * edgeLen / 2
        #         D[s, s] -= phiEdge * edgeLen / 2      # phi ujemne to od sasaiada - i - da plus dodamy do sadiada
        #         D[s, w] -= phiEdge * edgeLen / 2

    field.apply_bc_convectiveFlux(D, Rhs, phi.data)
    return D, Rhs


def eInt_implicit(mesh, matrixGen = lambda dims : np.zeros(shape=dims)):        # P
    P = matrixGen( (mesh.n, len(mesh.list_kr)) )      #jesli wywolujac eInt_implicit nie podamy to wywola lambda
    # wlacznie z WB dla WB -eLen
    for i, (k) in enumerate(mesh.list_kr):
        eLen = mesh.eLengths[i]
        P[k[2], i] = eLen
        if k[3] > -1:
            P[k[3], i] = - eLen
    return P

# zczytaj indeksy krawedzi nie brzegowych do listy index
def adjustPhi_eqSys(phiEdgeField):
    P = eInt_implicit(phiEdgeField.mesh)
    index = []
    for i, k in enumerate(phiEdgeField.mesh.list_kr):
        if k[3] != -1:
            index.append(i)             # indeksy krawedzi nie brzegowych (nie WB)
    return P, index

def adjustPhi(phiEdgeField):
    P, internIndex = adjustPhi_eqSys(phiEdgeField)
    Pp = P[:, internIndex]      # Pp to P bez warunkow brzegowych ich nie chcemy poprawiac na nich ma byc to co zadane
    F = P.dot(phiEdgeField.data)
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import cg

    M = Pp.dot(Pp.T)
    M = csr_matrix(M)
    Lambda = cg(M, F)[0]
    dPhi = - Pp.T.dot(Lambda)
    phiEdgeField.data[internIndex] += dPhi
    return 1

def adjustPhiold(phiEdgeField, matrixGen = lambda dims : np.zeros(shape=dims)):
    P = matrixGen((phiEdgeField.mesh.n, len(phiEdgeField.mesh.list_kr)))
    for i, (k) in enumerate(phiEdgeField.mesh.list_kr):
        eLen = phiEdgeField.mesh.eLengths[i]
        P[k[2], i] = eLen
        if k[3] > -1:
            P[k[3], i] = - eLen
    F = P.dot(phiEdgeField.data)
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import cg

    M = P.dot(P.T)
    M = csr_matrix(M)
    Lambda = cg(M, F)[0]
    dPhi = - P.T.dot(Lambda)
    phiEdgeField.data += dPhi
    return 2

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


def animate_contour_plot(framesDatas, sizeX=(0, 1), sizeY=(0, 1), dataRange=None, nLevels=10, skip=1, repeat=False, interval=5, diff=1., dt=1., nN=1., Tempa=-0., adj=0.):
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

    if len(framesDatas) == 0:   # frames data to T.data.reshape(n, m)
        raise Exception("Data frames number should be at least one")

    Ny, Nx = framesDatas[0].shape
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
    Aa = str(diff)
    Ba = str(dt)
    Ca = str(nN)
    Ta = str(Tempa)
    adj = str(adj)
    Aa = "diff: " + Aa + "   dt: "+ Ba + "   n: " + Ca + "   AdjPhi: " + adj
    cbar.ax.set_ylabel(Aa)

    if len(framesDatas) > 1:
        def animate(i):
            i = i * skip
            cs = plt.contourf(X, Y, framesDatas[i], ticks)
            cs.zmin = minD
            cs.zmmax = maxD
            plt.title('Frame %d' % (i + 1))
            return cs

        return animation.FuncAnimation(fig, animate, frames=len(framesDatas) / skip, interval=interval, repeat=repeat)

def source(mesh, cellData):
    import numpy as np
    return np.multiply(mesh.cells_areas, cellData)


