# macierz dla obliczen pedu bedzie sie zmieniala wiec trzeba bedzie ja modyfikowac z kroku na krok.

from copy import deepcopy

class fvMatrix:                   # do sparse co na przek, jakie wartosci (data) w jakim miejscu (indeses - kolumy)
    def __init__(self, meshOrMatrix):
                # kopie list
        if isinstance(meshOrMatrix, fvMatrix):           # czy meshOrMatrix jest obiektem fvMatrix
            self.data = deepcopy(meshOrMatrix.data)
            self.indices = deepcopy(meshOrMatrix.indices)       # kopia list
            self.diag = deepcopy(meshOrMatrix.diag)
            self.N = meshOrMatrix.N
            self.M = meshOrMatrix.M
        else:
            if isinstance(meshOrMatrix, int):            # jezeli ktos poda int jako rozmiar macierzy
                self.N = meshOrMatrix
                self.M = self.N
            elif isinstance(meshOrMatrix, tuple):           # (liczba, liczba)
                self.N, self.M = meshOrMatrix
            else:                                        # jesli argument wejsciowy to mesh
                self.N = meshOrMatrix.n
                self.M = self.N

            self.data = [list() for i in range(self.N)]
            self.indices = [list() for i in range(self.N)]
            self.diag = [0. for i in range(self.N)]

        # domyslnie nie ma tych wartosci dopoki ich nie utworzy
        self.__cache__ = None                            # nie istnieje obiekt cache
        self.vectorData = None
        self.vectorIndices = None
        self.rowLengths = None

        self.dtype = float


    @property   # mozna uzyc jako zmiennej
    def sparse(self):
        if self.__cache__ is None:          # gdy nie istnieje to:
            self.optimize()                 # (***) wywolaj optimize() ktora wpisuje wartosci do list
            from scipy.sparse import csr_matrix  # zaimportuj macierze rzadkie
            self.__cache__ = csr_matrix((self.vectorData, self.vectorIndices, self.rowLengths), shape=self.shape, dtype=float)

        return self.__cache__   # teraz juz utworzona wiec zapisuje ja (zeby nie tworzyc tego kilka razy)

    def optimize(self):   # (***)
        if self.vectorData is not None:
            return

        import numpy as np

        self.vectorData = list()        # tworzy listy ktore jeszcze nie zostaly utworzone w konstruktorze
        self.vectorIndices = list()

        for rowId, row in enumerate(self.data):
            dV = self.diag[rowId]
            if dV != 0:
                self.vectorData.append(dV)
            for val in row:
                self.vectorData.append(val)

        for rowId, row in enumerate(self.indices):
            if self.diag[rowId] != 0:
                self.vectorIndices.append(rowId)
            for val in row:
                self.vectorIndices.append(val)

        self.vectorIndices = np.array(self.vectorIndices, dtype=int)
        self.vectorData = np.array(self.vectorData, dtype=float)
        self.rowLengths = self.row_ptr()

    def matvec(self, X):

        self.optimize()

        import numpy as np

        N = len(self.data)

        row_ptr = self.rowLengths

        F = np.zeros(N,dtype=float)

        # for row in range(N):
        #     cStart = row_ptr[row]
        #     cEnd = row_ptr[row+1]
        #     globalColumns = self.vectorIndices[cStart:cEnd]
        #     F[row] = self.vectorData[ cStart:cEnd ].dot( X[ globalColumns ] )

        F = np.array([self.vectorData[row_ptr[row]:row_ptr[row+1]].dot(X[self.vectorIndices[row_ptr[row]:row_ptr[row+1]]]) for row in range(N)])

        return F

    def relax(self, coeff=0.5):
        for i in range(self.shape[0]):
            self.data[i] = [d * coeff for d in self.data[i]]
        self.reset_cache()

    def relax2(self, Rhs, Field, coeff):
        self.diag = [d / coeff for d in self.diag]
        for i in range(len(self.diag)):
            d = self.diag[i]
            self.diag[i] /= coeff
            Rhs[i] += (self.diag[i] - d) * Field.data[i]
            # Rhs[i] /= coeff
        self.reset_cache()

    def relax3(self, Rhs, Field, coeff):
        c2 = 1. - coeff
        for i in range(self.shape[0]):
            Rhs[i] -= sum([d * c2 * Field[index] for d, index in zip(self.data[i], self.indices[i])])
            self.data[i] = [d * coeff for d in self.data[i]]
        self.reset_cache()

    def dot(self, X):
        if X.shape[1] is None:
            return self.matvec(X)
        else:
            #implement matrix - matrix multiplication
            pass

    @property
    def T(self):
        pass #implent transposition

    def relax(self, coeff=0.5):
        for i in range(self.shape[0]):
            self.data[i] = [ d*coeff for d in self.data[i] ]
        self.reset_cache()


    @property
    def shape(self):
        return (self.N, self.M)

    # definiujemy operacje na macierzach :
    def __add__(self, other):
        mat = fvMatrix(self)
        mat.add(other)
        return mat

    def __sub__(self, other):
        mat = fvMatrix(self)
        mat.sub(other)
        return mat

    def __mul__(self, other):
        mat = fvMatrix(self)
        mat.mul(other)
        return mat

    def __neg__(self):
        mat = fvMatrix(self)
        mat.mul(-1.)
        return mat


# definiowanie operacji na macierzach (obiektach tej klasy)

    def reset_cache(self):
        self.__cache__ = None
        self.vectorData = None
        self.vectorIndices = None
        self.rowLengths = None


    def add(self, other):
        for row, (ind, da) in enumerate(zip(other.indices, other.data)):
            self.diag[row] += other.diag[row]
            for col, val in zip(ind, da):             # zip dziala na kilku listach
                self.addEntry(row, col, val)

        self.reset_cache()


    def sub(self, other):
        for row, (ind, da) in enumerate(zip(other.indices, other.data)):
            self.diag[row] -= other.diag[row]
            for col, val in zip(ind, da):
                self.addEntry(row, col, -val)

        self.reset_cache()

    def mul(self, coeff):
        for i in range(self.shape[0]):
            self.diag[i] *= coeff
            self.data[i] = [coeff * t for t in self.data[i]]

        self.reset_cache()


    def mulRowWise(self, rowCoeffs):
        for i, coeff in enumerate(rowCoeffs):
            self.diag[i] *= coeff
            self.data[i] = [coeff * t for t in self.data[i]]

        self.reset_cache()


    def addEntry(self, row, col, value):

        if row == col:                        # spr czy na przekatnej jezeli tak to dopisz do diag
            self.diag[row] += value

        elif col in self.indices[row]:
            colLocalId = -1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i

            self.data[row][colLocalId] += value

        else:
            self.data[row].append(value)
            self.indices[row].append(col)

        self.reset_cache()                    # cos zmienione w macierzy


    def setEntry(self, row, col, value):

        if row == col:                        # spr czy na przekatnej jezeli tak to wpisz do diag
            self.diag[row] = value

        elif col in self.indices[row]:
            colLocalId = -1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i
            self.data[row][colLocalId] = value

        else:
            self.data[row].append(value)            # append - dodaj do data[row] wartosc value
            self.indices[row].append(col)           # dodaj do indeksow w ktorych stoi wartosc kolumny w ktorej dana wartosc value

        self.reset_cache()

    def __setitem__(self, key, value):                  # !!!!!!!!! to jako dodawanie do macierzy D[w,s]
        self.setEntry(key[0], key[1], value)

    def __getitem__(self, item):
        row, col = item
        if row == col:
            return self.diag[row]

        if col in self.indices[row]:
            colLocalId = -1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i

            return self.data[row][colLocalId]
        else:
            return 0.


    def row_ptr(self):
        if self.rowLengths is not None:
            return self.rowLengths

        rowLengths = [0]
        sumLen = 0
        for rowId, d in enumerate(self.data):
            sumLen += len(d)
            if self.diag[rowId] != 0:
                sumLen += 1

            rowLengths.append(sumLen)
        return rowLengths

    @staticmethod                            # nie trzeba tworzyc obiektu aby urzyc tej metody
    def diagonal(mesh, diagValue=1.):

        if isinstance(mesh, int):            # gdy macierz
            N = mesh
        else:
            N = mesh.n

        mat = fvMatrix(N)

        if not hasattr(diagValue, "__iter__"):
            mat.diag = [diagValue for i in range(mat.shape[0])]
        elif len(diagValue) == mat.shape[0]:
            mat.diag = diagValue

        else:
            raise Exception("Can't assign vector of length "+str(len(diagValue)) +
                            " to diagonal of matrix "+str(mat.shape))

        return mat

    # mozna skeszowac i dac wersje bez diagonali ( przeniesc dla szybszego dzialania )
    def offdiagmul(self, U):               # liczy iloczyn danego wektora U i macierzy "masowej" bez el na diagonali
        import numpy as np
        H = [0. for i in range(self.shape[0])]

        for row, (ind, da) in enumerate(zip(self.indices, self.data)):
            for col, val in zip(ind, da):
                H[row] += U[col] * val          # val juz jest wart
        return np.array(H)


# mat[5,3] = 5.         To samo
# mat.addEntry(5,3,5.)
