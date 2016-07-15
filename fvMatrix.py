# macierz dla obliczen pedu bedzie sie zmieniala wiec trzeba bedzie ja modyfikowac z kroku na krok.


class fvMatrix:
    def __init__(self, meshOrMatrix):

        if isinstance(meshOrMatrix, fvMatrix):
            self.data = list(meshOrMatrix.data)
            self.indices = list(meshOrMatrix.indices)
        else:
            self.data = [list() for i in range(meshOrMatrix.n)]
            self.indices = [list() for i in range(meshOrMatrix.n)]

        self.__cache__ = None       # nie istnieje obiekt cahe
        self.vectorData = None
        self.vectorIndices = None
        self.rowLengths = None

        self.dtype = float


    @property   # mozna uzyc jako zmiennej
    def sparse(self):
        if self.__cache__ is None:          # gdy nie istnieje to:
            self.optimize()
            from scipy.sparse import csr_matrix
            self.__cache__ = csr_matrix((self.vectorData, self.vectorIndices, self.rowLengths), dtype=float)

        return self.__cache__   # teraz juz utworzona wiec zapisuje ja (zeby nie tworzyc tego kilka razy)

    def optimize(self):
        if self.vectorData is not None:
            return

        import numpy as np

        self.vectorData = list()
        self.vectorIndices = list()

        for row in self.data:
            for val in row:
                self.vectorData.append(val)

        for row in self.indices:
            for val in row:
                self.vectorIndices.append(val)

        self.vectorIndices = np.array(self.vectorIndices, dtype=int)
        self.vectorData = np.array(self.vectorData, dtype=float)
        self.rowLengths = self.row_ptr

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


    @property
    def shape(self):
        N = len(self.data)
        return (N,N)


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

# definiowanie operacji na macierzach (obiektach tej klasy)

    def reset_cache(self):
        self.__cache__ = None
        self.vectorData = None
        self.vectorIndices = None
        self.rowLengths = None


    def add(self, other):
        for row, (ind, da) in enumerate(zip(other.indices, other.data)):
            for col, val in zip(ind, da):
                self.addEntry(row, col, val)

        self.reset_cache()


    def sub(self, other):
        for row, (ind, da) in enumerate(zip(other.indices, other.data)):
            for col, val in zip(ind, da):
                self.addEntry(row, col, -val)

        self.reset_cache()

    def mul(self, coeff):
        for i in range(len(self.data)):
            self.data[i] = [coeff * t for t in self.data[i]]

        self.reset_cache()

    def addEntry(self, row, col, value):

        if col in self.indices[row]:
            colLocalId=-1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i

            self.data[row][colLocalId] += value
        else:
            self.data[row].append(value)
            self.indices[row].append(col)

        self.reset_cache()


    def setEntry(self, row, col, value):
        if col in self.indices[row]:
            colLocalId = -1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i
            self.data[row][colLocalId] = value
        else:
            self.data[row].append(value)
            self.indices[row].append(col)

        self.reset_cache()

    def __setitem__(self, key, value):
        self.setEntry(key[0], key[1], value)

    def __getitem__(self, item):
        row, col = item
        if col in self.indices[row]:
            colLocalId = -1
            for i, id in enumerate(self.indices[row]):
                if id == col:
                    colLocalId = i

            return self.data[row][colLocalId]
        else:
            return 0.


    @property
    def row_ptr(self):
        if self.rowLengths is not None:
            return self.rowLengths

        rowLengths = [0]
        sumLen = 0
        for d in self.data:
            sumLen += len(d)
            rowLengths.append(sumLen)
        return rowLengths

    @staticmethod   # nie trzeba tworzyc obiektu aby urzyc tej metody
    def diagonal(mesh, diagValue=1.):
        mat = fvMatrix(mesh)
        for c, (i, d) in enumerate(zip(mat.indices, mat.data)):
            d.append(diagValue)
            i.append(c)
        return mat




    # mat[5,3] = 5.
    # mat.addEntry(5,3,5.)
