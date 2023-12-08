from sympy import *
from sympy.polys.polytools import rem
from UEA import UEA
from utils import *
import itertools

class LieAlgebra:
    r"""
    Create an instance of custom finite dimenional Lie algebra with given
    variables serving as basis and definiton of Lie bracket.

    Explanation
    ===========

    This constructor takes two arguments:
    - a list of basis elements, nominatively list of objects of class
        sympy.Symbol with parameter `commutative=False` (important!)
    - a list of 2-tuples of the form ([x,y], expr) where `x` and `y`
        refer to variables used as basis elements and expr is the result
        of Lie bracket [x,y]. The function will automatically define
        [y,x] := -[x,y], so you do not have to declare all of the entries

    Note that the Lie bracket definition need to satisfy Jacobi identity.
    If it fails to do so, an exception is risen and LieAlgebra object will
    not be constructed.

    Examples
    ===========
    >>> from sympy import *
    >>> from LieAlgebra import *
    >>> x=[Symbol("x_"+str(i+1), commutative=False) for i in range(5)]
    >>> consts = [
        ([x[0],x[1]], 0), ([x[0],x[2]], 0), ([x[0],x[3]], 0),
        ([x[0],x[4]], -x[0]),  ([x[1],x[2]], 0), ([x[1],x[3]], x[0]),
        ([x[1],x[4]], 0), ([x[2],x[3]], x[1]), ([x[2],x[4]], x[2]),
        ([x[3],x[4]], -x[3])]
    >>> L = LieAlgebra(x, consts)
    >>> print(L)
    Lie algebra <x_1, x_2, x_3, x_4, x_5, >

    """

    def __init__(self, basis, struct_const, plain=False, **kwargs):
        self.basisElements = basis
        dim = len(basis)
        self.basis = [zeros(dim,1) for i in range(dim)]
        if plain: return
        for i in range(dim): self.basis[i][i] = 1
        self.struct_const = [ [] for i in range(dim)]
        for i in range(dim):
            self.struct_const[i] = [zeros(dim, 1) for j in range(dim)]
        for s in struct_const:
            out = splitIntoLinear(s[1])
            x, y = basis.index(s[0][0]), basis.index(s[0][1])
            for lin_arg in out:
                if lin_arg[1] == 0: continue
                z = basis.index(lin_arg[1])
                self.struct_const[x][y][z] = lin_arg[0]
                self.struct_const[y][x][z] = -lin_arg[0]
        self.parent = self
        self.dim = dim          #dimension of this subspace
        self.parentDim = dim    #length of parent algebra vectors - ie. dimension of topmost algebra
        # Make some assertions
        if "skip_check" in kwargs and not kwargs["skip_check"]:
            for x,y,z in itertools.product(self.basis, self.basis, self.basis):
                out = simplify(self.multiply(x, self.multiply(y,z)) + self.multiply(y, self.multiply(z,x)) + self.multiply(z, self.multiply(x,y)))
                assert (not any([i != 0 for i in out])),\
                    "Invalid Lie algebra, Jacobi identity not met!"
        assert any([not x.assumptions0["commutative"] for x in self.basisElements]),\
            "Some basis elements have commutative=True property"

    def __str__(self):
        """
        Pretty-print the LieAlgebra object.

        """
        out = "Lie algebra <"
        variables = Matrix(self.basisElements)
        for b in self.basis:
            out += (str((b.T * variables)[0]) + ", ")
        return out+">"

    @staticmethod
    def fromFile(f):
        """
        Load Lie algebra from file. Used for testing only.

        """
        loc = {}
        exec(open(f).read(), globals(), loc)
        return LieAlgebra(loc["x"], loc["consts"])

    @staticmethod
    def gl(n):
        """
        Create a LieAlgebra class instance of gl(n) general linear Lie algebra with respect
        to the standard basis.
        """
        N = n ** 2
        e = [Symbol("e"+str(i//n)+str(i%n), commutative=False) for i in range(N)]
        consts = []
        for i in range(N):
            for j in range(i+1, N, 1):
                I, J, K, L = i//n, i%n, j//n, j%n
                consts.append( ([e[i], e[j]], 
                    (1 if J == K else 0) * e[I * n + L] - (1 if L == I else 0) * e[K * n + J]))
        return LieAlgebra(e, consts, skip_check=True)


    def convertToExpr(self, x):
        return (x.T * Matrix(self.basisElements))[0]
    
    def isZero(self):
        """
        Return `true` if algebra is trivial, `false` otherwise.
        """
        return not self.basis

    def ad(self, x):
        """
        Return matrix of adjoint representation of x

        """
        return sum([Matrix.hstack(*[x[ind] * self.struct_const[ind][i] for i in range(self.dim)]) for ind in range(len(x))], zeros(self.dim))

    def multiply(self, x, y):
        """
        Lie bracket multiplication of [x,y]
        """
        out = zeros(self.dim, 1)
        for i in range(self.dim):
            for j in range(self.dim):
                mult = x[i]*y[j]
                mult *= self.struct_const[i][j]
                out += mult
                    
        return out

    def subalgebra(self, X, bugfix=False):
        """
        Create a new subalgebra of `self` spanned by set X

        """
        #BUGFIX change non-one element to 1 in vectors of x
        for x in range(len(X)):
            for i in range(len(X[x])):
                if bugfix and isinstance(X[x][i], Number) and X[x][i] not in (0,1):
                    X[x] *= 1/X[x][i]
                    break
        h = LieAlgebra(X, self.struct_const, plain=True)
        h.basisElements = self.basisElements
        h.basis = X
        h.parentDim = len(X)
        h.dim = self.parent.dim
        h.struct_const = self.struct_const
        h.parent = self.parent
        return h

    def ProductSpace(self, g, h):
        """
        Create new subalgebra of `self` spanned by basis of [g,h]
        """
        prodBasis = []
        for x in g.basis:
            for y in h.basis:
                prodBasis.append(self.multiply(x,y))
        return self.subalgebra(Matrix.hstack(*prodBasis).columnspace())

    def LowerCentralSeries(self):
        """
        Return a lower central series of `self`
        """
        L = [self]
        while not L[-1].isZero():
            L.append(self.ProductSpace(self, L[-1]))
        return L
    
    def DerivedSeries(self):
        """
        Return a derived series of `self`
        """
        L = [self]
        while not L[-1].isZero():
            L.append(self.ProductSpace(L[-1], L[-1]))
        return L

    def isAbelian(self):
        """
        Return `true` if the algebra is Abelian and `false` otherwise.
        """
        for x in self.basis:
            for y in self.basis:
                z = self.multiply(x,y)
                for i in z:
                    if i != 0: return False
        return True

    def isValidSubalgebra(self, X):
        """
        Return `true` is list of vectors X is a valid subalgebra of `self`,
        return `false` otherwise.

        """
        n = len(self.basis)
        for x in X:
            for y in X:
                if n in (Matrix.hstack(*X, self.multiply(x,y)).rref())[1]:
                    return False
        return True

    def SolvableRadical(self):
        """
        Returns a solvable radical of `self` ie. the maximal solvable ideal inside `self`
        in respect to inclusion relation.

        """
        L1 = self.ProductSpace(self, self)
        n,m = len(self.basis), len(L1.basis)
        A = Matrix(m, n, lambda j,i: trace(self.ad(self.basis[i])*L1.ad(L1.basis[j])))
        return self.subalgebra(A.nullspace())

    def NilRadical(self):
        """
        Returns a nilradical of `self` ie. the maximal nilpotent ideal inside `self`
        in respect to inclusion relation.

        """
        A = [self.ad(x) for x in self.basis]
        A1 = A.copy()
        k = len(A1)
        A = Matrix(k, k, lambda i,j: trace(A1[i] * A1[j]))
        outBasis = A.nullspace()
        radicalBasis = [zeros(self.dim) for i in range(len(outBasis))]  #Basis of a radical of (ad(x))* enveloping algebra
        for i in range(len(outBasis)):
            for j in range(len(A1)):
                radicalBasis[i] += outBasis[i][j] * A1[j]
        n,m = self.dim, len(A1)
        A = Matrix(len(self.basis), m, lambda i,j: trace(self.ad(self.basis[i]) * A1[j]))
        solution = A.nullspace()
        # Bugfix
        usedSymbols = sum([list((b.T * Matrix(self.basisElements))[0].free_symbols) for b in self.basis], [])
        usedSymbolsInd = [self.parent.basisElements.index(i) for i in usedSymbols]
        enlargedVectors = [zeros(self.dim, 1) for i in range(len(solution))]
        for i in range(len(solution)):
            for k,j in zip(usedSymbolsInd, range(len(solution[i]))):
                enlargedVectors[i][k] = solution[i][j]
        return self.subalgebra([v for v in enlargedVectors if v in self.basis])

    def isoCanonicAlgebra(self, symbol='y'):
        """
        Convert self into Lie algebra with a standard basis denoted by symbols `symbol`_i.

        """
        X = self.basis
        n = len(X)
        m = len(X[0])
        alpha = [Dummy("beta_"+str(i)) for i in range(n * len(X[0]))]
        A = Matrix(n,m, lambda i,j: alpha[i*n+j])
        equations = []
        for x, z in zip(X, matrix_to_vectors(eye(n))):
            out = A * x - z
            for i in out: equations.append(i)
        solution = solve(equations, alpha, set=True)
        # Substitute A with found solution
        A = A.subs([(_i,_j) for _i,_j in zip(solution[0], list(solution[1])[0])])
        # Substitute any parameter with 0
        A = A.subs([(_i,0) for _i in alpha])
        # Create new algebra, isomorphic to `self`
        y = [Dummy(symbol+"_"+str(i+1), commutative=False) for i in range(n)]
        y_Mat = Matrix(y)
        consts = []
        for i in range(n):
            for j in range(i + 1, n, 1):
                consts.append(([y[i],y[j]], ((A * self.multiply(X[i], X[j])).T * y_Mat)[0]))
        return (LieAlgebra(y, consts), A)

    def LeviDecomposition(self):
        """
        Return a Levi decomposition of the algebra as a 2-tuple of LieAlgebra objects:
        - a solvable radical of `self`
        - a semisimple algebra (so-called Levi subalgebra)

        """
        SR = self.SolvableRadical()
        R = SR.DerivedSeries()
        N = len(R)
        Y, map_ind = basisOfComplement(SR.basis, self.basis, returnIndices=True)
        map_BaseToY = {i:0 for i in range(self.dim)}
        for i,j in zip(map_ind, range(len(map_ind))):
            map_BaseToY[i] = j
        l = Y.copy()
        for t in range(N-1):
            V = basisOfComplement(R[t+1].basis, R[t].basis)
            #dimension of tau = dim(L) * dim(V)
            tau = [Symbol("tau_"+str(i+1)) for i in range(len(l)*len(V))]
            for j in range(len(Y)):
                l[j] += sum([tau[j*len(V)+i%len(V)] * V[i] for i in range(len(V))], zeros(self.dim, 1))
            equations = []
            equations_fix = []
            dummy_var = Matrix([Dummy(x.name) for x in self.basisElements])
            # Calculate the equations of the form [Y_i, Y_j] = sum(delta_k*Y_k)
            for I in range(len(Y)):
                for J in range(I+1, len(Y), 1):
                    left = self.multiply(l[I],l[J])         # left-hand side of the equation
                    delta = self.multiply(Y[I],Y[J])        # right-hand side of the equation
                    right = sum([delta[i]*l[map_BaseToY[i]] for i in map_ind], zeros(self.dim,1))
                    right_fix = zeros(self.dim, 1)
                    for i in set(range(self.dim))-set(map_ind):    #bugfix 
                        right_fix[i] += delta[i]
                    left -= right                           # transpose all terms to the left
                    # Apply modulo R[t] to the left-hand side
                    leftAsLinComb = expand((dummy_var.T*left)[0])
                    for r in R[t+1].basis:
                        vAsLinComb = (dummy_var.T*r)[0]
                        # In order to make this work, split the equation into part
                        # containing variables from R[t] and part not containing them
                        # We apply modulo R[t] to the former part only (SymPy bug??)
                        var = vAsLinComb.free_symbols
                        leftAsLinCombSplitWithVars = sum([leftAsLinComb.collect(var).coeff(x)*x for x in var])
                        leftAsLinCombSplitWithoutVars = simplify(leftAsLinComb - leftAsLinCombSplitWithVars)
                        remainder = rem(leftAsLinCombSplitWithVars, vAsLinComb)
                        leftAsLinComb = leftAsLinCombSplitWithoutVars + remainder
                    for i in range(len(self.basisElements)):
                        left[i] = leftAsLinComb.coeff(dummy_var[i])
                    equations += [i for i in left]
                    equations_fix += [i for i in left-right_fix]
            # Solve the equations with respect to tau
            solution = solve(equations, tau, set=True)
            if len(solution) >= 2 and not solution[1]: solution = solve(equations_fix, tau, set=True)  #bug fix
            for j in range(len(Y)):
                l[j] = l[j].subs([(_i,_j) for _i,_j in zip(solution[0], list(solution[1])[0])])
            # Substitute tau with sample values of alpha
            alpha = [Dummy("alpha_"+str(i)+str(t)) for i in range(len(tau))]
            for j in range(len(Y)):
                l[j] = l[j].subs([(_i,_j) for _i,_j in zip(tau, alpha)])
        return (SR, self.subalgebra(l))
