from sympy import *
from sympy.combinatorics.partitions import IntegerPartition
from sympy.utilities.iterables import permutations
from utils import pow_to_mul

class UEA:
    """
    Class for handling arithmetics in the universal enveloping algebras

    """

    def __init__(self, g):
        """
        Initialize U(g), where g is a LieAlgebra object.

        """
        self.g = g
        allSymbols = Matrix(g.basisElements)
        usedSymbols = [list((x.T * allSymbols)[0].free_symbols)[0] for x in g.basis]
        self.basisSymbols = usedSymbols
        self.bSymbolsMatrix = Matrix(usedSymbols).T
        self.bOrder = {self.bSymbolsMatrix[i]:i for i in range(self.bSymbolsMatrix.cols)}

    def __call__(self, x):
        return self.basisSymbols[x]
    
    def multiply(self, x, y):
        return (self.bSymbolsMatrix*self.g.multiply(x,y))[0]

    def permToMonomial(self, p):
        return Mul(*(Pow(self.bSymbolsMatrix[i], p[i]) for i in range(len(p))))

    def OrderedForm(self, poly):
        """
        Rewrite polynomial `poly` in the ordered form with respect to the order of
        the basis of the underlying Lie algebra. 

        """
        if isinstance(poly, Add):
            monos = poly.args
            return Add(*[self.OrderedForm(m) for m in monos])
        m = poly
        string = []
        if not isinstance(m, Pow):
            for arg in m.args:
                dePow = pow_to_mul(arg).args
                string += [arg] if len(dePow) == 0 else dePow
        else:
            string += pow_to_mul(m).args
        if not string: c = 0
        elif isinstance(string[0], Number):
            c = string[0]
            string = string[1:]
        else: c = 1
        n = len(string)
        for i in range(n-1):
            if self.bOrder[string[i]] > self.bOrder[string[i+1]]:
                first = Mul(*string[0:i], string[i+1], string[i], *string[i+2:])
                x = Matrix([(1 if j == self.bOrder[string[i]] else 0) for j in range(len(self.bOrder))])
                y = Matrix([(1 if j == self.bOrder[string[i+1]] else 0) for j in range(len(self.bOrder))])
                second = Mul(*string[0:i], self.multiply(x,y), *string[i+2:])
                return self.OrderedForm(c*first) + self.OrderedForm(c*expand(second))
        return m


    def basisElements(self, d):
        """
        Return set of basis monomials of `self` of degree `d`

        """
        n = self.bSymbolsMatrix.cols
        a = IntegerPartition([d])
        first = a.copy()
        p = set(permutations((d,)+(0,)*(n-1), n))
        while a.next_lex() != first:
            next = a.next_lex().args[1]
            if len(next) <= n:
                next += (0,)*(n - len(next))
            p.add(next)
            p = p.union(set(permutations(next, n)))
            a = a.next_lex()
        return {self.permToMonomial(perm) for perm in p if len(perm) <= n}

