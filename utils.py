from sympy import *

def summationZero(x):
    if isinstance(x, Number): return S(0)
    elif isinstance(x, Matrix):
        return zeros(*x.shape)

def linIndSubset(v):
    """
    Pick a linearly independent subset from among iterable
    of spanning vertical vectors `v` of the same length.

    """
    ind = Matrix.hstack(*v).rref()[1]
    return [v[i] for i in ind]

def solutionSpaceBasis(solution_space):
    for i in solution_space: k = i
    outBasis = []
    for i in k.free_symbols:
        diff = k.free_symbols.difference({i})
        s = [(i,1)]
        for j in diff: s.append((j,0))
        outBasis.append(Matrix(list(k.subs(s))))
    return outBasis

def basisOfComplement(v, w, returnIndices=False):
    complementedDim = len(v)
    ind = Matrix.hstack(*(v+w)).rref()[1]
    validInds = [i-complementedDim for i in ind if i >= complementedDim]
    out = [w[i] for i in validInds]
    if not returnIndices: return out
    return (out, validInds)

def isScalarMultiplicand(v,w):
    ind = Matrix.hstack(*(v+w)).rref()[1]
    return len(ind) == 1

def E(dim,i,j):
    out = zeros(dim)
    out[i,j] = 1
    return out

def splitIntoLinear(x):
    """
    Split x into linear combination of elements.

    """
    x = S(x)
    args = x.args
    if not args:
        out = [(1,x)]
    elif isinstance(x, Mul):
        return [args]
    else:
        out = []
        for s in args:
            arg = s.args
            if not arg: out.append((1,s))
            else: out.append(arg)
    return tuple(out)

def linear(f):
    """
    Decorator used to mark that single-argument f is a linear function

    """
    def inner(*args, **kwargs):
        insideClass = len(args) == 2
        if insideClass: x = args[1]
        else: x = args
        out = list(x.atoms(Add))
        if not x.atoms(Add):
            out = [x]
        if isinstance(x, Add):
            out = list(x.args)
        combX = []
        for p in out:
            if not p.args: combX.append((1,p))
            elif isinstance(p, Mul) and isinstance(p.args[0], Number):
                combX.append((p.args[0],Mul(*p.args[1:])))
            else:
                combX.append((1,p))
        out = Add(*[i[0]*f(i[1]) for i in combX]) if not insideClass else Add(*[i[0]*f(args[0], i[1]) for i in combX])
        return out
    return inner

def isMatrixZero(A):
    for a in A:
        if a != 0: return False
    return True

def areMatricesEqual(A,B):
    for i in range(A.cols*A.rows):
        if A[i] != B[i]: 
            return False
    return True

def pow_to_mul(expr):
    """
    Convert integer powers in an expression to Muls, like a**2 => a*a.
    """
    pows = list(expr.atoms(Pow))
    repl = zip(pows, (Mul(*[b]*e,evaluate=False) for b,e in (i.as_base_exp() for i in pows)))
    return expr.subs(repl)

def matrix_to_vectors(M):
    return [M.col(i) for i in range(M.cols)]

def equiv_vector(v):
    n = len(v)
    w = zeros(n,1)
    base_value = None
    nonzeroCount = 0
    for i in range(n):
        if v[i] != 0:
            if not base_value: base_value = v[i]
            nonzeroCount += 1
    if nonzeroCount == 1 and base_value:   #isomorphic to unitary vector
        for i in range(n):
            w[i] = v[i]/base_value
        return w
    return v

def normalize(v):
    for i in range(len(v)):
        if isinstance(v[i], Number) and v[i] not in (0,1):
            for j in range(len(v)):
                v[j] *= 1/v[i]
            return v
    return v