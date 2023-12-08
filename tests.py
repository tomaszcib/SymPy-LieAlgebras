from LieAlgebra import *
import time

def timing(f):
    def wrap(*args, **kwargs):
        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        print('    {:<14} {:>.3f}s'.format('found in', (time2-time1)*1.0))
        return ret
    return wrap

@timing
def calcRepresentation(L):
    return L.Representation()

# Launch the tests for files in 'examples'
for i in range(10):
    if i < 9: L = LieAlgebra.fromFile("examples/" + str(i+1) + ".py")
    else: L = LieAlgebra.gl(i-6)
    print("Test #%d: " % (i + 1), L)
    n = len(L.basis)

    print("a) Dimensions:")
    nil = L.NilRadical()
    SR, levi = L.LeviDecomposition()
    dim_nil = len(nil.basis)
    dim_sr = len(SR.basis)
    dim_levi = len(levi.basis)
    print('    {:<14} {:>4}'.format('g', n))
    print('    {:<14} {:>4}'.format('nil(g)', dim_nil))
    print('    {:<14} {:>4}'.format('sr(g)', dim_sr))
    print('    {:<14} {:>4}'.format('Levi subalg.', dim_levi))
    print("b) Series lengths:")
    DS = len(SR.DerivedSeries())
    LCS = len(nil.LowerCentralSeries())
    print('    {:<14} {:>4}'.format('DS', DS))
    print('    {:<14} {:>4}'.format('LCS', LCS))

    print("c) Representation")
    rep = calcRepresentation(L)
    dim_rep = list(rep.values())[0].cols
    print('    {:<14} {:>4}'.format('dim', dim_rep))
    
    notes = "    "
    if dim_nil == n: notes += " nilpotent"
    if dim_sr == n: notes += " solvable"
    if dim_levi == n: notes += " semisimple"
    if L.isAbelian(): notes += " abelian"
    if notes != "    ":
        print("d) Notes:")
        print(notes)
    
    print("\n")

