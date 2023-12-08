#6-dimensional solvable algebra; 4-dimensional nilradical
x = [Symbol("x_"+str(i), commutative=False) for i in range(6)]
consts = [
([x[0],x[1]], 0), ([x[0],x[2]], 0), ([x[0],x[3]], 0), ([x[0],x[4]], -3*x[0]), ([x[0],x[5]], 0),
([x[1],x[2]], 0), ([x[1],x[3]], 0), ([x[1],x[4]], -Rational(1,2)*x[1]), ([x[1],x[5]], 0),
([x[2],x[3]], 0), ([x[2],x[4]], 0), ([x[2],x[5]], 0),
([x[3],x[4]], -x[3]), ([x[3],x[5]], 0),
([x[4],x[5]], 0),
]
