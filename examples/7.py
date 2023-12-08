#4-dimensional algebra
x = [Symbol("x_"+str(i), commutative=False) for i in range(4)]
consts = [
([x[0],x[1]], x[1]), ([x[0],x[2]], x[1]+x[2]), ([x[0],x[3]], x[2]+x[3]),
([x[1],x[2]], 0), ([x[1],x[3]], 0),
([x[2],x[3]], 0)
]
