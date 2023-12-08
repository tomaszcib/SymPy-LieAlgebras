#5-dimensional Lie algebra, 5-dimensional nilradical
x = [Symbol("x_"+str(i+1), commutative=False) for i in range(5)]
consts = [
	([x[0],x[1]], 0), ([x[0],x[2]], 0), ([x[0],x[3]], 0), ([x[0],x[4]], 0),
	([x[1],x[2]], 0), ([x[1],x[3]], 0), ([x[1],x[4]], 0),
	([x[2],x[3]], x[1]), ([x[2],x[4]], x[0]),
	([x[3],x[4]], x[2]), 
]
