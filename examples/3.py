#6-dimensional algebra; 3-dimensional solvable radical, 3-dimensional Levi subalgebra

x = [Symbol("x_"+str(i+1), commutative=False) for i in range(6)]
consts = [
    ([x[0],x[1]], x[2]), ([x[0],x[2]], 0), ([x[0],x[3]], 0), ([x[0],x[4]], x[0]), ([x[0],x[5]], x[1]), 
    ([x[1],x[2]], 0), ([x[1],x[3]], -x[0]), ([x[1],x[4]], -x[1]), ([x[1],x[5]], 0),
    ([x[2],x[3]], 0), ([x[2],x[4]], 0), ([x[2],x[5]], 0),
    ([x[3],x[4]], 2*x[3]), ([x[3],x[5]], x[2]+x[4]),
    ([x[4],x[5]], 2*x[5])  
]
