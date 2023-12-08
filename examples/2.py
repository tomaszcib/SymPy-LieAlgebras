#3-dimensional semi-simple algebra (equal to sl(2,C))
x1,x2,h = [Symbol(s, commutative=False) for s in ("x_1", "x_2", "h")]
x = [x1,x2,h]
consts = [ ([h,x1], 2*x1), ([x1,x2], h), ([h,x2], -2*x2)]
