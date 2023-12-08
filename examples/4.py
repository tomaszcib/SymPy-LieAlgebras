#5-dimensional algebra; 5-dimensional nilpotent algebra Heisenberg algebra H_5
p = [Symbol("p_"+str(i+1), commutative=False) for i in range(2)]
q = [Symbol("q_"+str(i+1), commutative=False) for i in range(2)]
z = Symbol("z", commutative=False)
x = p+q+[z]
consts = [
([p[0],p[1]], 0), ([p[0],q[0]], z), ([p[0],q[1]], 0), ([p[0],z], 0),
([p[1],q[0]], 0), ([p[1],q[1]], z), ([p[1],z], 0),
([q[0],q[1]], 0), ([q[0],z], 0),
([q[1],z], 0),
]
