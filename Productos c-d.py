def cprod(x, y):
    """
    programa que devuelve el producto de x e y en los complejos 
    """
    return vector([x[0]*y[0] - y[1]*x[1], 
                   y[1]*x[0] + x[1]*y[0]
                  ])

def cconj(x):
    """
    programa que devuelve el conjugado de x en los complejos
    """
    return vector([x[0], -x[1]])

def qprod(x,y):
    """
    programa que devuelve el producto de x e y en los cuaterniones
    """
    a, b = x[:2], x[2:]
    c, d = y[:2], y[2:]
    
    z = zero_vector(SR,4)
    z[:2] = cprod(a, c) - cprod(cconj(d), b)
    z[2:] = cprod(d, a) + cprod(b, cconj(c))
    return z

def qconj(x):
    """
    programa que devuelve el conjugado de x en los cuaterniones
    """
    p = [-i for i in x[1:].list()]
    return vector([x[0]] + p)

def qgen(st): #programa para generar un octonion cualquiera
    z_list = [var(st + str(i)) for i in range(4)]
    z = vector(z_list)
    return z

def oprod(x,y):
    """
    programa que devuelve el producto de x e y en los octoniones
    """
    a, b = x[:4], x[4:]
    c, d = y[:4], y[4:]
    
    z = zero_vector(SR,8)
    z[:4] = qprod(a, c) - qprod(qconj(d), b)
    z[4:] = qprod(d, a) + qprod(b, qconj(c))
    return z

def oconj(x):
    """
    programa que devuelve el conjugado de x en los octoniones
    """
    p = [-i for i in x[1:].list()]
    return vector([x[0]] + p)

def ogen(st):
    z_list = [var(st + str(i)) for i in range(8)]
    z = vector(z_list)
    return z

e0,e1,e2,e3,e4,e5,e6,e7 = identity_matrix(8).columns()
o_basis = [e0,e1,e2,e3,e4,e5,e6,e7]

def sprod(x,y):
    """
    programa que devuelve el producto de x e y en los sedeniones
    """
    a, b = x[:8], x[8:]
    c, d = y[:8], y[8:]
    
    z = zero_vector(SR,16)
    z[:8] = oprod(a, c) - oprod(oconj(d), b)
    z[8:] = oprod(d, a) + oprod(b, oconj(c))
    return z

def sconj(x):
    """
    programa que devuelve el conjugado de x en los sedeniones
    """
    p = [-i for i in x[1:].list()]
    return vector([x[0]] + p)

def sgen(st):
    z_list = [var(st + str(i)) for i in range(16)]
    z = vector(z_list)
    return z

s_list = []
for j in range(16):
    s_list.append(var("s"+str(j)))
s_basis = identity_matrix(16).columns()

def octs_to_sed(u,v):
    """
    devuelve el sedenion (u,v) con u y v octoniones
    """
    z = zero_vector(SR,16)
    z[:8] = u
    z[8:] = v
    return z


def o_conm(x,y):
    """
    conmutador en los octoniones
    """
    return oprod(x,y)-oprod(y,x)

def o_asoc(x,y,z):
    """
    asociador en los octoniones
    """
    p = oprod(x,y)
    q = oprod(y,z)
    return oprod(p,z)-oprod(x,q)

def s_conm(x,y):
    """
    conmutador en los sedeniones
    """
    return sprod(x,y)-sprod(y,x)

def s_asoc(x,y,z):
    """
    asociador en los sedeniones
    """
    p = sprod(x,y)
    q = sprod(y,z)
    return sprod(p,z)-sprod(x,q)