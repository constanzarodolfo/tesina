def octs_to_Z(u,v,x,y):
    """
    devuelve un elemento de R32 dados los octoniones u, v, x e y
    """
    z = zero_vector(SR,32)
    z[:16] = octs_to_sed(u,v)
    z[16:] = octs_to_sed(x,y)
    return z

def Z_to_octs(p):
    """
    devuelve una lista de 4 octoniones dado un elemento de R32
    """
    return[p[:8],p[8:16],p[16:24],p[24:]]

def Zprod(u,v,x,y):
    """
    devuelve el producto en los sedeniones de (u,v) y (x,y)
    """
    z1 = octs_to_sed(u,v)
    z2 = octs_to_sed(x,y)
    return sprod(z1,z2)

def Zprod2(p):
    """
    devuelve el producto en los sedeniones de los sedeniones dados por las primeras 16 entradas
    de p y las 16 últimas
    """
    x = p[:16]
    y = p[16:]
    return sprod(x,y)

def act_g2_Z(A,p):
    """
    dado un elemento p de R32, devuelve (Ap1,Ap2,Ap3,Ap4)
    """
    z = zero_vector(SR,32)
    z[:8] = A*(p[:8])
    z[8:16] = A*(p[8:16])
    z[16:24] = A*(p[16:24])
    z[24:] = A*(p[24:])
    return z