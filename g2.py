def Mgen(st,n): 
    """
    devuelve una matriz genérica de tamaño n
    """
    mlist=[]
    for i in range(n):
        for j in range(n):
            mlist.append(st+str(i)+str(j))
            var(st+str(i)+str(j))
    A = matrix(SR,n,n,mlist)
    return A

D = Mgen('d',8)

g2_eqs = flatten([(D * oprod(u,v) - oprod(D * u, v) - oprod(u, D * v)).list() for u in o_basis for v in o_basis])

D_g2_gen = D.subs(solve(g2_eqs, D.list(),algorithm='sympy'))

D_vars = [ent for ent in D_g2_gen.list() if len(str(ent)) == 3]

def D_vars_eqs(i):
    L = []
    for j in range(14):
        if j!=i:
            L.append(D_vars[j]==0)
        else:
            L.append(D_vars[j]==1)
    return L

basis_aux = identity_matrix(14).columns() 

g2_basis1 = []

for i in range(14):
    Q = D_g2_gen.subs(D_vars_eqs(i))
    g2_basis1 = g2_basis1+[Q]

def g2_esc(D): 
    """
    devuelve el vector de coordenadas de D en la base g2_basis1
    """
    return vector([D[5,2],D[5,3],D[5,4],D[6,1],
                   D[6,2],D[6,3],D[6,4],D[6,5],
                   D[7,1],D[7,2],D[7,3],D[7,4],
                   D[7,5],D[7,6]])
    
def coord(D, i): 
    """
    devuelve la coordenada i-ésima de D en la base g2_basis1
    """
    return g2_esc(D)[i]

def adj(X):
    """
    Devuelve la matriz de ad_X, donde X es una matriz de g2
    """
    M = zero_matrix(14,14)
    for j in range(14): #quiero calcular la columna j-ésima, para ello calculo [X,e_j] y lo escribo en escalares de la base
        C = X*g2_basis1[j]-g2_basis1[j]*X 
        C_esc = g2_esc(C) 
        for i in range(14): #el elemento (i,j) de la matriz, es la i-ésima coordenada de [X,e_j]
            M[i,j] = C_esc[i]
    return M

def ad1(k):
    """
    Devuelve la matrix de ad_Ei, en donde Ei es el i-esimo
    elto de la base de g2
    """
    A = zero_matrix(14,14)
    for j in range(14):
        B = g2_basis1[k] * g2_basis1[j] - g2_basis1[j] * g2_basis1[k]
        B_coord = g2_esc(B)
        for i in range(14):
            A[i,j] = B_coord[i]
    return A

## ad_basis1 = [ad1(k) for k in range(14)]

def killing(X,Y): 
    return (adj(X) * adj(Y)).trace()

B_mat = matrix(14,14, [killing(u,v) for u in g2_basis1 for v in g2_basis1])

Biinv = -B_mat

def g2_pi(x,y):
    """
    producto interno en g2 dado por la geometria biinvariante
    """
    return g2_esc(x)*Biinv*g2_esc(y)

def g2_norm(x):
    return sqrt(g2_pi(x,x))

g2_bo = [g2_basis1[0]]
for i in range(13):
    k= i+1
    h = g2_basis1[k] - (sum((g2_pi(g2_basis1[k],g2_bo[j])/g2_pi(g2_bo[j],g2_bo[j]))*g2_bo[j] for j in range(i)))
    g2_bo = g2_bo + [h]
    
g2_basis = [v/g2_norm(v) for v in g2_bo]

C1_list = g2_esc(g2_basis[0]).list()
for i in range(13):
    C1_list = C1_list + g2_esc(g2_basis[i+1]).list()
    
C1 = matrix(SR,14,14,C1_list)
C = transpose(C1)
Inv_C = C.inverse()

def g2_esc_bon(x):
    x_esc1 = g2_esc(x)
    x_esc_bon = Inv_C*g2_esc(x)
    return x_esc_bon

def bracket(x, y):
    return x*y - y*x