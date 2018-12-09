# Representations of SL(2,Z/p^nZ)

# Parameters p = prime, n = power of the prime, D = Delta, k = integer (0 \leq k \leq n-1), sigma = integer ((sigma,p) = 1)

# Parameters are specified here, the next program will list the primitive characters based on these parameters. A representation is specified by choosing these + the primitive character

p = 7
n = 2
k = 2
D = 3
sigma = 9

Del = p^k*D

Rn = Integers(p^n)
Rn_k = Integers(p^(n-k))
G = [(a,b) for a in Rn for b in Rn_k ]

def norm(u):
    U = ZZ(u[0])^2 + Del*ZZ(u[1])^2 
    return Rn(U)  
    
def traceconj(u,v):
    U = ZZ(u[0])*ZZ(v[0]) + Del*ZZ(u[1])*ZZ(v[1]) 
    return Rn(U)  

def mult(u,v):
    U1 = ZZ(u[0])*ZZ(v[0]) - Del*ZZ(u[1])*ZZ(v[1])
    U2 = ZZ(u[0])*ZZ(v[1]) + ZZ(u[1])*ZZ(v[0])
    U = (Rn(U1),Rn_k(U2))
    return U






#Representations

# Notes: This code returns representations of the form R_k(Delta, sigma). The representations are giant and haven't been reduced to the irreducible representations determined by characters chi

F = UniversalCyclotomicField()
i = F.zeta(4)
Rp = Integers(p)

# legendre symbol    
def legendre(a):
    b = Rp(a)^((p-1)/2)
    if b == -1:
        return -1
    else:
        return ZZ(b)


# character determined by sigma, input integer a and returns value of character
def esig(a):
    el = F.zeta(p^n)^(sigma*a)
    return el
 

# Epsilon is constant that occurs in representation of w. 

def epsilonfinder():
    if legendre(-1) == 1 and k%2 == 1:
        epsilon = 1
    elif legendre(-1) == -1 and k%2 == 1:
        epsilon = -i
    else:
        epsilon = (-1)^n
    return epsilon 

eps = epsilonfinder()  


#Representations are specified by the actions of matrices A,B,w



# Action of A
# Input is a in (Z/p^nZ)^times which encodes the matrix [[a,0],[0,a^{-1}]]

def ActionA(a):
    M = []
    for l in G:
        V = []
        for j in G:          
            if (a*l[0],a*l[1]) == j:               
                V.append(legendre(a)^k)
            else:
                V.append(0)
        M.append(V)    
    return matrix(M).transpose()
        


# Action of B
# Input is b in (Z/p^nZ) which encodes the matrix [[1,b],[0,1]]

def ActionB(b):
    M = []
    for l in G:
        V = []
        for j in G:            
            if l == j:
                modulus = b * norm(j)                
                V.append(esig(ZZ(modulus)))
            else:
                V.append(0)
        M.append(V)
    return matrix(M)



# Action of w
# precomputed constants

squareroot_of_p = i*sum(legendre(a)*F.zeta(p)^a for a in (1..p-1)) if p%4 == 3 else sum(legendre(a)*F.zeta(p)^a for a in (1..p-1))

constant = squareroot_of_p^(k -2*n)*legendre(D)^(n-k)*legendre(sigma)^k*eps


# No input, w = [[0,-1],[1,0]]. 

def Actionw():
    M = []
    for l in G:        
        V = []
        for j in G:           
            Value = esig(-2*ZZ(traceconj(j,l)))
            V.append(Value*constant)
        M.append(V)
    return matrix(M)

    













#Group Decomposition

S = SL(2,Rn)

#Bruhat matrices/Decomposition

def diag(a):
    A = S([[a,0],[0,a^(-1)]])
    return A

def uptri(b):
    A = S([[1,b],[0,1]])
    return A

w = S([[0,-1],[1,0]])



# This program takes a matrix over SL(2,Z/p^nZ) and breaks it down into components. It returns an ordered tuple of matrices that multiply out to the original matrix.  

def decomp(A):
    AG = A.list()    
    if AG[1][0]%p == 0:        
        Bruhat = (w, uptri(-AG[1][0]*(AG[0][0])^(-1)), w, diag(-AG[0][0]), uptri(AG[0][1]*AG[0][0]^(-1)))        
    else:
        Bruhat = (uptri(AG[0][0]*AG[1][0]^(-1)), w, diag(AG[1][0]), uptri(AG[1][1]*AG[1][0]^(-1)))      
    return Bruhat


# This program takes an arbitrary matrix over SL(2,Z/p^nZ) and returns the corresponding representation. We should probably switch this to doing more precomputations.


def repdecomp(A):
    AG = A.list()    
    if AG[1][0]%p == 0:        
        Bruhat = Actionw()*ActionB(-AG[1][0]*(AG[0][0])^(-1))*Actionw()*ActionA(-AG[0][0])*ActionB(AG[0][1]*AG[0][0]^(-1))         
    else:
        Bruhat = ActionB(AG[0][0]*AG[1][0]^(-1))*Actionw()*ActionA(AG[1][0])*ActionB(AG[1][1]*AG[1][0]^(-1))      
    return Bruhat


# Examples. Use S.list() to produce the list of elements in SL(2, Z/p^nZ). We can select the i^th element in the list by doing S.list()[i] and check its order by S.list()[i].order(). The program repdecomp(S.list()[i]) will produce the matrix associated to this element. Alterntatively you can produce your own elements in SL(2, Z/p^nZ) by doing S([[a,b],[c,d]]).

#repdecomp(S.list()[27])

