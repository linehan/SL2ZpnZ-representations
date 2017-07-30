import sys
from sage.all import *

# Representations of SL(2,Z/p^nZ)

# Parameters p = prime, n = power of the prime, D = Delta, k = integer (0 \leq k \leq n-1), sigma = integer ((sigma,p) = 1)

# Parameters are specified here, the next program will list the primitive characters based on these parameters. A representation is specified by choosing these + the primitive character

p = 5
n = 2
k = 1
D = 9
sigma = 7

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















# Finding primitive Characters

# Description: Here I assume C is a cyclic group C = <c> with order K = |C|. Every character is given by c -> e^(2pi*i*m/K) for some m in Z/KZ. The program returns the values for m that generate primitive characters.


# Finding generator <c> = C

C = [ j for j in G if norm(j) == 1 ]
K = len(C)
R0 = Integers(K)

def ord(j):
    t = j
    order = 1
    while t != (Rn(1),Rn_k(0)):
        t = mult(t,j)
        order = order + 1
    return order

def gen(C):
    t = 0
    while ord(C[t]) != K:
        t = t+1
    return C[t]

c = gen(C)


# This is the subgroup C_n-1

R1 = Integers(p^(n-1))
R2 = Integers(p^(n-1-k)) 
CL = [j for j in C if R1(j[0]) == R1(1) and R2(j[1]) == R2(0)]



# Testing if character is primitive

# Takes a character chi and a group G and returns a list with [chi(g),g]
def characterlist(chi,A):
    P = []
    for j in A:
        t = c
        a = R0(chi)
        while t != j:
            t = mult(t,c)
            a = a + chi
        P.append((j,a))
    return P
 

#returns a list char_prim with all primitive characters   
char_prim = []
for j in R0:
    B = characterlist(j,CL)
    if all( b[1] == 0 for b in B) == False:
        char_prim.append(j)

print char_prim




















# Finding the orbit space for the representation. For this set chi to be any primitive character (use program above to find these)

chi = 4


# Orbits without dependence on character

# Returns a list of tuples [alpha,n] where alpha is a representative for the orbit and n is the size of the orbit.
def orbits(G):
    T = G
    orbits = [ ]
    while len(T) != 0:
        k = T[0]
        orb_k = []
        for c in C:
            P = mult(k,c)
            if P not in orb_k:
                orb_k.append(P)
        orbits.append([orb_k[0], len(orb_k)])
        T0 = [t for t in T if t not in orb_k]
        T = T0 
    return orbits

    
# Orbits with dependence on character chi

# Returns a list of tuples [alpha,n] where alpha is a representative for the orbit and n is the size of the orbit.
def orbits_chi(G):
    T = G
    orbitschar = [ ]
    while len(T) != 0:
        k = T[0]
        orb_k1 = []
        orb_k2 = []
        for c in C:
            P = mult(k,c)
            if P not in orb_k1:
                orb_k1.append(P)
        for c in characterlist(chi,C):
            P = mult(k,c[0])
            if (P,c[1]) not in orb_k2:
                orb_k2.append((P,c[1]))
        if len(orb_k1) == len(orb_k2):
            orbitschar.append([orb_k1[0], len(orb_k1)])
        T0 = [t for t in T if t not in orb_k1]
        T = T0 
    return orbitschar

orbits_chi(G)















#Representations

# Notes:  Base field Q(zeta) still really slow. I might try magma

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
 
 
# character determined by chi on C    
def echar(a):
    chart = len(C)
    el = F.zeta(chart)^(a)
    return el

# Epsilon is constant that occurs in representation of w. 

def epsilonfinder():
    if legendre(-1) == 1 and k%2 == 1:
        epsilon = 1
    elif legendre(-1) == -1 and k%2 == 1:
        epsilon = -i
    else:
        epsilon = -1^n
    return epsilon 

eps = epsilonfinder()  


#Representations are specified by the actions of matrices A,B,w



# Action of A
# Input is a in (Z/p^nZ)^times which encodes the matrix [[a,0],[0,a^{-1}]]

def ActionA(a):
    M = []
    for l in orbits_chi(G):
        V = []
        for j in orbits_chi(G):
            Chars = characterlist(chi,C)
            H = [mult(c[0],j[0]) for c in Chars]            
            if (a*l[0][0],a*l[0][1]) in H:
                g = H.index((a*l[0][0],a*l[0][1]))                
                V.append(legendre(a)^k*echar(ZZ(Chars[g][1])))
            else:
                V.append(0)
        M.append(V)    
    return matrix(M).transpose()
        


# Action of B
# Input is b in (Z/p^nZ) which encodes the matrix [[1,b],[0,1]]

def ActionB(b):
    M = []
    for l in orbits_chi(G):
        V = []
        for j in orbits_chi(G):            
            if l == j:
                modulus = b * norm(j[0])                
                V.append(esig(ZZ(modulus)))
            else:
                V.append(0)
        M.append(V)
    return matrix(M)



# Action of w
# precomputed constants

squareroot_of_p = i*sum(legendre(a)*F.zeta(p)^a for a in (1..p-1)) if p%4 == 3 else sum(legendre(a)*F.zeta(p)^a for a in (1..p-1))

constant = squareroot_of_p^(k -2*n)*legendre(D)^(n-k)*legendre(sigma)^k*eps

Chars = characterlist(chi,C)


# No input, w = [[0,-1],[1,0]]. This calculation putters out at p = 7. 

def Actionw():
    M = []
    for l in orbits_chi(G):        
        H = [mult(c[0],l[0]) for c in Chars]
        V = []
        for j in orbits_chi(G):           
            Sum1 = 0
            for h in H:
                g = H.index(h)
                Value = esig(-2*ZZ(traceconj(j[0],h)))*echar(ZZ(Chars[g][1]))
                Sum1 = Value + Sum1
            V.append(Sum1*constant)
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

a = S([[1, 1][0,1]]);
print(a);
print(repdecomp(a));
#repdecomp(S.list()[27])


