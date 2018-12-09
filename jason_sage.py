##!/usr/bin/env sage
#!/usr/bin/sage -python

import sys
from sage.all import *

# Representations of SL(2,Z/p^nZ)


# Parameters are specified here. A representation is specified by choice of 
# these parameters + a primitive character.
#
# Parameters:
#       p     = prime, 
#       n     = power of the prime, 
#       D     = Delta, 
#       k     = integer (0 \leq k \leq n-1), 
#       sigma = integer ((sigma,p) = 1)
#
#p = 5
#n = 2
#k = 1
#D = 9
#sigma = 7

p = 5
n = 2 
k = 1 
D = 9
sigma = 7 

Del = (p**k)*D

print("Running SL2REP with p=%d n=%d k=%d D=%d sigma=%d Del=%d" % (p, n, k, D, sigma, Del));

print("Constructing the group G...");
Z    = IntegerRing();
Rn   = IntegerModRing(p**n)
Rn_k = IntegerModRing(p**(n-k))

G = [(a,b) for a in Rn for b in Rn_k ]

def norm(u):
    x = Z(u[0])**2 + Del*Z(u[1])**2 
    return Rn(x)  
    
def traceconj(u,v):
    x = Z(u[0])*Z(v[0]) + Del*Z(u[1])*Z(v[1]) 
    return Rn(x)  

def mult(u,v):
    w1 = Z(u[0])*Z(v[0]) - Del*Z(u[1])*Z(v[1])
    w2 = Z(u[0])*Z(v[1]) + Z(u[1])*Z(v[0])
    w  = (Rn(w1),Rn_k(w2))
    return w 

        #for key in ["p","n","k","sigma","delta_prime"]: 
            #if key not in kwargs:
                #print("You must supply TanakaGroup with %s" % (key));
                #exit();
            #else:
                #self.__setattr__(key, kwargs[key]);


#def isprime(n):
    #''' check if integer n is a prime '''

    ## make sure n is a positive integer
    #n = abs(int(n))

    ## 0 and 1 are not primes
    #if n < 2:
        #return False

    ## 2 is the only even prime number
    #if n == 2: 
        #return True    

    ## all other even numbers are not primes
    #if not n & 1: 
        #return False

    ## range starts with 3 and only needs to go up 
    ## the square root of n for all odd numbers
    #for x in range(3, int(n**0.5) + 1, 2):
        #if n % x == 0:
            #return False

    #return True


    


#def check_tanaka_params(p, n, k, sigma, delta_prime): 

    #assert isinstance(p, int),   "Parameter 'p' must be an integer";
    #assert isinstance(n, int),   "Parameter 'n' must be an integer";
    #assert isinstance(k, int),   "Parameter 'k' must be an integer";

    #assert p > 0,                "Parameter 'p' must be positive";
    #assert n > 1,                "Parameter 'n' must be > 1";
    #assert k >= 0,               "Parameter 'k' must be >= 0";
    #assert k <= n,               "Parameter 'k' must be <= n";

    #assert isprime(p),           "Parameter 'p' must be prime!"; 

    #assert sigma % p != 0,       "Parameter 'sigma' must be != 0 (mod p)";
    #assert delta_prime % p != 0, "Parameter 'delta_prime' must be != 0 (mod p)";

    #return {"p":p, "n":n, "k":k, "sigma":sigma, "delta_prime":delta_prime};

    

#class TanakaMonoid(): 
    #""" A monoid under multiplication, a group under addition. """

    #def __init__(self, p, n, k, sigma, delta_prime):
            
        #self.param = check_tanaka_params(p, n, k, sigma, delta_prime);

        #self.delta = (self.param.p**self.param.k)*self.param.delta_prime;

        #self.Z    = IntegerRing();
        #self.Rn   = IntegerModRing(self.param.p**self.param.n);
        #self.Rn_k = IntegerModRing(self.param.p**(self.param.n-self.param.k));
        ## TODO: Should these not be Zn and Zn_k ?

        #self.members = [ (a,b) for a in self.Rn for b in self.Rn_k ];

    #def __iter__(self):
        #""" Allow user to iterate over the object like a list of tuples """
        #return iter(self.members);

    #def size(self, g=False): 
        #""" 'Order' is a little too suggestive perhaps. """
        #return len(self.members);

    #def norm(self, u):
        #""" Norm is just traceconj(u,u) = <u,u> in Tanaka's notation """
        #x = self.Z(u[0])**2 + self.delta*self.Z(u[1])**2;
        #return self.Rn(x);
    
    #def traceconj(self, u,v):
        #x = self.Z(u[0])*self.Z(v[0]) + self.delta*self.Z(u[1])*self.Z(v[1]);
        #return self.Rn(x);

    #def mult(self, u, v):
        #w1 = self.Z(u[0])*self.Z(v[0]) - self.delta*self.Z(u[1])*self.Z(v[1]);
        #w2 = self.Z(u[0])*self.Z(v[1]) + self.Z(u[1])*self.Z(v[0]);
        #w  = (self.Rn(w1),self.Rn_k(w2));
        #return w;

    #def get_subgroup(self, subset):
        #return TanakaSubgroup(
            #subset      = subset,
            #p           = self.param.p, 
            #n           = self.param.n, 
            #k           = self.param.k, 
            #sigma       = self.param.sigma, 
            #delta_prime = self.param.delta_prime,
        #);

    #def subgroup_C(self):
        #return self.get_subgroup([g for g in self.members if self.norm(g) == 1]);


#class TanakaSubgroup(TanakaMonoid):

    #def __init__(self, subgroup):
        #self.members   = subgroup;
        #self.generator = None;

    #def order(self, g=None):
        #if g == None:
            #order = len(self.members);
        #else:
            #order = 1;
            #t = g;
            ## BUG: If C is not cyclic (we have not proven it always is), then 
            ## there is some value such that this loop will never halt. 
            #while t != (self.Rn(1), self.Rn_k(0)):
                #t = self.mult(t, g);
                #order = order + 1;

        #return order

    #def get_generator(self):
        #if self.generator != None:
            ## Cache it
            #return self.generator;
        #else:
            #size = self.order();
            #t = 0 
            #while self.order(self.members[t]) != size:
                #t = t+1
            #return self.members[t];

#K = len(C)

#print(C);

## R0 = Z/KZ, ring of integers modulo K
#R0 = Integers(K)


## TODO (BUG): If C is NOT cyclic, then there is some value such that ord() 
## will loop forever. 
#def ord(j):
    #t = j
    #order = 1
    #while t != (Rn(1),Rn_k(0)):
        #t = mult(t,j)
        #order = order + 1
    #return order

#def gen(C):
    #t = 0 
    #while ord(C[t]) != K:
        #t = t+1
    #return C[t]

#print("Finding a generator for C...");

#c = gen(C)

#print("C = <(%d,%d)>" % (c[0],c[1]));


## This is the subgroup C_n-1

#print("Generating subgroup C_{n-1}...");
#R1 = Integers(p**(n-1))
#R2 = Integers(p**(n-1-k)) 
#CL = [j for j in C if R1(j[0]) == R1(1) and R2(j[1]) == R2(0)]



## Testing if character is primitive

## Takes a character chi and a group G and returns a list with [chi(g),g]
#def characterlist(chi,A):
    #P = []
    #for j in A:
        #t = c
        #a = R0(chi)
        #while t != j:
            #t = mult(t,c)
            #a = a + chi
        #P.append((j,a))

    #return P
 

#print("Generating primitive character table...");
##returns a list char_prim with all primitive characters   
#char_prim = []
#for j in R0:
    #B = characterlist(j,CL)
    #if all( b[1] == 0 for b in B) == False:
        #char_prim.append(j)

##print char_prim




#GG = TanakaGroup(params={
        #p = 5,
        #n = 2,
        #k = 1,
        #delta_prime = 9,
        #sigma = 7
#});


#for u in GG:
    #print(u);


#exit();





# Finding primitive Characters

# Description: 
#       Here I assume C is a cyclic group C = <c> of order K = |C|. 
#       Every character is given by c -> e^(2pi*i*m/K) for some m in Z/KZ. 
#       The program returns the values for m that generate primitive 
#       characters.


# Finding generator <c> = C
print("Constructing subgroup C...");
C = [ j for j in G if norm(j) == 1 ]
K = len(C)

print(C);



# R0 = Z/KZ, ring of integers modulo K
R0 = Integers(K)


# TODO (BUG): If C is NOT cyclic, then there is some value such that ord() 
# will loop forever. 
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

print("Finding a generator for C...");

c = gen(C)

print("C = <(%d,%d)>" % (c[0],c[1]));


# This is the subgroup C_n-1

print("Generating subgroup C_{n-1}...");
R1 = Integers(p**(n-1))
R2 = Integers(p**(n-1-k)) 
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
 

print("Generating primitive character table...");
#returns a list char_prim with all primitive characters   
char_prim = []
for j in R0:
    B = characterlist(j,CL)
    if all( b[1] == 0 for b in B) == False:
        char_prim.append(j)

#print char_prim

#print characterlist(4, C); 
#exit();




















# Finding the orbit space for the representation. 
# For this set chi to be any primitive character 
# (use program above to find these)

chi = 4


# Orbits without dependence on character

# Returns a list of tuples [alpha,n] where alpha is a representative 
# for the orbit and n is the size of the orbit.


# JASON QUESTION:
# Why is this outputting a tuple? where do we use the 
# size of the orbit? i.e. the second value in the tuple?
#
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

# Returns a list of tuples [alpha,n] where alpha is a representative 
# for the orbit and n is the size of the orbit.
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

#orbits_chi(G)















#Representations

# Notes:  Base field Q(zeta) still really slow. I might try magma
print("Generating fields for representations...");

F = UniversalCyclotomicField()
i = F.zeta(4)   # <--- JASON: is this 4 the chosen chi value?
    # F.zeta(n) is an alias for F.gen(n), returning standard nth root of unity
    # In fact, F.gen(n,k) is supported, giving kth power of standard nth root
    # of unity. This may be faster than what we do below in echar()?
Rp = Integers(p)

# legendre symbol    
def legendre(a):
    b = Rp(a)**((p-1)/2)
    if b == -1:
        return -1
    else:
        return Z(b)


# character determined by sigma, input integer a and returns value of character
def esig(a):
    el = F.zeta(p**n)**(sigma*a)
    return el
 
 
# character determined by chi on C     --- JASON: How is 'chi' determining anything?
def echar(a):
    chart = len(C)
    el = F.zeta(chart)**(a)
    return el

# Epsilon is constant that occurs in representation of w. 
def epsilonfinder():
    if legendre(-1) == 1 and k%2 == 1:
        epsilon = 1
    elif legendre(-1) == -1 and k%2 == 1:
        epsilon = -i
    else:
        epsilon = -1**n
    return epsilon 

print("Finding epsilon to use for 'w'...");
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
                # Can we not multiply the entire matrix by legendre(a)**k
                # once we are all done? it's just a scalar.
                V.append(legendre(a)**k*echar(Z(Chars[g][1])))
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
                V.append(esig(Z(modulus)))
            else:
                V.append(0)
        M.append(V)
    return matrix(M)



# Action of w
# precomputed constants
print("Computing tables for action of 'w'...");

squareroot_of_p = i* sum(legendre(a)*F.zeta(p)**a for a in range(1,p)) if p%4 == 3 else sum(legendre(a)*F.zeta(p)**a for a in range(1, p))

constant = squareroot_of_p**(k -2*n)*legendre(D)**(n-k)*legendre(sigma)**k*eps

Chars = characterlist(chi,C)


# No input, w = [[0,-1],[1,0]]. This calculation putters out at p = 7. 

def ActionW():
    M = []
    for l in orbits_chi(G):        
        H = [mult(c[0],l[0]) for c in Chars]
        V = []
        for j in orbits_chi(G):           
            Sum1 = 0
            for h in H:
                g = H.index(h)
                Value = esig(-2*Z(traceconj(j[0],h)))*echar(Z(Chars[g][1]))
                Sum1 = Value + Sum1
            V.append(Sum1*constant)
        M.append(V)
    return matrix(M)


    def W(self):
        """ 
        Action of W 
        Return:

        NOTE
        The result is cached.
        """
        M = [];
        for o1 in self.orbit_reps:        
            H = [self.G.mult_fast(c, o1) for (c, cx) in self.X];
            V = [];
            for o2 in self.orbit_reps:
                Sum1 = 0
                for h in H:
                    g = H.index(h)
                    # pull things out of this sum
                    Value = self._e_sigma(-2*(self.G.traceconj_fast(o2,h))) \
                          * self._e_chi((self.X[g][1]));
                    Sum1 = Value + Sum1;
                V.append(Sum1*self.W_constant)
            M.append(V)
        self.W_cached = matrix(M);




















#Group Decomposition

print("getting group");
S = SL(2,Rn)

#Bruhat matrices/Decomposition

def diag(a):
    A = S([[a,0],[0,a**(-1)]])
    return A

def uptri(b):
    A = S([[1,b],[0,1]])
    return A

w = S([[0,-1],[1,0]])



# This program takes a matrix over SL(2,Z/p^nZ) and breaks it down into components. It returns an ordered tuple of matrices that multiply out to the original matrix.  

def decomp(A):
    AG = A.list()    
    if AG[1][0]%p == 0:        
        Bruhat = (w, uptri(-AG[1][0]*(AG[0][0])**(-1)), w, diag(-AG[0][0]), uptri(AG[0][1]*AG[0][0]**(-1)))        
    else:
        Bruhat = (uptri(AG[0][0]*AG[1][0]**(-1)), w, diag(AG[1][0]), uptri(AG[1][1]*AG[1][0]**(-1)))      
    return Bruhat


# This program takes an arbitrary matrix over SL(2,Z/p^nZ) and returns the corresponding representation. We should probably switch this to doing more precomputations.


def repdecomp(A):
    AG = A.list()    
    if AG[1][0]%p == 0:        
        Bruhat = ActionW()*ActionB(-AG[1][0]*(AG[0][0])**(-1))*ActionW()*ActionA(-AG[0][0])*ActionB(AG[0][1]*AG[0][0]**(-1))         
    else:
        Bruhat = ActionB(AG[0][0]*AG[1][0]**(-1))*ActionW()*ActionA(AG[1][0])*ActionB(AG[1][1]*AG[1][0]**(-1))      
    return Bruhat


# Examples. Use S.list() to produce the list of elements in SL(2, Z/p^nZ). We can select the i^th element in the list by doing S.list()[i] and check its order by S.list()[i].order(). The program repdecomp(S.list()[i]) will produce the matrix associated to this element. Alterntatively you can produce your own elements in SL(2, Z/p^nZ) by doing S([[a,b],[c,d]]).

#a = S([[1, 1],[0,1]]);
#print("repdecomping");
#print(repdecomp(a));
#repdecomp(S.list()[27])

def test_output():
    return (ActionB(3), ActionA(4), ActionW(), constant);

def test_primitives():
    return char_prim;



###### TESTS #####

#def TEST_C(C):

    #fail = 0;

    ## C IS A GROUP?
    #for i in range(0, len(C)):
        #found_inverse = 0;
        #for j in range(0, len(C)):
            #w = mult(C[i],C[j]);

            #found_inverse = (found_inverse or w == (1,0));

            #if w not in C:
                #print("[FAIL] Test group-ness of subgroup C"); 
                #print("       mult((%d,%d),(%d,%d)) not in C" % (C[i][0],C[i][1],C[j][0],C[j][1]));
                #fail = 1;
            #elif j == (len(C)-1) and found_inverse == 0:
                #print("[FAIL] Test group-ness of subgroup C"); 
                #print("       (%d,%d) has no inverse?" % (C[i][0],C[i][1]));
                #fail = 1;

    #print("[PASS] Test group-ness of subgroup C"); 

    ## C IS ABELIAN?
    #for i in range(0, len(C)):
        #for j in range(0, len(C)):
            #if mult(C[i],C[j]) != mult(C[j],C[i]):
                #print("[FAIL] Test abelian-ness of subgroup C"); 
                #fail = 1;

    #print("[PASS] Test abelian-ness of subgroup C"); 

    ## C HAS THE PREDICTED ORDER?
    #if k == 0:
        #if len(C) == p**(n-1)*(p - (-1*Del/p)):
            #print("[PASS] Test order of subgroup C"); 
        #else:
            #print("[FAIL] Test order of subgroup C"); 
            #fail = 1;
    #elif 1 <= k and k <= n-1:
        #if len(C) == 2*p**(n-k):
            #print("[PASS] Test order of subgroup C"); 
        #else:
            #print("[FAIL] Test order of subgroup C"); 
            #fail = 1;
    #else:
        #print("[FAIL] Test order of subgroup C"); 
        #print("       MALFORMED parameter k=%d" % (k));
        #fail = 1;

    #return not fail;


#TEST_C(C);
