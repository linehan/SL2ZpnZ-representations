#!/usr/bin/env sage

import sys
from sage.all import *


"""
Helper functions
"""
def xgcd(a,b):
    x, prevx = 0, 1;
    y, prevy = 1, 0;

    while b:
        q = a/b
        x, prevx = prevx - q*x, x
        y, prevy = prevy - q*y, y
        a, b = b, a % b
    
    return a, prevx, prevy


"""
Compute the modular multiplicative inverse
"""
def modinv(a, m):
    g, x, y = xgcd(a, m)
    
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m


"""
maps K^* -> C

DENOTE the base field as K and the quadratic extension as L. 
For example, if K = F_p, L = F_{p^2}.
"""
def jay(z, char, p, L_gen):
    s = 0;

    # 0,...,p^2-1  (L -- we filter out non-units manually to get L^*)
    for n in range(0, p**2):

        # For each element of L
        l = L_gen**n;

        # For each element of (L^*) 
        if l.is_unit() == False:
            continue;

        # This is the field norm N:L->K. 
        norm = l**(p + 1);
        if norm == z:
            # The denominator 'p' is correct (not p-1). See paper.
            ch = complex(exp(2.0*pi*I*int(l.trace())/int(p)));
            s = s + ch * char[n];

    return s / p;


def chi(p, j, k):
    z = complex(exp(2.0*pi*I*float(j)*(float(k)/float(p-1))));
    z = complex(round(z.real, 10), round(z.imag, 10));
    return z;


def principal_characters(p, g):
    """
    principal_characters()
    ----------------------
    Create a matrix containing characters on the subgroup B
    @p    : Order of K 
    @g    : A generator of K*
    Return: @p x @p complex matrix with 
            col j, row k contains the value of character chi_i on element g^k.

    NOTES:
        In Rockmore/Lafferty, corresponds to
            \psi_j(a^k), 
        where 'a' is a generator of (Z/pZ)*

    BEWARE: 
        In both Rockmore-Lafferty and Rockmore-Maslen, the enumeration
        bounds are not quite correct. 
                
        R-L gives 0<j<p-1, 
        R-M gives 0<=j<=p-1,
        
        The correct bounds should be (I think) 1<=j<=p-1.
    """
    # NOTE
    # Row 0, Col 0 will always be 0
    # The characters are defined for (Z/pZ)* only.
    M = matrix(CDF, p);

    for j in range(1, p):
        for k in range(1, p):
            M[int(mod(g**k, p)), j] = chi(p, j, k); 

    return M;


def nondecomposable_characters(p):
    """
    The discrete series representations are induced from a certain
    kind of character which we compute here.

    Let K = F_p,     the finite field of order p (prime).
    Let L = F_{p^2}, the (unique) quadratic extension of K.

    There is a correspondence 

        nondecomposable           discrete series 
        characters of L     <=>   representations
                                  for GL_2(p)

    We want these nondecomposable characters.

    1. Fix a generator 'g' of L 
    2. For each element of L* (written g^j), define a character 

        ch_j: L* -> C,
        
       storing its values on each element g^n of L* as:

        ch[n,j],

       where j in {1,...,p^2-1}, 
             n in {1,...,p^2-1}. 

       The character we choose is:

         exp(2*pi*i*j*n*(1.0/(p^2-1)))
    """
    # The unique quadratic extension of F_p is F_{p^2}.
    pp = p**2;

    # Construct the field F_{p^2}
    K = GF(pp, name="g", modulus="primitive");

    # Obtain a generator for F_{p^2}
    print("waiting for generator...");
    g = K.gen();
    print("got generator");

    # We will store the nondecomposable characters as 
    # column vectors of a (p^2)x(p^2) complex matrix.
    # here 'CDF' means 'Complex Double Format'
    ch = matrix(CDF, pp);

    # TODO: Is this where things are messing up? See R+L p121 col2 parag2.
    # Do we need tr(j) rather than j?

    # For efficiency's sake, we will exponentiate this constant value. 
    char_base = complex(exp(2.0*pi*I*(1.0/(float(pp)-1.0))));

    for j in range(1, pp+1):
        # We add the j at the outer loop. 
        char_j = char_base**j;
        for n in range(1, pp+1):
            # The final character value.
            ch[n-1, j-1] = char_j**n;

    # Remove all decomposables
    # For the logic behind this, see the brief write-up
    # by Reyes, which I will soon modify for us.
    c = 0;
    for i in range(1, pp):
        if mod(p*i, pp-1) == i:
            print("is decomposable");
            # Set column 'i' to be a vector of zeros.
            ch[:,i-1] = vector([0]*pp)
            c = c + 1;
            
    print("removed "+str(c)+" decomposables");
    return ch;



def principal_characters2(p, g):
    chars = []; 
    for j in range(0, p):
        Xj = [];
        for k in range(0, p):
            Xj.append((int(mod(g**k,p)), j*k));
        chars.append(Xj);

    return chars;

def nondecomposable_characters2(p):
    pp = p**2;
    K = GF(pp, name="g", modulus="primitive");
    g = K.gen();

    chars = [];

    for j in range(1, pp+1):
        Xj = [];
        for n in range(1, pp+1):
            Xj.append((n, j*n));

        if mod(p*j, pp-1) == j:
            print("%d is decomposable" % (j));
        else:
            chars.append(Xj);

    return chars;

class DiscreteSeriesRepresentation():
    def __init__(self, p, g, character):
        self.p = p; # order of field K = Z/pZ
        self.g = g; # generator for K

        self.char = character;

        # L = quadratic extension of K (only used to calculate W)
        self.L = GF(self.p**2, name="g", modulus="primitive");

        # Generator for L 
        self.L_gen = self.L.gen();

    def U(self, b):
        R = matrix(CDF, self.p-1);

        base = complex(exp(2.0*pi*I/self.p));

        # 1,...,p-1 (K^*)
        for x in range(1, self.p):
            j = mod(x*b, self.p);
            R[x-1, x-1] = base ** int(j); 
        return R; 

    def T(self, a):
        R     = matrix(CDF, self.p-1);
        a_inv = modinv(a, self.p);

        for x in range(1, self.p):
            # e_x --> e_{a^{-2}x}
            new_x = mod((a_inv**2)*x, self.p);
            R[x-1, new_x-1] = self.char[a_inv]; 

        return R; 

    def W(self):
        pp = self.p**2;

        M = matrix(CDF, self.p-1);
        J = matrix(CDF, self.p-1);

        # Set up the diagonal
        for i in range(1, self.p):
            # This is as written in Rockmore+Lafferty
            #M[i-1,i-1] = 1.0 / self.char[i];
            # Correction to typo in Rockmore+Lafferty
            M[i-1,i-1] = self.char[modinv(i, self.p)];

        for x in range(1, self.p):
            # 1,...,p-1 (K^*)
            for y in range(1, self.p):
                z = mod(x*y, self.p);
                J[x-1,y-1] = jay(z, self.char, self.p, self.L_gen);

        return J * M; 

class PrincipalSeriesRepresentation():
    def __init__(self, p, g, character):
        self.p = p; # order
        self.g = g; # generator

        # Character should be a vector, though perhaps we
        # can support functions too.
        self.char = character;

    def U(self, x):
        """
        U() -- Represent element in U < SL_2(Z/pZ)
        ----------------------------------------
        @x    : Integer 1,...,p specifying element of U
        Return: Permutation matrix size (p+1)^2 

        NOTES
        - U isomorphic to K^+ (K considered as additive grp.)
        - Members of U can be specified by a single scalar @x
        """

        R = matrix(CDF, self.p+1);

        # 0,...,p-1
        for u in range(0, self.p):
            # e_u ~> e_{u-a}
            R[u, int(mod(u - x, self.p))] = 1;

        # e_infty ~> e_infty
        R[self.p, self.p] = 1;

        # Change to Rockmore-Lafferty basis
        return R;

    def T(self, x):
        """
        T() -- Represent element in T < SL_2(Z/pZ)
        ------------------------------------------
        @x    : Integer 1,...,p-1 specifying element of T
        Return: Complex matrix size (p+1)^2

        NOTES
        -T isomorphic to K^* (unit group of K)
        """

        R = matrix(CDF, self.p+1);

        # Scalar used for e_infty 
        scalar_e_infty = self.char[x];

        # Scalar used for e_u
        scalar_e_u = self.char[modinv(x, self.p)];

        # 0,...,p-1
        for u in range(0, self.p):
            # e_u ~> e_{(x^2)*u}
            R[u, mod((x**2)*u, self.p)] = scalar_e_u;

        # e_infty ~> e_infty 
        R[self.p, self.p] = scalar_e_infty;

        return R;

    def W(self):
        """
        W() -- Represent element W in SL_2(Z/pZ)
        ----------------------------------------
        Return: Complex matrix size (p+1)^2
        """

        R = matrix(CDF, self.p+1);

        # e_0 ~> e_infty
        R[0, self.p] = 1;

        # 1,...,p-1
        for u in range(1, self.p):
            # e_u ~> e_{-u^{-1}}
            R[u, mod(-1*modinv(u, self.p), self.p)] = self.char[u];

        # e_infty ~> e_0 
        R[self.p, 0] = self.char[mod(-1, self.p)];

        return R; 


def getit():
    p = 5;
    g = 2;

    ch1 = principal_characters(p, g);
    ch2 = nondecomposable_characters(p);
    r  = PrincipalSeriesRepresentation(p, g, ch.column(2));

    return r.T(3);


Xpc = principal_characters2(3, 2);
Xnd = nondecomposable_characters2(3);
print("Nondecomp.");
for x in Xnd:
    print(x);
print("Principal");
for x in Xpc:
    print(x);
