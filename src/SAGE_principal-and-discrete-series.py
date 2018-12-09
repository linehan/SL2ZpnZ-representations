#!/usr/bin/env sage

import sys
from sage.all import *

"""
WARNINGS:

(1.)
        The Python function mod(x,y) returns a value which is
        interpreted as being a member of that modular object.

        EXAMPLE
                a = mod(x,y);

                if a == 0, then a-1 == y, NOT 0.

        SOLUTION
                a = int(mod(x,y)) will make it a normal integer.


(2.)
        The Python function range(x,y) returns a list of values
                x,x+1,...,y-1,
        i.e., it does not include the value y. 

        EXAMPLE
                range(0,5) => [0,1,2,3,4]

FACTS:
        P.S. Reps have dimension p+1 -- p-1 elements of Z/pZ*, e_0, and e_infty
        D.S. Reps have dimension p-1 -- p-1 elements of Z/pZ* 

        P.S. Convention is g^1,g^2,...,g^p-1, since g^p-1 = g^0 
        D.S. Convention is g^0,g^1,...,g^p-2, since there is no e_0.

"""




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



def jay(char, gqe, p, z):

        # DENOTE the base field as K and the quadratic extension as L. 
        # 
        # For example, if K = F_p, L = F_{p^2}.

        pp = p**2; # order of quad. extension
        s  = 0;
        c  = 0;

        # NEED this to be over the unit group (L)*?
        # Yes, but order doesn't matter because you're summing to
        # get a number -- so just check if it's a unit, you don't
        # have to enumerate the unit group.
        for n in range(0, pp):

                t = gqe**n;

                if t.is_unit() == False:
                        continue;

                # This is the field norm for L->K. 
                # It is surjective onto (K)* (see paper).
                norm = t**(p + 1);

                if norm == z:
                        
                        # This is the trace of 't' from L->K.
                        tr  = (t**p) + t;

                        # This is the trace of 'tr' from K->F_p. (when is same?) 
                        tr2 = tr.trace();

                        # NOTE: Are there unit group shenanigans going on??
                        #print(tr, tr2);

                        # See paper. The denominator 'p' is correct.
                        ch = complex(exp(2.0*pi*I*float(tr2)/float(p)));

                        s = s + ch * char[n];

        return s / p;


"""
A nontrivial character K^+ -> C must be fixed for the computation of the
discrete series. This is that character.
"""
#def chi_prime(j, p):
        #return complex(exp(2.0*pi*I*float(j)/float(p)));

#def chi_prime_power(j, q):
        #return complex(exp(2.0*pi*I*float(j.trace())/float(q)));



#def permutation_from_primitive_root(p, r):
        
        ## Matrix really only needs to be (p-1)x(p-1), but we extend by ident
        ## because this is what we will eventually be working with.
        #M = matrix(p+1);

        ## From 1,...,p-1, all the powers of the unit group (F_p)*
        #for i in range(1, p):
                #M[i-1, mod(r**i, p)-1] = 1;

        ## Extend by ident
        #M[p-1,p-1] = 1;
        #M[p,p]     = 1;

        #return M



#def permutation_princ_series_basis(p):

        ## The initial (p-1)x(p-1) elements of the basis are 
        ## joined by two additional elements e_0, e_{\infty}. 
        #M = matrix(p+1);

        #n_evens = floor((p-1)/2);

        #for i in range(1, p):
                #if mod(i,2) == 0:
                        #M[(i/2) - 1, i - 1] = 1;
                #else:
                        ## Float prevents stupid integer division to 0 for 1/2
                        #M[n_evens + ceil(float(i)/float(2)) - 1, i - 1] = 1;

        ## Extend by ident
        #M[p-1,p-1] = 1;
        #M[p,p]     = 1;

        #return M

def generator_permutation(p, g):
        M1 = matrix(p);
        #
        # Permute according to chosen generator 'g'
        #       0,   1,   ..., p-2 
        #       
        #       g^0, g^1, ..., g^{p-2}
        #
        # NOTE must shift all left by 1 (not done here)
        # in order to avoid empty first col.
        #
        for i in range(0, p-1):
                print(i, mod(g**i, p));
                M1[i, mod(g**i, p)] = 1;

        return M1;


#def principal_series_cob(p, g):
        ## There are p+1 basis elements:
        ##       1.) 0,...,p-2 corresp. to g^0,...,g^{p-2} (elements of (F_p)*)
        ##       2.) The basis element e_0
        ##       3.) The basis element e_{infty}

        #M1 = matrix(p+1);
        #M2 = matrix(p+1);

        ## Extend by ident, e_0, e_{\infty} don't move.
        #M1[p-1,p-1] = M2[p-1,p-1] = 1;
        #M1[p,p]     = M2[p,p]     = 1;

        ##
        ## Permute according to chosen generator 'g'
        ##       0,   1,   ..., p-2 
        ##       
        ##       g^0, g^1, ..., g^{p-2}
        ##
        #for i in range(0, p-1):
                #M1[i, mod(g**i, p)-1] = 1;

        ##
        ## Permute according to Rockmore/Lafferty 
        ##       0, 1, 2, ..., p-2 
        ##
        ##       2, 4, 6, ..., p-3, 0, 1, 3, 5, ..., p-2
        ##
        #n_evens = (p-1)/2;

        #for i in range(0, p-1):
                #if i == 0:
                        #M2[n_evens-1,0] = 1;

                #elif mod(i,2) == 0:
                        #print(i, i/2);
                        #M2[(i/2) - 1, i] = 1;

                #else:
                        ## Float prevents stupid integer division to 0 for 1/2
                        #M2[n_evens + floor(float(i)/float(2)), i] = 1;

        ## Extend by ident
        #M2[p-1,p-1] = 1;
        #M2[p,p]     = 1;

        #return M2 * M1;

#
# I Believe these maps to be correct.
#

def map1(p, g):
        # Maps k -> g^k

        M = matrix(p+1);

        M[0,0] = 1;
        M[p,p] = 1;

        for k in range(1, p):
                #print(k, mod(g**k,p));
                M[k, mod(g**k, p)] = 1; 

        return M;

def map2(p):
        # Maps to Rockmore-Lafferty even/odd order 

        M = matrix(p+1);

        M[p-1, 0] = 1;
        M[p,p]    = 1;

        n_evens = (p-1)/2;

        for k in range(1, p):
                if mod(k,2) == 0:
                        M[(k/2) - 1, k] = 1;
                else:
                        #print(k, n_evens + floor(float(k)/float(2)));
                        M[n_evens + floor(float(k)/float(2)), k] = 1;
        return M;

def principal_series_cob(p,g):
        return map2(p) * map1(p,g);

#e_3 = matrix(8,1);
#e_3[6,0] = 1;

#print(e_3);
#print("\n");

#cob = principal_series_cob(7,3);

#print(cob * e_3);

#print(principal_series_cob(7,3));
#exit();



#print(map1(7,3) * map2(7));
##print("\n");
##print(map2(7) * map1(7,3));
#exit();

#def principal_series_cob3(p, g):

        #M = matrix(p+1);

        #n_evens = (p-1)/2;

        #for i in range(0, p-1):
                #if i < n_evens:
                        #print(i, 2*(i+1), mod(g**(2*(i+1)), p));
                        #M1[i, mod(g**(2*(i+1)), p)] = 1; 
                #else:
                        #print(i, 2*(i-n_evens)+1, mod(g**(i-(n_evens - 1)), p));
                        #n = int(mod(2*(i-n_evens)+1, p));

                        #M1[i, mod(g**n, p)] = 1; 

        #return M1;

#def principal_series_cob2(p, g):

        #M1 = matrix(p+1);

        #M1[0,p-1] = 1;
        #M1[p,p]   = 1;

        #n_evens = (p-1)/2;
        #print(n_evens);

        #for i in range(0, p-1):
                #if i < n_evens:
                        #print(i, 2*(i+1), mod(g**(2*(i+1)), p));
                        #M1[i, mod(g**(2*(i+1)), p)] = 1; 
                #else:
                        #print(i, 2*(i-n_evens)+1, mod(g**(i-(n_evens - 1)), p));
                        #n = int(mod(2*(i-n_evens)+1, p));

                        #M1[i, mod(g**n, p)] = 1; 

        #return M1;

#def principal_series_cob2(p, g):

        #M1 = matrix(p+1);

        #M1[0,p-1] = 1;
        #M1[p,p]   = 1;

        #n_evens = (p-1)/2;
        #print(n_evens);

        #for i in range(0, p-1):
                #if i < n_evens:
                        #print(i, 2*(i+1), mod(g**(2*(i+1)), p));
                        #M1[i, mod(g**(2*(i+1)), p)] = 1; 
                #else:
                        #print(i, 2*(i-n_evens)+1, mod(g**(i-(n_evens - 1)), p));
                        #n = int(mod(2*(i-n_evens)+1, p));

                        #M1[i, mod(g**n, p)] = 1; 

        #return M1;


#print(principal_series_cob2(7,3));
#exit();






def discrete_series_cob(p, g):
        M = matrix(p-1);

        # g^0, g^1, ..., g^{p-2}        (g^{p-1} = g^0 (mod p))
        for k in range(0, p-1):
                M[k, mod(g**k, p)] = 1;

        return M

"""
principal_characters()
----------------------
Create a matrix containing characters on the subgroup B
@p    : Order of K 
@g    : A generator of K*
Return: @px@p complex matrix described below

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
def chi(p, j, k):
        z = complex(exp(2.0*pi*I*float(j)*(float(k)/float(p-1))));
        z = complex(round(z.real, 10), round(z.imag, 10));

        return z;


def principal_characters(p, g):

        # NOTE
        # Row 0, Col 0 will always be 0
        # The characters are defined for (Z/pZ)* only.
        M = matrix(CDF, p);

        for j in range(1, p):
                for k in range(1, p):
                        M[int(mod(g**k, p)), j] = chi(p, j, k); 

        return M;


"""
There is a correspondence 

    nondecomposable characters          discrete series 
    of the unique quadratic      <=>    representations
    extension of F_p                    for GL_2(p)

Hence we want these nondecomposable characters.

So we first fix a generator 'g' of F_{p^2}, and then
for each element of (F_{p^2})*, written as g^j, define 
a character

    ch_j: (F_{p^2})* -> C,
    
storing its values on each element g^n of (F_{p^2})* as:

    ch[n,j],

where j in {1,...,p^2-1}, 
      n in {1,...,p^2-1}. 

The character we choose is the usual one:

    exp(2*pi*i*j*n*(1.0/(p^2-1)))
"""
def nondecomposable_characters(p):

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


class DSRep():

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                self.char = character;

                self.basis = [ mod(g**k, p) for k in range(1,p) ];

        ####################################################### DISCRETE SERIES

        def U(self, u):
                R = matrix(CDF, self.p-1);

                for x in range(1, self.p):
			j = mod(x*u, self.p);
                        R[x-1, x-1] = complex(exp(2.0*pi*I*int(j)/self.p));

                return R; 

        def T(self, t):
                R = matrix(CDF, self.p-1);

                scalar = self.char[modinv(t, self.p)];

                for x in range(1, self.p):
                        # e_x --> e_{t^{-2}x}
                        new_x = mod((modinv(t,self.p)**2)*x, self.p);

                        R[x-1, new_x-1] = scalar;

                return R; 

        def W(self):
                #ch = nondecomposable_characters(p);

                pp = self.p**2;

                M = matrix(CDF, self.p-1);
                J = matrix(CDF, self.p-1);

                # TODO: The instantiation of this field and the determination
                # of a generator could be done in __init__().

                # Make this field
                K = GF(pp, name="g", modulus="primitive");

                # Make a generator.
                g = K.gen();

                # Set up the diagonal
                for i in range(1, self.p):
                        M[i-1,i-1] = 1.0 / self.char[i];

                # Every element of K^*
                for x in range(1, self.p):

                        # Every element of K^*
                        for y in range(1, self.p):

                                z = mod(x*y, self.p);

                                J[x-1,y-1] = jay(self.char, g, self.p, z);

                return J * M;

        def W_inv(self):
		return self.W().inverse();




class DiscreteSeriesRepresentation():

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                self.char = character;

                # Discrete series cob
                #self.cob = permutation_from_primitive_root(self.p, self.g);
                #self.cob = discrete_series_cob(p, g);

                self.basis = [ mod(g**k, p) for k in range(1,p) ];

                #self.power_of = list();

                #for k, x in enumerate(self.basis):
                        #self.power_of[x] = k;

        ####################################################### DISCRETE SERIES

        """
        Don't need to pass to field extension, no dependence on char 
        """
        #def U(self, u):

                #R = matrix(CDF, self.p+1);

                #for i in range(1, self.p+1):
			#j = mod(i*u, self.p);

                        ## The fixed irreducible character \chi in Rockmore paper 
                        ## WE DONT EVEN PERMUTE!
                        #R[i-1, i-1] = complex(exp(2.0*pi*I*int(j)/self.p));

                ## Extend mx by ident
                #R[self.p+1-1, self.p+1-1] = 1;
                        
                #return self.cob * R * self.cob.inverse();


        """
        Don't need to pass to field extension.
        """
        #def T(self, t):

                #R = matrix(CDF, self.p+1);

                #scalar_e_u;

                #for x in range(1, self.p):
                        ## e_x --> e_{a^{-2}x}
                        #new_x = mod((modinv(t,self.p)**2)*x, self.p);

                        #R[x - 1, new_x - 1] = 1;

                ## Scalar goes here

                #return self.cob * R * self.cob.inverse();

        #####################

        def U(self, u):

                R = matrix(CDF, self.p-1);

                for x in range(1, self.p):
			j = mod(x*u, self.p);
                        R[x-1, x-1] = complex(exp(2.0*pi*I*int(j)/self.p));

                return R; # Don't believe we need COB here (see paper)

        def T(self, t):

                R = matrix(CDF, self.p-1);

                scalar = self.char[modinv(t, self.p)];

                for x in range(1, self.p):
                        # e_x --> e_{t^{-2}x}

                        new_x = mod((modinv(t,self.p)**2)*x, self.p);

                        R[x-1, new_x-1] = scalar;

                return R; # Then do we need COB here?

        def W(self):
                #ch = nondecomposable_characters(p);

                pp = self.p**2;

                M = matrix(CDF, self.p-1);
                J = matrix(CDF, self.p-1);

                # TODO: The instantiation of this field and the determination
                # of a generator could be done in __init__().

                # Make this field
                K = GF(pp, name="g", modulus="primitive");

                # Make a generator.
                g = K.gen();

                # Set up the diagonal
                for i in range(1, self.p):
                        M[i-1,i-1] = 1.0 / self.char[i];

                # Every element of K^*
                for x in range(1, self.p):

                        # Every element of K^*
                        for y in range(1, self.p):

                                z = mod(x*y, self.p);

                                J[x-1,y-1] = jay(self.char, g, self.p, z);

                return J * M;
                # TODO: NO NEED FOR COB HERE?

        """
        Must pass to field extension
        """
        #def W(self):
                ##ch = nondecomposable_characters(p);

                #pp = self.p**2;

                #M = matrix(CDF, self.p);
                #J = matrix(CDF, self.p);

                ## Make this field
                #K = GF(pp, name="g", modulus="primitive");

                ## Make a generator.
                #g = K.gen();

                ## Set up the diagonal
                #for i in range(1, self.p+1):
                        #M[i-1,i-1] = 1.0 / self.char[i,1];


                ## Every element of K*
                #for x in range(0,self.p):

                        ## Every element of K*
                        #for y in range(0, self.p):

                                #z = mod(x*y, self.p);

                                #J[x,y] = jay(self.char[:,1], g, self.p, z);

                #return J * M;
                # TODO: NO NEED FOR COB HERE?


        def W_inv(self):
		return self.W().inverse();


class PrincipalSeriesRepresentation():

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                # Character should be a vector, though perhaps we
                # can support functions too.
                self.char = character;

                # Principal series cob
                self.cob = principal_series_cob(self.p, self.g);

                #self.p2v = [ mod(self.g**i, self.p) for i in range(1, self.p)]
                #self.v2p = [ 0 for i in range(0, self.p) ];

                #for i in range(0, self.p):
                        #self.v2p[mod(self.g**i, self.p)] = i;

                #print(self.p2v);
                #print(self.v2p);

        """
	U()
	---
	Generate the representation of an element in the subgroup U of SL_2(p)

	@a    : An integer from 1,...,p specifying the element of U
	Return: A (p+1)x(p+1) matrix of complex values

	NOTES
		-It is sufficient to specify an integer 1,...,p 
		 in order to specify the corresponding element of U. 

        	-The subgroup U is isomorphic to K^+ (K considered as 
		 an additive group) 
        """
        def U(self, a):

                R = matrix(CDF, self.p+1);

                for u in range(0, self.p):
                        # Send e_u --> e_{u-a}
                        R[u, int(mod(u - a, self.p))] = 1;

                # Send e_infty --> e_infty
                R[self.p, self.p] = 1;

                # Shift into the re-ordered basis relative to a fixed root of unity.
                return self.cob * R * self.cob.inverse();

        """
	T()
	---
	Generate the representation of an element in the subgroup T of SL_2(p)

	@t    : An integer from 1,...,p-1 specifying the element of T.
	Return: A (p+1)x(p+1) matrix of complex values

	NOTES
		-The subgroup T is isomorphic to K^* (the unit group of K)
        """
        def T(self, t):

                R = matrix(CDF, self.p+1);

                # Scalar used for e_infty 
                scalar_e_infty = self.char[t];

                # Scalar used for e_u
                scalar_e_u = self.char[modinv(t, self.p)];

                for u in range(0, self.p):
                        # Send e_u --> e_{(t^2)*u}
                        R[u, mod((t**2)*u, self.p)] = scalar_e_u;

                # Send e_infty --> e_infty 
                R[self.p, self.p] = scalar_e_infty;

                return self.cob * R * self.cob.inverse();

	"""
	W()
	---
	Generate the representation of the 'w' element in SL_2(p) 
	
	Return: A (p+1)x(p+1) matrix of complex values.
	"""
        def W(self):

                R = matrix(CDF, self.p+1);

                # Each element 1,...,p-1
                for u in range(1, self.p):
                        # Send e_u --> e_{-u^{-1}}
                        R[u, mod(-1*modinv(u, self.p), self.p)] = self.char[u];

                # Send e_0 --> e_infty
                R[0, self.p] = 1;

                # Send e_infty => e_0 * scalar 
                R[self.p, 0] = self.char[mod(-1, self.p)];

                return self.cob * R * self.cob.inverse();

	"""
	W_inv()
	-------
	Generate the representation of the 'w^{-1}' element in SL_2(p) 
	
	Return: A (p+1)x(p+1) matrix of complex values.
	"""
        def W_inv(self):
                return self.W().inverse();
                

def EXPPP():
        ch = nondecomposable_characters(7);
        r = DiscreteSeriesRepresentation(7, 3, ch.column(1));

        A = matrix(CDF, 7-1);

        for u in range(0,7):
                a = r.U2(u);
                print(a.round(2));
                print("\n");
                A += a;

        print(A.round(2));
        print("\nT\n");

        for t in range(1,7):
                a = r.T2(t);
                print(a.round(2));
                print("\n");

        print("W\n");
        w = r.W2();
        print(w.round(2));
        print("\n");


#EXPPP();
#exit();

def second_largest_eigenvalue(matrix):

        evals = matrix.eigenvalues();

        first  = complex(0);
        second = complex(0);

        for i in range(0, len(evals)):
                x = complex(evals[i]);
                if abs(x) > abs(second):
                        if abs(x) > abs(first):
                                second = first;
                                first = x;
                        else:
                                second = x;

        return second;


def primes(n):
    """ Returns  a list of primes < n """
    sieve = [True] * n
    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)/(2*i)+1)
    return [2] + [i for i in xrange(3,n,2) if sieve[i]]


PRIMES = [
	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 
	71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
];

SMALL_PRIMES = [
	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 
	71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
];

SMALL_PRIMES_AND_GENERATORS = [
	(3, 2), (5, 2), (7, 3), (11, 2), (13, 2), (17, 3)#, (19, 2), (23, 5), 
	#(29, 2), (31, 3), (37, 2), (41, 6), (43, 3), (47, 5), (53, 2), (59, 2), 
	#(61, 2), (67, 2), (71, 7), (73, 5), (79, 3), (83, 2), (89, 3), (97, 5), 
	#(101, 2), (103, 5), (107, 2), (109, 6), (113, 3), (127, 3)
];

PRIMES_AND_GENERATORS = [
	(3, 2), (5, 2), (7, 3), (11, 2), (13, 2), (17, 3), (19, 2), (23, 5), 
	(29, 2), (31, 3), (37, 2), (41, 6), (43, 3), (47, 5), (53, 2), (59, 2), 
	(61, 2), (67, 2), (71, 7), (73, 5), (79, 3), (83, 2), (89, 3), (97, 5), 
	(101, 2), (103, 5), (107, 2), (109, 6), (113, 3)#, (127, 3)
];



def rep(R, p, a, b, c, d):
    """ 
    Implements Bruhat decomposition 
    @a    :
    @b    :
    @c    :
    @d    :
    Return:
    """
    ring = IntegerModRing(p);

    a = ring(a);
    b = ring(b);
    c = ring(c);
    d = ring(d);

    print(int(d*c**(-1)));

    if c % p == 0:
        bruhat = R.W()               \
               * R.U(int(-c*(a)**(-1)))   \
               * R.W()               \
               * R.T(int(-a))             \
               * R.U(int(b*a**(-1)));
    else:
        bruhat = R.U(int(a*c**(-1)))      \
               * R.W()               \
               * R.T(int(c))              \
               * R.U(int(d*c**(-1)));
    
    return bruhat; 

def EXPERIMENT_g3_disc():
	LD = [];
	LP = [];

	g3_csv = open("g3_spectrum_disc.csv", "w");
	g3_csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");

	g3_line = 0;

        for p, g in PRIMES_AND_GENERATORS:
        #for p, g in [(3,2), (5,2), (7,3)]:

                print("in prime:"+str(p)+" generator:"+str(g));

		ch = nondecomposable_characters(p);

                # PRINCIPAL SERIES ALL eigenvalues
                for i in range(1, p):
		        R = DiscreteSeriesRepresentation(p, g, ch.column(i));

                        elem_1 = rep(R, p, 1, 1, -1, 0);
                        elem_2 = rep(R, p, 0, -1, 1, 1);
                        elem_3 = R.W();
                        elem_4 = R.W().inverse();

                        g3 = elem_1 + elem_2 + elem_3 + elem_4;
                        
                        for evalue, e, multiplicity in g3.eigenvectors_right(): 
				g3_csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					g3_line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				g3_line += 1;

EXPERIMENT_g3_disc();
exit();


def EXPERIMENT_g3():
	LD = [];
	LP = [];

	g3_csv = open("g3_spectrum.csv", "w");
	g3_csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");

	g3_line = 0;

        for p, g in PRIMES_AND_GENERATORS:
        #for p, g in [(3,2), (5,2), (7,3)]:

                print("in prime:"+str(p)+" generator:"+str(g));

		chp = principal_characters(p, g);

                # PRINCIPAL SERIES ALL eigenvalues
                for i in range(1, p):
		        R = PrincipalSeriesRepresentation(p, g, chp.column(i));

                        elem_1 = rep(R, p, 1, 1, -1, 0);
                        elem_2 = rep(R, p, 0, -1, 1, 1);
                        elem_3 = R.W();
                        elem_4 = R.W().inverse();

                        g3 = elem_1 + elem_2 + elem_3 + elem_4;
                        
                        for evalue, e, multiplicity in g3.eigenvectors_right(): 
				g3_csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					g3_line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				g3_line += 1;

EXPERIMENT_g3();
exit();


#chd = nondecomposable_characters(13);
#exit();

def EXPERIMENT_random_spectrum():

	csv = open("random_spectrum3.csv", "w");

	csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");

	line = 0;

        for p, g in PRIMES_AND_GENERATORS:

                if p <= 11:
                        continue;

                print("in prime:"+str(p)+" generator:"+str(g));

		chp = principal_characters(p, g);

                # PRINCIPAL SERIES ALL eigenvalues
                for i in range(1, p):
		        r = PrincipalSeriesRepresentation(p, g, chp.column(i));

                        m = r.T(11)*r.U(7);

                        for evalue, e, multiplicity in m.eigenvectors_right(): 
				csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				line += 1;

EXPERIMENT_random_spectrum();
exit();


def EXPERIMENT_sle_growth_princ():

        csv = open("sle_growth_princ.csv", "w");

        csv.write("line, p, evalue_scaled, evalue_real_part, evalue_imag_part\n");

        line = 0;

        for p, g in PRIMES_AND_GENERATORS:

                print("in prime:"+str(p)+" generator:"+str(g));

		ch = principal_characters(p, g);

                best = 0;



                # PRINCIPAL SERIES
                for i in range(1, p):
		        r = PrincipalSeriesRepresentation(p, g, ch.column(i));

                        m = r.U(1)+r.U(p-1)+r.W()+r.W().inverse();
                        
                        ev = second_largest_eigenvalue(m);

                        if abs(ev) > abs(best):
                                best = ev;

                csv.write("%d, %d, %05f, %05f, %05f\n" % (
                        line, 
                        p, 
                        abs(best.real/4.0),
                        best.real,
                        best.imag
                ));

                line += 1;


EXPERIMENT_sle_growth_princ();
exit();

def EXPERIMENT_sle_growth_disc():

        for p, g in SMALL_PRIMES_AND_GENERATORS:

                print("in prime:"+str(p)+" generator:"+str(g));

		ch = nondecomposable_characters(p);

                best = 0;

	        csv = open("sle_growth_disc.csv", "w");

	        csv.write("line, p, evalue_scaled, evalue_real_part, evalue_imag_part\n");

                line = 0;

                # DISCRETE SERIES
                for i in range(1, p):
		        r = DiscreteSeriesRepresentation(p, g, ch.column(i));

                        m = r.U(1)+r.U(p-1)+r.W()+r.W().inverse();
                        
                        ev = second_largest_eigenvalue(m);

                        if abs(ev) > abs(best):
                                best = ev;

                csv.write("%d, %d, %05f, %05f, %05fi\n" % (
                        line, 
                        p, 
                        float(abs(best.real))/4.0,
                        best.real,
                        best.imag
                ));

                line += 1;

EXPERIMENT_sle_growth_disc();
exit();


def EXPERIMENT_spectrum_2():
	LD = [];
	LP = [];

	g1_csv = open("g1_spectrum_disc.csv", "w");

	g1_csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");

	g1_line = 0;

        for p, g in SMALL_PRIMES_AND_GENERATORS:

                print("in prime:"+str(p)+" generator:"+str(g));

		chp = nondecomposable_characters(p);

                # DISCRETE SERIES ALL eigenvalues
                for i in range(1, p):
		        r = DiscreteSeriesRepresentation(p, g, chp.column(i));

			W    = r.W();
			Winv = W.inverse();

                        g1 = r.U(1)       + r.U(p-1)     + W + Winv;
                        #g2 = r.U((p+1)/2) + r.U((p-1)/2) + W + Winv;
                        
                        for evalue, e, multiplicity in g1.eigenvectors_right(): 
				g1_csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					g1_line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				g1_line += 1;

EXPERIMENT_spectrum_2();
exit();



def EXPERIMENT_spectrum_1():
	LD = [];
	LP = [];

	g1_csv = open("g1_spectrum.csv", "w");
	g2_csv = open("g2_spectrum.csv", "w");

	g1_csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");
	g2_csv.write("line, p, multiplicity, evalue_scaled, evalue_real_part, evalue_imag_part\n");

	g1_line = 0;
	g2_line = 0;

        for p, g in PRIMES_AND_GENERATORS:

                print("in prime:"+str(p)+" generator:"+str(g));

		chp = principal_characters(p, g);

                # PRINCIPAL SERIES ALL eigenvalues
                for i in range(1, p):
		        r = PrincipalSeriesRepresentation(p, g, chp.column(i));

			W    = r.W();
			Winv = W.inverse();

                        g1 = r.U(1)       + r.U(p-1)     + W + Winv;
                        g2 = r.U((p+1)/2) + r.U((p-1)/2) + W + Winv;
                        
                        for evalue, e, multiplicity in g1.eigenvectors_right(): 
				g1_csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					g1_line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				g1_line += 1;

                        for evalue, e, multiplicity in g2.eigenvectors_right(): 
				g2_csv.write("%d, %d, %d, %05f, %05f, %05fi\n" % (
					g2_line, 
					p, 
					multiplicity,
					float(evalue.real())/4.0,
					evalue.real(),
					evalue.imag()
                                ));
				g2_line += 1;
#EXPERIMENT_spectrum_1();
#exit();
		

def EXPERIMENT_QUA3():

	LD = [];
	LP = [];

        for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:

		k = GF(p, modulus="primitive");
		g = int(k.gen());
                print("in prime:"+str(p)+" generator:"+str(g));

		chp = principal_characters(p, g);
		#chd = nondecomposable_characters(p);

                # PRINCIPAL SERIES second largest eigenvalues
                #for i in range(1, p):
			#r = PrincipalSeriesRepresentation(p, g, chp.column(i));

                        #m = r.U(1)+r.U(p-1)+r.W() + r.W().inverse();
                        
                        #ev = second_largest_eigenvalue(m);

                        #LP.append((float(abs(ev))/4.0, p));

                # PRINCIPAL SERIES ALL eigenvalues
                for i in range(1, p):
		        r = PrincipalSeriesRepresentation(p, g, chp.column(i));

                        #m = r.U(1)+r.U(p-1)+r.W() + r.W().inverse();
                        m = r.U((p+1)/2)+r.U((p-1)/2)+r.W() + r.W().inverse();
                        
                        for ev in m.eigenvalues(): 
                                LP.append((float(ev.real())/4.0, p));

		
	GP = list_plot(LP);
	GP.save("spectrum_cayley_g2_princ_all.pdf");

#EXPERIMENT_QUA3();
#exit();

def EXPERIMENT_QUA2():

	L  = [];
	LD = [];
	L2 = [];

        for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:

		k = GF(p, modulus="primitive");
		g = int(k.gen());
                print("in prime:"+str(p)+" generator:"+str(g));

		ch  = principal_characters(p, g);
		#chd = nondecomposable_characters(p);

                best = 0;
                bestd= 0;

                # PRINCIPAL SERIES
                for i in range(1, p):
		        r = PrincipalSeriesRepresentation(p, g, ch.column(i));

                        m = r.U(1)+r.U(p-1)+r.W() + r.W().inverse();
                        
                        ev = second_largest_eigenvalue(m);

                        if abs(ev) > abs(best):
                                best = ev;

                # DISCRETE SERIES
                #for i in range(1, p):
			#r = DiscreteSeriesRepresentation(p, g, chd.column(i));

                        #m = r.U2(1)+r.U2(p-1)+r.W2() + r.W2().inverse();
                        
                        #ev = second_largest_eigenvalue(m);

                        #if abs(ev) > abs(bestd):
                                #bestd = ev;

                #LD.append(abs(bestd));
                L2.append(float(abs(best))/4.0);

		
	G = list_plot(L2);
	G.save("basis_g1_secondlargest_track_princ.pdf");

	#G2 = list_plot(LD);
	#G2.save("basis_g1_secondlargest_disc.pdf");

EXPERIMENT_QUA2();
exit();

def EXPERIMENT_QUA():

	L = [];
	L2 = [];

        #for p in primes(100):
        for p in range(7,8):
		k = GF(p, modulus="primitive");
		g = int(k.gen());

		ch = principal_characters(p, g);

                # TUwU part
                for t in range(1, p):
                        for u in range(0, p):
                                # Character values determine reps
                                for i in range(1, p):
                                        print(str(t)+":"+str(u)+":"+str(i));
		                        r = PrincipalSeriesRepresentation(p, g, ch.column(i));
                                        s = r.T(t)*r.U(u)*r.W()*r.U(u);

                                        ev = second_largest_eigenvalue(s);

		                        L.append(ev);
		                        L2.append(abs(ev));

                # TU part
                for t in range(1, p):
                        for u in range(0, p):
                                # Character values determine reps
                                for i in range(1, p):
                                        print(str(t)+":"+str(u)+":"+str(i));
		                        r = PrincipalSeriesRepresentation(p, g, ch.column(i));
                                        s = r.T(t)*r.U(u);

                                        ev = second_largest_eigenvalue(s);

		                        L.append(ev);
		                        L2.append(abs(ev));

		
	G = list_plot(L);
	G.save("thinger.pdf", xmin=-1.1, ymin=-1.1, xmax=1.1, ymax=1.1);

	G2 = list_plot(L2);
	G2.save("thinger2.pdf");

EXPERIMENT_QUA();
exit();


def EXPERIMENT_1():

	L = [];
	L2 = [];

        for p in primes(100):
		k = GF(p, modulus="primitive");
		g = int(k.gen());

		ch = principal_characters(p, g);

		#for u in range(1, p):
		# ONE IS AS GOOD AS ANY OTHER I SUPPOSE - EVALUES DONT CHANGE
		# BASED ON WHAT CHARACTER IS CHOSEN.
		r  = PrincipalSeriesRepresentation(p, g, ch.column(1));

		# Generating set
		a = r.U(1);
		b = r.U(mod(-1, p));
		w    = r.W();
		winv = r.W_inv();

		L.append(second_largest_eigenvalue(a));
		L.append(second_largest_eigenvalue(b));
		L.append(second_largest_eigenvalue(w));
		L.append(second_largest_eigenvalue(winv));

		L2.append(abs(second_largest_eigenvalue(a)));

		
	G = list_plot(L);
	G.save("thinger.pdf", xmin=-1.1, ymin=-1.1, xmax=1.1, ymax=1.1);

	G2 = list_plot(L2);
	G2.save("thinger2.pdf");

EXPERIMENT_1();
exit();


def TEST_PRINC_U():
        p = 7;
        g = 3;

	ch = principal_characters(p, g);

        r = PrincipalSeriesRepresentation(p, g, ch.column(1));

        A = matrix(p+1);

        for k in range(0, p):
                a = r.U(k);
                A += a;

        print(A);

def TEST_PRINC_T():
        p = 7;
        g = 3;

	ch = principal_characters(p, g);

        r = PrincipalSeriesRepresentation(p, g, ch.column(1));

        for a in range(1, p):
                t = r.T(a);
                print(t.round(2));
                print("\n");

def TEST_PRINC_W():
        p = 7;
        g = 3;

	ch = principal_characters(p, g);

        r = PrincipalSeriesRepresentation(p, g, ch.column(1));

        print(r.W().round(2));
        print("\n");
        print(r.W_inv().round(2));
        print("\n");


def SCHUR_CHECK():
        p = 7;
        g = 3;

	ch = principal_characters(p, g);
        r = PrincipalSeriesRepresentation(p, g, ch.column(2));

        rw = r.W();

        M1 = matrix(p+1, p+1);
        M2 = matrix(p+1, p+1);

        for t in range(1,p):
                for u in range(0,p):
                        for v in range(0,p):
                                M1 += r.T(t)*r.U(u)*rw*r.U(v);

        for t in range(1,p):
                for u in range(0,p):
                        M2 += r.T(t)*r.U(u);

        print((M1 + M2).round(2));
        # Should be all 0's

                                
SCHUR_CHECK();           
exit();


def TEST_HOMOMORPHISM():
        p = 7;
        g = 3;

	ch = principal_characters(p, g);

        r = PrincipalSeriesRepresentation(p, g, ch.column(1));

        ############################# U HOMOMORPHISM

        U1 = matrix(2,2);
        U1[0,0] = 1;
        U1[0,1] = 2;
        U1[1,0] = 1;
        U1[1,1] = 0;

        U2 = matrix(2,2);
        U1[0,0] = 1;
        U1[0,1] = 4;
        U1[1,0] = 1;
        U1[1,1] = 0;

        ua = r.U(2);
        ub = r.U(4);

        q = ua*ub;

        z = r.U(6);
        print(ua*ub);
        print("\n");
        print(z);
        print("\n");
        print(ua * ub == z);

        ############################# T HOMOMORPHISM

        T1 = matrix(2,2);
        T1[0,0] = 2;
        T1[0,1] = 0;
        T1[1,0] = -2;
        T1[1,1] = 0;

        T2 = matrix(2,2);
        T2[0,0] = 4;
        T2[0,1] = 0;
        T2[1,0] = -4;
        T2[1,1] = 0;

        ta = r.T(2);
        tb = r.T(4);

        q = ta*tb;

        z = r.T(1);
        print((ta*tb).round(5));
        print("\n");
        print(z);
        print("\n");
        print((ta * tb).round(5) == z);

        ############################# W HOMOMORPHISM (??)

         



#TEST_HOMOMORPHISM();
#exit();
#TEST_PRINC_U();
TEST_PRINC_T();
#TEST_PRINC_W();
exit();


#def EXPERIMENT_2(p):

            #k = GF(p, modulus="primitive");
            #g = int(k.gen());

            #print(p, g);

	    #ch = principal_characters(p, g);

            #for j in range(1, 2):
                #r = PrincipalSeriesRepresentation(p, g, ch.column(j));

                ##print(ch.round(2));

                #A = matrix(p+1);

                #for k in range(0, p):
                    #a = r.U(k);
                    #A += a;

                    ##print("\n");
                    ##print(a.round(2));
                    ##print("\n");

                #print(A);



	    # Generating set
            #a = r.U(1);
            #a = r.T(1);
            ##b = r.U(mod(-1, p));
            ##w    = r.W();
            ##winv = r.W_inv();

            #print("\n");
            #print(a.round(2));
            #print("\n");
            #print(b.round(2));
            #print("\n");
            #print(w.round(2));
            #print("\n");
            #print(winv.round(2));

EXPERIMENT_2(7);
exit();



def experiment(p, g):

        for j in range(1, p):
                print("\nseries "+str(j));
                for k in range(1, p):

                        val = mod(g**k, p);

                        # These floats are necessary in Python 2.7 in order
                        # to make the division be non-truncating. This was
                        # fixed in Python 3.x
                        res = complex(exp(2.0*pi*I*float(j)*(float(k)/float(p-1))));

                        res = complex(round(res.real, 10), round(res.imag, 10));

                        print(val, res, res*res);

experiment(7,3);
exit();





ch = nondecomposable_characters(7);
r = DiscreteSeriesRepresentation(7, 3, ch.column(1));

m = r.T(2);
print(m.round(2));

exit();


EXPERIMENT_1();
exit();



ch = principal_characters(7,3);

r = PrincipalSeriesRepresentation(7, 3, ch.column(1));
m = r.U(3);

print(m.round(2));
print("\n");

exit();



print(ch.ncols());



exit();

#exit();

for i in range(1,7):
        r = PrincipalSeriesRepresentation(7, 3, ch);
        m = r.T(i);
        print(m.eigenvalues());
        print("\n");

exit();


r = PrincipalSeriesRepresentation(7, 3, chi_1);

R3 = r.T(3);

print(R3);

                
#cob = get_disc_series_cob(7,3);
#R3 = disc_rep_T(cob, 7, 3);
#R2 = disc_rep_T(cob, 7, 2);
#R6 = disc_rep_T(cob, 7, 6);
#print(R3.round(2));
#print("\n");
#print(R2.round(2));
#print("\n");
#print((R3*R2).round(2));
#print("\n");
#print(R6.round(2));
#print("\n");
