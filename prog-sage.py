#!/usr/bin/env sage

import sys
from sage.all import *

def permutation_from_primitive_root(p, r):
        
        # Matrix really only needs to be (p-1)x(p-1), but we extend by ident
        # because this is what we will eventually be working with.
        M = matrix(p+1);

        # From 1,...,p-1, all the powers of the unit group (F_p)*
        for i in range(1, p):
                M[i-1, mod(r**i, p)-1] = 1;

        # Extend by ident
        M[p-1,p-1] = 1;
        M[p,p]     = 1;

        return M


def permutation_princ_series_basis(p):

        # Matrix really only needs to be (p-1)x(p-1), but we extend by ident
        # because this is what we will eventually be working with.
        M = matrix(p+1);

        n_evens = floor((p-1)/2);

        for i in range(1, p):
                if mod(i,2) == 0:
                        M[(i/2) - 1, i - 1] = 1;
                else:
                        # Float prevents stupid integer division to 0 for 1/2
                        M[n_evens + ceil(float(i)/float(2)) - 1, i - 1] = 1;

        # Extend by ident
        M[p-1,p-1] = 1;
        M[p,p]     = 1;

        return M


def get_princ_series_cob(p,r):
        pp = permutation_from_primitive_root(p, r);
        pr = permutation_princ_series_basis(p);

        return pr * pp;


def get_disc_series_cob(p,r):
        return permutation_from_primitive_root(p, r);


def princ_rep_U(cob, p, a):

        rep = matrix(p+1);

        # The subgroup U is isomorphic to Z/pZ considered as a group
        # additively.
        #
        # 'u' over every element in finite field Z/pZ  
        # (remember range goes 1,..,p+1-1=p)
        # TODO: 
        # should this be additive??? Or multiplicative???)
        for u in range(1, p+1):
                
                # Shifted value of the 'u'
                new_u = mod(u - a, p);

                if new_u == 0:
                        # i.e. if 'u' == 'a' then we
                        # send it to the position of e_0,
                        # which is column 'p'.
                        new_u = p;

                rep[u-1, new_u-1] = 1;

        # e_infty is always stable.
        rep[p+1 - 1, p+1 - 1] = 1;

        # Shift into the re-ordered basis relative to a fixed root of unity.
        return cob * rep * cob.inverse();

                        
def princ_rep_T(cob, p, t):

        rep = matrix(p+1);

        # The subgroup T is isomorphic to (Z/pZ)*
        for u in range(1,p+1):

                # The shifted 'u'
                new_u = mod((t**2)*u, p);

                if new_u == 0:
                        # As above
                        new_u = p;

                rep[u - 1, new_u - 1] = 1;

        # e_infty is always stable
        rep[p+1 - 1, p+1 - 1] = 1;

        return cob * rep * cob.inverse();


############################################################### DISCRETE SERIES

"""
Helper functions
"""

#def field_norm(F, x):
        ## x should be a member of F
        #return x**(F.order() + 1);

#def field_trace(F, n, x):
        #return 0;

def xgcd(a,b):
        x, prevx = 0, 1;
        y, prevy = 1, 0;

        while b:
                q = a/b
                x, prevx = prevx - q*x, x
                y, prevy = prevy - q*y, y
                a, b = b, a % b
        
        return a, prevx, prevy

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

                        s = s + ch * char[n,0];

        return s / p;

"""
A nontrivial character K^+ -> C must be fixed for the computation of the
discrete series. This is that character.
"""
def chi_prime(j, p):
        return complex(exp(2.0*pi*I*float(j)/float(p)));

def chi_prime_power(j, q):
        return complex(exp(2.0*pi*I*float(j.trace())/float(q)));


def permutation_from_primitive_root(p, r):
        
        # Matrix really only needs to be (p-1)x(p-1), but we extend by ident
        # because this is what we will eventually be working with.
        M = matrix(p+1);

        # From 1,...,p-1, all the powers of the unit group (F_p)*
        for i in range(1, p):
                M[i-1, mod(r**i, p)-1] = 1;

        # Extend by ident
        M[p-1,p-1] = 1;
        M[p,p]     = 1;

        return M

def permutation_princ_series_basis(p):

        # Matrix really only needs to be (p-1)x(p-1), but we extend by ident
        # because this is what we will eventually be working with.
        M = matrix(p+1);

        n_evens = floor((p-1)/2);

        for i in range(1, p):
                if mod(i,2) == 0:
                        M[(i/2) - 1, i - 1] = 1;
                else:
                        # Float prevents stupid integer division to 0 for 1/2
                        M[n_evens + ceil(float(i)/float(2)) - 1, i - 1] = 1;

        # Extend by ident
        M[p-1,p-1] = 1;
        M[p,p]     = 1;

        return M

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
        g = K.gen();

        # We will store the nondecomposable characters as 
        # column vectors of a (p^2)x(p^2) complex matrix.
        # here 'CDF' means 'Complex Double Format'
        ch = matrix(CDF, pp);

        for j in range(1, pp+1):
                for n in range(1, pp+1):
                        # Note that g^n will probably be a polynomial of
                        # some kind, so the indexing of the matrix must 
                        # be done with the power (n) of the generator.
                        ch[n-1, j-1] = complex(exp(2.0*pi*I*j*n*(1.0/(float(pp)-1.0))));

        # Now remove all decomposables
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

class DiscreteSeriesRepresentation(p, g, character):

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                self.character = character;

                # Discrete series cob
                self.cob = permutation_from_primitive_root(self.p, self.g);

        ####################################################### DISCRETE SERIES

        """
        Don't need to pass to field extension, no dependence on char 
        """
        def disc_rep_U(cob, p, a):

                rep = matrix(CDF, p+1);

                for i in range(1, p):
                        rep[i-1, i-1] = chi_prime(mod(i*a, p), p);

                # Extend mx
                rep[p-1, p-1]     = 1;
                rep[p+1-1, p+1-1] = 1;
                        
                return cob * rep * cob.inverse();

        """
        Don't need to pass to field extension.
        """
        def disc_rep_T(cob, p, t):

                rep = matrix(CDF, p+1);

                for u in range(1,p):

                        # The shifted 'u'
                        new_u = mod((modinv(t,p)**2)*u, p);

                        if new_u == 0:
                                new_u = p;

                        rep[u - 1, new_u - 1] = 1;

                # Extend
                rep[p-1, p-1]     = 1;
                rep[p+1-1, p+1-1] = 1;

                # Scalar goes here

                return cob * rep * cob.inverse();

        """
        Must pass to field extension
        """
        def disc_rep_W(p):

                ch = nondecomposable_characters(p);

                pp = p**2;

                M = matrix(CDF, p);
                J = matrix(CDF, p);

                # Make this field
                K = GF(pp, name="g", modulus="primitive");

                # Make a generator.
                g = K.gen();

                # Set up the diagonal
                for i in range(1, p+1):
                        M[i-1,i-1] = 1.0 / ch[i,1];

                # Every element of K*
                for x in range(0,p):

                        # Every element of K*
                        for y in range(0, p):

                                z = mod(x*y, p);

                                J[x,y] = jay(ch[:,1], g, p, z);

                return J * M;

class PrincipalSeriesRepresentation(p, g, character):

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                self.character = character;

                # Principal series cob
                self.cob = permutation_princ_series_basis(self.p)
                         * permutation_from_primitive_root(self.p, self.g);

        ###################################################### PRINCIPAL SERIES

        def princ_rep_U(a):

                R = matrix(CDF, self.p+1);

                # The subgroup U is isomorphic to Z/pZ considered as a group
                # additively.
                #
                # 'u' over every element in finite field Z/pZ  
                # (remember range goes 1,..,p+1-1=p)
                # TODO: 
                # should this be additive??? Or multiplicative???)
                for u in range(1, self.p+1):
                        
                        # Shifted value of the 'u'
                        new_u = mod(u - a, self.p);

                        if new_u == 0:
                                # i.e. if 'u' == 'a' then we
                                # send it to the position of e_0,
                                # which is column 'p'.
                                new_u = self.p;

                        R[u-1, new_u-1] = 1;

                # e_infty is always stable.
                R[self.p+1 - 1, self.p+1 - 1] = 1;

                # NO CHARACTER SCALAR, PURE PERMUTATION 

                # Shift into the re-ordered basis relative to a fixed root of unity.
                return self.cob * R * self.cob.inverse();

        def princ_rep_T(t):

                R = matrix(self.p+1);

                # The subgroup T is isomorphic to (Z/pZ)*
                for u in range(1,self.p+1):

                        # The shifted 'u'
                        new_u = mod((t**2)*u, self.p);

                        if new_u == 0:
                                # As above
                                new_u = self.p;

                        R[u - 1, new_u - 1] = 1;

                # e_infty is always stable
                R[self.p+1 - 1, self.p+1 - 1] = 1;

                # ADD CHARACTER IN HERE

                return self.cob * R * self.cob.inverse();

        def princ_rep_W(w):
                # FILL IN
                
        ####################################################### DISCRETE SERIES


        """
        Don't need to pass to field extension, no dependence on char 
        """
        def disc_rep_U(cob, p, a):

                rep = matrix(CDF, p+1);

                for i in range(1, p):
                        rep[i-1, i-1] = chi_prime(mod(i*a, p), p);

                # Extend mx
                rep[p-1, p-1]     = 1;
                rep[p+1-1, p+1-1] = 1;
                        
                return cob * rep * cob.inverse();

        """
        Don't need to pass to field extension.
        """
        def disc_rep_T(cob, p, t):

                rep = matrix(CDF, p+1);

                for u in range(1,p):

                        # The shifted 'u'
                        new_u = mod((modinv(t,p)**2)*u, p);

                        if new_u == 0:
                                new_u = p;

                        rep[u - 1, new_u - 1] = 1;

                # Extend
                rep[p-1, p-1]     = 1;
                rep[p+1-1, p+1-1] = 1;

                # Scalar goes here

                return cob * rep * cob.inverse();

        """
        Must pass to field extension
        """
        def disc_rep_W(p):

                ch = nondecomposable_characters(p);

                pp = p**2;

                M = matrix(CDF, p);
                J = matrix(CDF, p);

                # Make this field
                K = GF(pp, name="g", modulus="primitive");

                # Make a generator.
                g = K.gen();

                # Set up the diagonal
                for i in range(1, p+1):
                        M[i-1,i-1] = 1.0 / ch[i,1];

                # Every element of K*
                for x in range(0,p):

                        # Every element of K*
                        for y in range(0, p):

                                z = mod(x*y, p);

                                J[x,y] = jay(ch[:,1], g, p, z);

                return J * M;



cob = get_disc_series_cob(7,3);

R3 = disc_rep_T(cob, 7, 3);
R2 = disc_rep_T(cob, 7, 2);
R6 = disc_rep_T(cob, 7, 6);


print(R3.round(2));
print("\n");
print(R2.round(2));
print("\n");
print((R3*R2).round(2));
print("\n");
print(R6.round(2));
print("\n");

