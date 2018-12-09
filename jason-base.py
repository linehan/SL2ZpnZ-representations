#!/usr/bin/env sage

import sys
from sage.all import *

"""
WARNINGS:

(1.)
        The Python function mod(x,y) returns a member of Z/yZ,
        not just an integer value. 

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




"""
        maps K^* -> C
"""
def jay(z, char, p, L_gen):

        # DENOTE the base field as K and the quadratic extension as L. 
        # 
        # For example, if K = F_p, L = F_{p^2}.

        s = 0;

        # 0,...,p^2-1  (L -- we filter out non-units manually to get L^*)
        for n in range(0, p**2):

                # For each element of L
                l = L_gen**n;

                # For each element of (L^*) 
                if l.is_unit() == False:
                        continue;

                # This is the field norm N:L->K. 
                # It is surjective onto K^* (see Rockmore+Lafferty).
                #
                # In a finite Galois extension L=GF(q^n) of K=GF(q),
                # for all a \in L, the field norm can be written 
                #
                #       N(a^q) = N(a) = a^{(q^n-1)/(q-1)}. 
                #
                # In our case, n=2, hence
                #
                #       N(a^q) = N(a) = a^{(q^2-1)/(q-1)} 
                #                     = a^(q+1).
                #
                norm = l**(p + 1);

                if norm == z:
                        #tr  = (l**p) + l; # Another way to compute the trace 
                        tr = l.trace();

                        # The denominator 'p' is correct (not p-1). See paper.
                        ch = complex(exp(2.0*pi*I*int(tr)/int(p)));

                        s = s + ch * char[n];

        return s / p;




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


def map1(p, g):
        # Maps k -> g^k

        M = matrix(p+1);

        M[0,0] = 1;
        M[p,p] = 1;

        for k in range(1, p):
                #print(k, mod(g**k,p));
                M[k, mod(g**k, p)] = 1; 

        print("k -> g^k");
        print(M);
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
        print("even/odd order");
        print(M);
        return M;

def principal_series_cob(p,g):
        return map2(p) * map1(p,g);



class DiscreteSeriesRepresentation():

        def __init__(self, p, g, character):
                self.p = p; # order of field K = Z/pZ
                self.g = g; # generator for K

                self.char = character;

                self.basis = [ mod(g**k, p) for k in range(1,p) ];

                # L = quadratic extension of K (only used to calculate W)
                self.L = GF(self.p**2, name="g", modulus="primitive");

                # Generator for L 
                self.L_gen = self.L.gen();

        def U(self, b):

                R = matrix(CDF, self.p-1);

                # For efficiency's sake, pre-compute the character base.
                base = complex(exp(2.0*pi*I/self.p));

                # 1,...,p-1 (K^*)
                for x in range(1, self.p):
                        j = mod(x*b, self.p);
                        R[x-1, x-1] = base ** int(j); 
                        #R[x-1, x-1] = complex(exp(2.0*pi*I*int(j)/self.p));


                return R; 

        def T(self, a):

                R = matrix(CDF, self.p-1);

                a_inverse = modinv(a, self.p);

                for x in range(1, self.p):
                        # e_x --> e_{a^{-2}x}
                        new_x = mod((a_inverse**2)*x, self.p);
                        R[x-1, new_x-1] = self.char[a_inverse]; 

                return R; # TODO need COB or not?

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

                # 1,...,p-1 (K^*)
                for x in range(1, self.p):
                        # 1,...,p-1 (K^*)
                        for y in range(1, self.p):
                                z = mod(x*y, self.p);
                                J[x-1,y-1] = jay(z, self.char, self.p, self.L_gen);

                # This is the circulant matrix we wanted (it is circulant)
                #print("J");
                #print(J.round(2));
                #print("\n");

                ## This is the diagonal matrix we wanted.
                #print("M");
                #print(M.round(2));
                #print("\n");

                return J * M; # TODO COB or not?

        def W_inv(self):
                # Can just invert the representation and be good? 
                # Because gp. homomorphism?
                # 
                # TODO: Check how close to being singular
                return self.W().inverse();



class PrincipalSeriesRepresentation():

        def __init__(self, p, g, character):
                self.p = p; # order
                self.g = g; # generator

                # Character should be a vector, though perhaps we
                # can support functions too.
                self.char = character;

                # Principal series change-of-basis matrix 
                # relative to a fixed root of unity.
                self.cob = principal_series_cob(self.p, self.g);


        """
        U() -- Represent element in U < SL_2(Z/pZ)
        ----------------------------------------
        @x    : Integer 1,...,p specifying element of U
        Return: Permutation matrix size (p+1)^2 

        NOTES
        - U isomorphic to K^+ (K considered as additive grp.)
        - Members of U can be specified by a single scalar @x
        """
        def U(self, x):

                R = matrix(CDF, self.p+1);

                # 0,...,p-1
                for u in range(0, self.p):
                        # e_u ~> e_{u-a}
                        R[u, int(mod(u - x, self.p))] = 1;

                # e_infty ~> e_infty
                R[self.p, self.p] = 1;

                # Change to Rockmore-Lafferty basis
                return self.cob * R * self.cob.inverse();

        """
        T() -- Represent element in T < SL_2(Z/pZ)
        ------------------------------------------
        @x    : Integer 1,...,p-1 specifying element of T
        Return: Complex matrix size (p+1)^2

        NOTES
        -T isomorphic to K^* (unit group of K)
        """
        def T(self, x):

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

                # Change to Rockmore-Lafferty basis
                return self.cob * R * self.cob.inverse();

        """
        W() -- Represent element W in SL_2(Z/pZ)
        ----------------------------------------
        Return: Complex matrix size (p+1)^2
        """
        def W(self):

                R = matrix(CDF, self.p+1);

                # e_0 ~> e_infty
                R[0, self.p] = 1;

                # 1,...,p-1
                for u in range(1, self.p):
                        # e_u ~> e_{-u^{-1}}
                        R[u, mod(-1*modinv(u, self.p), self.p)] = self.char[u];

                # e_infty ~> e_0 
                R[self.p, 0] = self.char[mod(-1, self.p)];

                # Change to Rockmore-Lafferty basis
                return self.cob * R * self.cob.inverse();


        def W_inv(self):
                # I think can just invert the representation and be good? 
                # because group homomorphism?
                # 
                # Should mind check how close is to being singular... don't
                # want a mess.
                return self.W().inverse();


def CHECKER():
        p = 5;
        g = 2;

        ch = principal_characters(p, g);
        r  = PrincipalSeriesRepresentation(p, g, ch.column(2));

        M = r.U(3);

        print(M.round(2));


CHECKER();
exit();

def SCHUR_CHECK():
        p = 7;
        g = 3;

        #
        # Princ series 
        #

        ch = principal_characters(p, g);
        r  = PrincipalSeriesRepresentation(p, g, ch.column(2));

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


        #
        # Discrete series turn
        #

	chp = nondecomposable_characters(p);

        print(chp.round(2));

        r = DiscreteSeriesRepresentation(p, g, chp.column(32));

        rw = r.W();

        M1 = matrix(CDF, p-1);
        M2 = matrix(CDF, p-1);

        for t in range(1,p):
                for u in range(0,p):
                        for v in range(0,p):
                                M1 += r.T(t)*r.U(u)*rw*r.U(v);

        for t in range(1,p):
                for u in range(0,p):
                        M2 += r.T(t)*r.U(u);

        print("\n");
        print((M1 + M2).round(2));
        # again should be all 0's

SCHUR_CHECK();
exit();

                                
def disc_w():
        p = 7;
        g = 3;

	chp = nondecomposable_characters(p);

        r = DiscreteSeriesRepresentation(p, g, chp.column(2));

        x = r.U(1);
        y = r.U(p-1); # p-1 == -1 mod p
        w = r.W();

        winv = w.inverse();

        print(winv.round(2));
        print("\n");
        print(w.round(2));
        print("\n");
        print((w*winv).round(2));
        print("\n");
        print("\n");
        
        print((x*y).round(2));
        print("\n");
        print((x**p).round(2));
        print("\n");
        print((y**p).round(2));

disc_w();
exit();


def DISC_SERIES_BREAKDOWN():
        p = 7;
        g = 3;

	chp = nondecomposable_characters(p);

        print(chp.round(2));

        r = DiscreteSeriesRepresentation(p, g, chp.column(32));

        print(r.U(1).round(2));       
        print("\n");
        print(r.U(p-1).round(2));
        print("\n");
        print(r.W().round(2));

        #for t in range(1,p):
                #print(r.T(t).round(2));
                #print("\n");

DISC_SERIES_BREAKDOWN();
exit();



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
