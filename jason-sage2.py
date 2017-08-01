#!/usr/bin/sage -python

import sys
from sage.all import *


# Helpful static functions that use the classes below. 
###############################################################################

#def _get_orbits(group, subset):
    ## Computes the orbits of the action of a group on a subset of
    ## that group.
    ##
    ## Returns a list of tuples [alpha,n] where alpha is a representative 
    ## for the orbit and n is the size of the orbit.
    #orbits = [];
    #G      = group.members;

    #while len(G) != 0:
        #u     = G[0]; # get first thing in here
        #orb_u = [];
        #for v in subset:
            #uv = group.mult(u,v);
            #if uv not in orb_u:
                #orb_u.append(uv);
        #orbits.append([orb_u[0], len(orb_u)]);
        #G = [g for g in G if g not in orb_u];

    #return orbits


# Returns a list of tuples [alpha,n] where alpha is a representative 
# for the orbit and n is the size of the orbit.

#
#    TO OBTAIN A BASIS FOR IndVx:
#    ----------------------------
#    1. Let C act on G by left multiplication and 
#       decompose the set G into orbits. 
#    
#    2. Let b_i be a set of representatives for these 
#       orbits.
#
#    3. Define a basis for IndVx as 
#
#           e_i = \sum_{c\in C} x(c)(b_i * c)
#
#    We must be careful here: G is a monoid so it is possible to have
#
#        b_i * c_1 = b_i * c_2,
#
#    in which case if 
#
#        x(c_1) != x(c_2), 
#
#    we simply define e_i = 0.

def _get_orbit_representatives_IndVx(G, C, x):
    """
    _get_orbit_representatives_IndVx()
    ----------------------------------
    Compute the orbits of the action of C on G by left-multiplication 

    @G: The monoid G
    @C: The subgroup C<G
    @x: The index of a primitive/principal character 

    Return: List of b_i in G, for each [b_i] in the orbit space G/C. 
    """
    orbit_reps = [];
    G_elements = G.members; # Copying this isn't efficient 

    while len(G) != 0:
        g      = G_elements[0];
        orb_g  = [];
        orb_xg = [];

        for c in C:
            gc = G.mult(g, c);
            if gc not in orb_g:
                orb_g.append(gc);

        for (c, xc) in C.get_character(x):
            gc = group.mult(g, c);
            if (gc, xc) not in orb_xg:
                orb_xg.append((gc, xc));

        if len(orb_g) == len(orb_xg):
            """
            This happens precisely when there is v1, v2 in subgroup
            such that for this u, 
                u*v1 == u*v2 
            and
                chi(v1) != chi(v2).

            In such a case, there will be fewer elements in orb_u
            than in orb_xu.

            This implies that f(u) = 0 (see our paper draft), and
            this leads to badness ... ? So we do not count it.

            (THIS IS SORT OF CRUDE -- CONSIDER REFIT)
            """
            orbit_reps.append([orb_g[0], len(orb_g)])
        else:
            print("Orbit length not equal!");

        G_elements = [u for u in G_elements if u not in orb_g];

    return orbit_reps;


def _get_chi(j, G, generator):
    """
    Let C be the complex numbers.

    A character X:G->C on a multiplicative abelian group 

        G = <g>, |G| = n

    must map the generator g to one of the complex n-th 
    roots of unity (there are n): 

        e^(2*pi*i*(j/n)) : 0 <= j <= n-1.

    The choice of j determines the values of X on any 
    other element u = g^k of G:

        X(g^k) = (X(g))^k = e^(2*pi*i*(j/n)*k),

    and every complex n-th root of unity corresponds to
    a character on G. Hence we can index the characters 
    of G by this choice of j.

    If we let B = e^(2*pi*i*(1/n)), we can see that 

        X_j(g^k) = B^(j*k),

    and hence if we know B, we can store a character X_j
    as a list of tuples:

        (u = g^k, j*k),

    and recover the complex value later in the computation. 
    """
    if j >= G.get_order():
        print("Character index 'j' is out of range.");

    R0 = Integers(G.get_order());

    chi = [];
    gen = generator;

    for u in G:
        v   = gen;
        jk  = R0(j); # Why do we do this here. Ask Ben.
        while v != u:
            v  = G.mult(v, gen);
            jk = jk + j;  

        chi.append((u, jk));

    return chi;


###############################################################################

class TanakaParams():

    def isprime(self, n):
        ''' check if integer n is a prime '''

        n = abs(int(n))     # make sure n is a positive integer

        if n < 2:
            return False    # 0 and 1 are not primes

        if n == 2: 
            return True     # 2 is the only even prime number

        if not n & 1: 
            return False    # all other even numbers are not primes

        # range starts with 3 and only needs to go up 
        # the square root of n for all odd numbers
        for x in range(3, int(n**0.5) + 1, 2):
            if n % x == 0:
                return False

        return True

    def __init__(self, p, n, k, sigma, delta_prime):
        assert isinstance(p, int),   "Parameter 'p' must be an integer";
        assert isinstance(n, int),   "Parameter 'n' must be an integer";
        assert isinstance(k, int),   "Parameter 'k' must be an integer";

        assert p > 0,                "Parameter 'p' must be positive";
        assert n > 1,                "Parameter 'n' must be > 1";
        assert k >= 0,               "Parameter 'k' must be >= 0";
        assert k <= n,               "Parameter 'k' must be <= n";

        assert self.isprime(p),      "Parameter 'p' must be prime!"; 

        assert sigma % p != 0,       "Parameter 'sigma' must be != 0 (mod p)";
        assert delta_prime % p != 0, "Parameter 'delta_prime' must be != 0 (mod p)";

        self.p           = p;
        self.n           = n;
        self.k           = k;
        self.sigma       = sigma;
        self.delta_prime = delta_prime;

###############################################################################

class TanakaMonoid(): 
    """ A monoid under multiplication, a group under addition. """

    def __init__(self, p, n, k, sigma, delta_prime):
        self.param = TanakaParams(p, n, k, sigma, delta_prime);

        self.delta = (self.param.p**self.param.k)*self.param.delta_prime;

        self.Z    = IntegerRing();
        self.Rn   = IntegerModRing(self.param.p**self.param.n);
        self.Rn_k = IntegerModRing(self.param.p**(self.param.n-self.param.k));
        # TODO: Should these not be Zn and Zn_k ?

        self.members = [ (a,b) for a in self.Rn for b in self.Rn_k ];

    def __iter__(self):
        """ Allow user to iterate over the object like a list of tuples """
        return iter(self.members);

    def size(self, g=False): 
        """ 'Order' is a little too suggestive perhaps. """
        return len(self.members);

    def norm(self, u):
        """ Norm is just traceconj(u,u) = <u,u> in Tanaka's notation """
        """ Note we promote the elements to full integers, then demote
            them again into Rn at the end. 
        """
            be enforced to be elements
        x = self.Z(u[0])**2 + self.delta*self.Z(u[1])**2;
        return self.Rn(x);
    
    def traceconj(self, u,v):
        x = self.Z(u[0])*self.Z(v[0]) + self.delta*self.Z(u[1])*self.Z(v[1]);
        return self.Rn(x);

    def mult(self, u, v):
        w1 = self.Z(u[0])*self.Z(v[0]) - self.delta*self.Z(u[1])*self.Z(v[1]);
        w2 = self.Z(u[0])*self.Z(v[1]) + self.Z(u[1])*self.Z(v[0]);
        w  = (self.Rn(w1),self.Rn_k(w2));
        return w;

    def get_subgroup_C(self):
        return TanakaCyclicAbelianSubgroup(
            monoid=self, 
            subset=[g for g in self.members if self.norm(g) == 1]
        );

    def get_orbits_of_action_on(self, subgroup):
        if self.param != subgroup.monoid.param:
            print("Not my subgroup!");

        return _get_orbits(self, subgroup);

    def get_orbits_of_action_on_with_char(self, subgroup, j):
        if self.param != subgroup.monoid.param:
            print("Not my subgroup!");

        return _get_orbit_representatives_IndVx(self, subgroup, j);

        

###############################################################################



class TanakaCyclicAbelianSubgroup(TanakaMonoid):
    """ 
    Used for the subgroups CL < C < G 
    """
    def __init__(self, monoid, subset, generator=None):
        TanakaMonoid.__init__(self, 
            p           = monoid.param.p, 
            n           = monoid.param.n, 
            k           = monoid.param.k, 
            sigma       = monoid.param.sigma, # TODO: sigma not used in classes 
            delta_prime = monoid.param.delta_prime
        );

        self.monoid               = monoid;
        self.members              = subset;

        self.generator            = generator;
        self.characters           = {};
        self.primitive_characters = [];

        # Are these the R_0, R_1 in Tanaka's paper (p.125) ?
        self.R0 = Integers(self.get_order());
        self.R1 = Integers(self.param.p**(self.param.n-1));
        self.R2 = Integers(self.param.p**(self.param.n-1-self.param.k));


    def get_order(self, element=None):
        """ WARNING: Halts only if the subgroup is cyclic. """ 
        if element == None:
            return len(self.members);
        else:
            u = element;
            o = 1;
            # BUG: Will not halt if group not cyclic (no proof yet in case of C) 
            while u != (self.Rn(1), self.Rn_k(0)):
                u = self.mult(u, element);
                o = o + 1;

        return o; 


    def get_generator(self):
        """ WARNING: Halts only if the subgroup is cyclic. """ 
        if self.generator == None:
            t = 0 
            while self.get_order(self.members[t]) != self.get_order():
                t = t+1
            self.generator = self.members[t];
            
        return self.generator;


    def get_character(self, chi):
        if chi not in self.characters:
            self.characters[chi] = _get_chi(chi, self, self.get_generator());

        return self.characters[chi];


    def get_primitive_characters(self):
        """ This operation works only for the particular CL<C... not class-y """ 
        if len(self.primitive_characters) == 0:

            R0  = self.R0;
            R1  = self.R1;
            R2  = self.R2;

            # The subgroup CL < C (also called C_{n-1})
            # Obviously C cyclic abelian => CL cyclic abelian. 
            CL = TanakaCyclicAbelianSubgroup(
                monoid    = self.monoid, 
                generator = self.get_generator(),
                subset    = [
                    j for j in self if R1(j[0]) == R1(1) and R2(j[1]) == R2(0)
                ], 
            );

            for chi in R0:
                if all(u[1] == 0 for u in CL.get_character(chi)) == False:
                    self.primitive_characters.append(chi);

        return self.primitive_characters;



###############################################################################

def ActionA(a):
    M = []
    # We are building the matrix one entry at a time by looping over
    # the indices l,j. 
    #
    # A(a)phi(u) = (a/p)^k phi(ua) where (./p) is the Legendre symbol
    # 
    # phi(ua) is a matrix so we need to apply to each of its entries.
    # apparently these are indexed by orbits_chi.
    #
    # orbits_chi(G) = [[alpha,n]] where alpha is a rep of the orbit and
    # n is the size of the orbit.
    #
    # then get the list of characters for C under chi.
    #
    # Then form H = {cj : c\in C, j represents an orbit in orbit_chi}
    #
    # Then if one of these is (a*l1, a*l2) = (l1,l2)*a,
    #
    # look up the character value under that, do some stuff to it, and
    # make that the value at (l,j).
    #
    for l in orbits_chi(G):
        V = []
        for j in orbits_chi(G):
            Chars = characterlist(chi,C)
            H = [mult(c[0],j[0]) for c in Chars]            
            #
            # If there is c in C such that
            #   c*j = a*l,
            # Then (it is unique? or it doesn't matter which) 
            #   a*l in H
            # and set
            #   A(a)_(l,j) = (a/p)^k * echar(chi(a*l))
            #
            if (a*l[0][0],a*l[0][1]) in H:
                g = H.index((a*l[0][0],a*l[0][1]))                
                # and isn't chars[g][1] going to give the character value?
                # what do we do with echar then on the VALUE of the character?
                V.append(legendre(a)**k*echar(Z(Chars[g][1])))
            else:
                V.append(0)
        M.append(V)    
    return matrix(M).transpose()

###############################################################################


G = TanakaMonoid(
        p = 5,
        n = 2,
        k = 1,
        delta_prime = 9,
        sigma = 7
);


C = G.get_subgroup_C();

ch_4    = C.get_character(4);
ch_prim = C.get_primitive_characters();

orbits  = G.get_orbits_of_action_on(C);

orbits2  = G.get_orbits_of_action_on_with_char(C, 4);

print(orbits);
print(orbits2);


#print(ch_4);
#print(ch_prim);


exit();


