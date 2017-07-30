#!/usr/bin/sage -python

import sys
from sage.all import *


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
        return TanakaSubgroup(
            monoid=self, 
            subset=[g for g in self.members if self.norm(g) == 1]
        );

###############################################################################
def get_chi(j, G, generator):
    """
    If G = <g> is a multiplicative group, 
    the characters are defined:

        chi_j(g^k) = exp(2*pi*i*(1/|G|))^{j*k}

    where 0 <= j <= |G|-1.

    ----

    We don't represent characters as functions in this program. 

    Each character is indexed by an integer j (0 <= j <= |G|-1).

    The value of a character on an element u = g^k in G is 
    determined by the fixed base B = exp(2*pi*i*(1/|G|)), 
    raised to the power j*k. 

    Hence to store character values, we only store the value j*k,
    adding the base later.

    We can then represent each character as a list of tuples:

        (u, jk),

    where u = g^k is an element of G, and j is the character index. 

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


class TanakaCyclicAbelianSubgroup(TanakaMonoid):
    """ 
    Used for the subgroups CL < C < G 
    """

    def __init__(self, monoid, subset, generator=None):
        TanakaMonoid.__init__(self, 
            p           = monoid.param.p, 
            n           = monoid.param.n, 
            k           = monoid.param.k, 
            sigma       = monoid.param.sigma, 
            delta_prime = monoid.param.delta_prime
        );

        self.monoid               = monoid;
        self.members              = subset;

        self.generator            = generator;
        self.characters           = {};
        self.primitive_characters = [];

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
            self.characters[chi] = get_chi(chi, self, self.get_generator());

        return self.characters[chi];


    def get_primitive_characters(self):
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



G = TanakaMonoid(
        p = 5,
        n = 2,
        k = 1,
        delta_prime = 9,
        sigma = 7
);


for u in G:
    print(u);


C = G.get_subgroup_C();

cc = C.get_character(4);

ch = C.get_primitive_characters();

print(ch);
print(cc);


exit();







            #g = self.get_generator();
            #for j in self.members:
                #t = g;
                #a = self.R0(chi);
                #while t != j:
                    #t = self.mult(t, g);
                    #a = a + chi; 
                #self.characters[chi].append((j, a));

        #return self.characters[chi];

    #def get_primitive_characters(self):
        #if len(self.primitive_characters) == 0:
            ## Though we are using CL < C here, any generator for C
            ## is also a generator for CL, so the generator can be
            ## the same.
            #g = self.get_generator();
            #for chi in self.R0:
                #B = [];
                #for j in self.CL:
                    #t = g
                    #a = chi;
                    #while t != j:
                        #t = self.mult(t, g);
                        #a = a + chi; 
                    #B.append((j, a));

                #if all(b[1] == 0 for b in B) == False:
                    #self.primitive_characters.append(chi);

        #return self.primitive_characters;
