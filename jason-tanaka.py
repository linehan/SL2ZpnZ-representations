#!/usr/bin/sage -python

import sys
from sage.all import *

###############################################################################
# BACKEND CLASSES
###############################################################################

class TanakaParams():
    """
    Validate parameters used to instantiate Tanaka objects.
    """
    def __init__(self, p, n, k, delta_prime, sigma):
        assert isinstance(p, int),   "Parameter 'p' must be an integer";
        assert isinstance(n, int),   "Parameter 'n' must be an integer";
        assert isinstance(k, int),   "Parameter 'k' must be an integer";

        assert p > 0,                "Parameter 'p' must be positive";
        assert n > 1,                "Parameter 'n' must be > 1";
        assert k >= 0,               "Parameter 'k' must be >= 0";
        assert k <= n,               "Parameter 'k' must be <= n";

        # is_prime() is from SAGE.
        assert is_prime(p),          "Parameter 'p' must be prime!"; 
        assert sigma % p != 0,       "Parameter 'sigma' must be != 0 (mod p)";

        assert delta_prime % p != 0, "Parameter 'delta_prime' must be != 0 (mod p)";

        self.p           = p;
        self.n           = n;
        self.k           = k;
        self.delta_prime = delta_prime;
        self.sigma       = sigma;


class TanakaMonoid(): 
    """
    Implements the monoid G from the text.
    """
    def __init__(self, params, members=None):
        if not isinstance(params, TanakaParams):
            raise TypeError("Expected TanakaParams");

        self.param = params;
        self.delta = (self.param.p**self.param.k)*self.param.delta_prime;

        self.Rn   = IntegerModRing(self.param.p**self.param.n);
        self.Rn_k = IntegerModRing(self.param.p**(self.param.n-self.param.k));

        if members == None:
            self.members = [ (a,b) for a in self.Rn for b in self.Rn_k ];
        else:
            self.members = members;

        # Provides indexes for each of the elements of this structure,
        # and also for elements (characters) of its associated character
        # group, when isomorphic (e.g. in TanakaCyclicAbelianSubgroup).
        self.R0 = Integers(self.order());
        
    def __iter__(self):
        return iter(self.members);
    def __str__(self):
        return str(self.members);
    def __repr__(self):
        pass;
    def __getitem__(self, i):
        return self.members[i];

    def order(self): 
        return len(self.members);

    def norm(self, u):
        x = ZZ(u[0])**2 + self.delta*ZZ(u[1])**2;
        return self.Rn(x);
    
    def traceconj(self, u, v):
        x = ZZ(u[0])*ZZ(v[0]) + self.delta*ZZ(u[1])*ZZ(v[1]);
        return self.Rn(x);

    def mult(self, u, v):
        w1 = ZZ(u[0])*ZZ(v[0]) - self.delta*ZZ(u[1])*ZZ(v[1]);
        w2 = ZZ(u[0])*ZZ(v[1]) + ZZ(u[1])*ZZ(v[0]);
        w  = (self.Rn(w1), self.Rn_k(w2));
        return w;

    def traceconj_fast(self, u, v):
        x = (u[0])*(v[0]) + self.delta*(u[1])*(v[1]);
        return self.Rn(x);

    def mult_fast(self, u, v):
        w1 = (u[0])*(v[0]) - self.delta*(u[1])*(v[1]);
        w2 = (u[0])*(v[1]) + (u[1])*(v[0]);
        w  = (self.Rn(w1), self.Rn_k(w2));
        return w;
        

class TanakaCyclicAbelianSubgroup(TanakaMonoid):
    """
    Implements the subgroups C and CL from the text. 
    """
    def __init__(self, parent, members, generator=None):
        TanakaMonoid.__init__(self, params=parent.param, members=members);
        self.parent = parent;
        self.gen    = generator;

    def order_of(self, element):
        g = element;
        o = 1;
        while g != (self.Rn(1), self.Rn_k(0)):
            g = self.mult(g, element);
            o = o + 1;
        return o; 

    def generator(self):
        if self.gen == None:
            t = 0 
            while self.order_of(self.members[t]) != self.order():
                t = t+1
            self.gen = self.members[t];
        return self.gen;


class TanakaCharacter():
    """
    A TanakaCharacter is a list of tuples in a certain format.
    """
    def __init__(self, domain, x):
        if not isinstance(domain, TanakaCyclicAbelianSubgroup):
            raise TypeError("Domain must be TanakaCyclicAbelianSubgroup");

        data = [];
        gen  = domain.generator();

        for g in domain:
            h = gen;
            a = domain.R0(x); # Why do we cast with R0 here? To truncate?
            while h != g:
                h = domain.mult(h, gen);
                a = a + x;
            data.append((g, a));

        self.data = data;

    def __iter__(self):
        return iter(self.data);
    def __str__(self):
        return str(self.data);
    def __getitem__(self, i):
        return self.data[i];


class TanakaRepSpace():
    """
    This class implements R_k(delta_prime, sigma) from the text. 
    """
    def __init__(self, p, n, k, delta_prime, sigma):
        self.param = params = TanakaParams(p, n, k, delta_prime, sigma);

        R1 = Integers(p**(n-1));
        R2 = Integers(p**(n-1-k));

        self.G = TanakaMonoid(params);

        self.C = TanakaCyclicAbelianSubgroup(
            parent  = self.G, 
            members = [g for g in self.G if self.G.norm(g) == 1]
        );

        self.CL = TanakaCyclicAbelianSubgroup(
            parent    = self.G,
            generator = self.C.generator(), # is this generator in CL?
            members = [g for g in self.C if R1(g[0])==R1(1) and R2(g[1])==R2(0)],
        );

        self.F  = UniversalCyclotomicField();
        self.Rp = Integers(p); 

        self.W_constant = self._W_constant();
        self.W_cached   = None;

        self.base_esigma = self.F.zeta(self.param.p**self.param.n)**(self.param.sigma);
        self.base_echi   = self.F.zeta(self.C.order());

    def _W_constant(self):
        i = self.F.zeta(4); 
        if self._legendre(-1) == 1 and self.param.k%2 == 1:
            epsilon = 1;
        elif self._legendre(-1) == -1 and self.param.k%2 == 1:
            epsilon = -i;
        else:
            epsilon = -1**n;

        tmp_sum = sum(self._legendre(a)*self.F.zeta(self.param.p,a) for a in range(1,self.param.p));
        if self.param.p % 4 == 3:
            sqrt_p = i * tmp_sum; 
        else:
            sqrt_p = tmp_sum;

        # The mighty constant
        return sqrt_p**(self.param.k-2*self.param.n) \
             * self._legendre(self.param.delta_prime)**(self.param.n-self.param.k) \
             * self._legendre(self.param.sigma)**self.param.k*epsilon; 

    def _legendre(self, a):
        b = self.Rp(a)**((self.param.p-1)/2)
        if b == -1:
            return -1;
        else:
            return ZZ(b);

    def _esigma(self, a):
        return self.base_esigma**a; 

    def _echi(self, a):
        return self.base_echi**a;

    def _get_orbits(self):
        orbit_reps = [];
        G_members  = list(self.G.members); # Copying this isn't my fav. thing.. 
        X          = dict(self.X);         # Dict works better than tuple here. 

        while len(G_members) != 0:
            g          = G_members[0];
            orbit_of_g = {};
            degenerate = False;

            for c in self.C:
                gc = self.G.mult(g, c);
                if gc not in orbit_of_g:
                    orbit_of_g[gc] = c; 
                else: 
                    d = orbit_of_g[gc]; # gc == gd
                    if X[c] != X[d]:
                        degenerate = True;
                        break;

            if not degenerate:
                orbit_reps.append(self.G.mult(g, self.C[0]));

            for u in orbit_of_g:
                G_members.remove(u);

        return orbit_reps;

    def set_primitive_character(self, chi):
        self.X          = TanakaCharacter(self.C, chi);
        self.orbit_reps = self._get_orbits();

    def get_primitive_characters(self):
        """ Note returns list of character indices, NOT TanakaCharacter-s """
        primitives = [];
        for x in self.C.R0:
            if not all(cx == 0 for (c, cx) in TanakaCharacter(self.CL, x)):
                primitives.append(x);
        return primitives;

    def A(self, a):
        M = [];
        for o1 in self.orbit_reps:
            V = [];
            for o2 in self.orbit_reps:
                H = [self.G.mult(c, o2) for (c, xc) in self.X];
                # TODO: Put this scalar mult on a class ?
                if (a*o1[0], a*o1[1]) in H:
                    g   = H.index((a*o1[0], a*o1[1]));                
                    val = self._legendre(a)**self.param.k*self._echi((self.X[g][1]));
                    V.append(val);
                else:
                    V.append(0);
            M.append(V);
        return matrix(M).transpose()

    def B(self, b):
        M = [];
        for o1 in self.orbit_reps:
            V = []
            for o2 in self.orbit_reps:            
                if o1 == o2:
                    modulus = b * self.G.norm(o2);
                    V.append(self._esigma(ZZ(modulus)))
                else:
                    V.append(0)
            M.append(V)
        return matrix(M)

    def W(self):
        if self.W_cached == None:
            M = [];
            for o1 in self.orbit_reps:        
                H = [self.G.mult_fast(c, o1) for (c, cx) in self.X];
                V = [];
                for o2 in self.orbit_reps:
                    Sum1 = 0
                    for h in H:
                        g = H.index(h)
                        # pull things out of this sum
                        Value = self._esigma(-2*(self.G.traceconj_fast(o2,h))) \
                              * self._echi((self.X[g][1]));
                        Sum1 = Value + Sum1;
                    V.append(Sum1*self.W_constant)
                M.append(V)
            self.W_cached = matrix(M);
        return self.W_cached;

    def rep(self, a, b, c, d):
        """ Implements Bruhat decomposition """
        ring = Integers(self.param.p**self.param.n);

        a = ring(a);
        b = ring(b);
        c = ring(c);
        d = ring(d);

        if c % self.param.p == 0:
            bruhat = self.W()               \
                   * self.B(-c*(a)**(-1))   \
                   * self.W()               \
                   * self.A(-a)             \
                   * self.B(b*a**(-1));
        else:
            bruhat = self.B(a*c**(-1))      \
                   * self.W()               \
                   * self.A(c)              \
                   * self.B(d*c**(-1));
        
        return bruhat; 


class TanakaSystem():
    def __init__(self, p, n):
        self.p = p;
        self.n = n;

    def representation_space(self, k, delta_prime, sigma):
        return TanakaRepSpace(self.p, self.n, k, delta_prime, sigma);


###############################################################################
# CALLER INTERFACE
###############################################################################

def get_representations_of(a, b, c, d, p, n):
    """
    This will ultimately be one of the main interfaces.
    Will do something like below.
    """
    pass;


###############################################################################
# TESTING
###############################################################################
#T = TanakaSystem(p = 5, n = 2);
#R = T.representation_space(k=1, delta_prime=9, sigma=7);

#T = TanakaSystem(p = 17, n = 2);
#R = T.representation_space(k=1, delta_prime=1, sigma=1);

#X = R.get_primitive_characters();

#R.set_primitive_character(X[0]);

#print(R.W());

T = TanakaSystem(p=5, n=2);
R = T.representation_space(k=1, delta_prime=9, sigma=7);

X = R.get_primitive_characters();

R.set_primitive_character(4);

#from jason_sage import test_output

#out = test_output();

#print(R.B(3) == out[0]);
#print(R.A(4) == out[1]);
#print(R.W()  == out[2]);

R.B(3);
R.A(4);
R.W();


#print(R.rep(1, 3, 0, -1));
#print(R.B(3));
#print(R.rep(1, 3, 0, -1) == R.B(3));

#for x in X:
    #R.set_primitive_character(x);
    #print(R.B(3));

