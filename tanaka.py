#!/usr/bin/sage -python
"""
 INTRODUCTION
 ------------ 
 This program implements the construction of the complex
 representations of the groups 

   SL_2(Z/(p^n)Z)      (p: odd prime, n: positive integer)

 using the method pioneered by Kloosterman in 1946 and 
 extended by Tanaka in 1967. 
 
 The notation introduced in [KLOO46] and later in [TANA67] 
 is used where appropriate to aid in cross-reference of
 this program with the theory.

 A paper [TODO] detailing the aspects of the computation 
 and reviewing the techniques in [KLOO46] and [TANA67]
 is slated to accompany the release of this program.

 References
 ----------
 1. [KLOO46] Kloosterman, H. (1946). The Behavior of General Theta Functions 
 Under the Modular Group and the Characters of Binary Modular Congruence 
 Groups. II. Annals of Mathematics, 47(3), second series, 317-375. 
 doi:10.2307/1969082

 2. [TANA67] Tanaka, S. (1967). Irreducible representations of the binary 
 modular congruence groups mod p^l. J. Math. Kyoto Univ. 7, no. 2, 123-132. 
 doi:10.1215/kjm/1250524272. 
"""

import sys
from sage.all import *

# GLOBAL STATE 
###############################################################################
UCF = UniversalCyclotomicField();

# CLASSES
###############################################################################
class TanakaParams():
    """
    Validate parameters used to instantiate Tanaka objects.
    
    NOTE
    Validation is implemented as a class (rather than a function) to make
    it faster/simpler for receivers to explicitly verify inputs by testing 

        isinstance(P, TanakaParams).

    """
    def __init__(self, p, n, k, D, s):
        """
        @p: odd prime
        @n: the power that @p will be raised to
        @k: TODO
        @D: TODO
        @s: TODO
        """
        assert isinstance(p, int),   "Parameter 'p' must be an integer";
        assert isinstance(n, int),   "Parameter 'n' must be an integer";
        assert isinstance(k, int),   "Parameter 'k' must be an integer";

        assert p > 0,                "Parameter 'p' must be positive";
        assert p != 2,               "Parameter 'p' must be odd";
        assert n >= 1,               "Parameter 'n' must be > 1";
        assert k >= 0,               "Parameter 'k' must be >= 0";
        assert k <= n,               "Parameter 'k' must be <= n";

        # is_prime() is from SAGE.
        assert is_prime(p),          "Parameter 'p' must be prime!"; 
        assert s % p != 0,           "Parameter 's' must be != 0 (mod p)";
        assert D % p != 0,           "Parameter 'D' must be != 0 (mod p)";

        self.p = p;
        self.n = n;
        self.k = k;
        self.D = D;
        self.s = s;


class TanakaMonoid(): 
    """
    Implement the monoid G = Z/(p^n)Z x Z/(p^(n-k))Z introduced in [TANA67].
    """
    def __init__(self, params, members=None):
        """
        Instantiate a TanakaMonoid object
        @params : TanakaParams object
        @members: List of tuples (a,b) in Z/(p^n)Z x Z/(p^(n-k))Z.
        """
        if not isinstance(params, TanakaParams):
            raise TypeError("Expected TanakaParams");

        self.params = params;
        self.delta  = (params.p**params.k)*params.D;

        self.Zpn  = IntegerModRing(params.p**params.n);
        self.Zpnk = IntegerModRing(params.p**(params.n-params.k));

        if members == None:
            self.members = [(a,b) for a in self.Zpn for b in self.Zpnk];
        else:
            self.members = members;

    def __iter__(self):
        return iter(self.members);
    def __str__(self):
        return str(self.members);
    def __getitem__(self, i):
        return self.members[i];

    def order(self): 
        return len(self.members);

    def norm(self, g):
        """
        Implements the norm map defined on the monoid. 
        @g    : Member of the monoid Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Value in Z/(p^n)Z
        """
        x = ZZ(g[0])**2 + self.delta*ZZ(g[1])**2;
        return self.Zpn(x);
    
    def traceconj(self, u, v):
        """
        Implements the trace conjugate defined on the monoid.
        @u    : Member of the monoid
        @v    : Member of the monoid
        Return: Value in Z/(p^n)Z
        """
        x = ZZ(u[0])*ZZ(v[0]) + self.delta*ZZ(u[1])*ZZ(v[1]);
        return ZZ(self.Zpn(x));

    def mult(self, u, v):
        """
        Implements the special multiplication defined on the monoid.
        @u    : Member of the monoid
        @v    : Member of the monoid
        Return: Member of the monoid
        """
        u1 = ZZ(u[0]); u2 = ZZ(u[1]);
        v1 = ZZ(v[0]); v2 = ZZ(v[1]);

        w1 = u1*v1 - self.delta*u2*v2;
        w2 = u1*v2 + u2*v1;

        # TODO: This is casting to a kind of "TanakaMonoidMember"
        return (self.Zpn(w1), self.Zpnk(w2));

    def cyclic_abelian_subgroup(self, members, generator=None):
        """
        Instantiate a subgroup derived from the monoid.
        @members  : A list of members of the subgroup derived from self.members 
        @generator: A generating member, to ensure consistent indexing 
        Return    : TanakaCyclicAbelianSubgroup object.
        """
        return TanakaCyclicAbelianSubgroup(
                params    = self.params,
                members   = members, 
                generator = generator
        );


class TanakaCyclicAbelianSubgroup(TanakaMonoid):
    """
    Implements the subgroups C and CL from the text. 
        
    CAUTION 
    No explicit test for cyclic or abelian properties is 
    performed! If the group is not cyclic, its methods 
    will fail to halt. 
    """
    def __init__(self, params, members, generator=None):
        """
        Instantiate a TanakaCyclicAbelianSubgroup object.
        @params   : TanakaParams object 
        @members  : List of members derived from the parent. 
        @generator: Generator in @members, to ensure consistent indexing.
        """
        TanakaMonoid.__init__(self, params=params, members=members);
        self._generator = generator;

    def order_of(self, member):
        """
        Determine the order of a member of the group.
        @member: An element from self.members
        Return : Integer order of that element
        """
        g = member;
        o = 1;
        while g != (self.Zpn(1), self.Zpnk(0)):
            g = self.mult(g, member);
            o = o + 1;
        return o; 

    def generator(self):
        """
        Select a generator for the group.
        Return: A member of the group which generates it.
        """
        if self._generator == None:
            t = 0 
            while self.order_of(self.members[t]) != self.order():
                t = t+1
            self._generator = self.members[t];
        return self._generator;


class CyclicCharacter():
    """
    Implement a character of the form e^(2*pi*i*(1/n)*j).

    NOTES
    If C is a cyclic group of order m, then any character on C
    must correspond to one of the m many primitive m-th roots of 
    unity, each given as an integer power of e^(2*pi*i(1/|C|)).
    Then we can simply index each character by providing an integer
    corresponding to this integer power.
    """
    def __init__(self, n, power=1, exact=True):
        self.base = UCF.zeta(n)**power;
        if exact == False:
            self.base = complex(self.base.real(), self.base.imag());

    def eval(self, j):
        return self.base**ZZ(j);


class TanakaCharacter(CyclicCharacter):
    """
    Implement a character on a TanakaCyclicAbelianSubgroup.

    NOTES
    The domain C of a character X is always finite, so 
    the character is stored as a list of tuples

        (c, xc) in C x Z,

    where c is the element of C, and xc, the image of c under X,
    is an integer value from which the actual complex value of 
    the character can be recovered (avoids loss of precision 
    during intermediate steps). 
    
    If C=<g>, c=g^j is an element of C, then xc = j. 
    """
    def __init__(self, C, x, exact=True):
        """
        Instantiate a TanakaCharacter object.
        @C: TanakaCyclicAbelianSubgroup object.
        @x: Integer index of the character. 
        """
        if not isinstance(C, TanakaCyclicAbelianSubgroup):
            raise TypeError("Domain must be TanakaCyclicAbelianSubgroup");

        CyclicCharacter.__init__(self, C.order(), exact=exact);

        data  = [];
        gen   = C.generator();
        order = C.order();

        for c in C:
            tmp = gen;
            pwr   = 1;
            while tmp != c:
                tmp = C.mult(tmp, gen);
                pwr = pwr + 1;
            data.append((c, mod(x*pwr, order)));

        self.data = data;

    def __iter__(self):
        return iter(self.data);
    def __str__(self):
        return str(self.data);
    def __getitem__(self, i):
        return self.data[i];


class TanakaRepSpace():
    """
    Implement the representation space R_k(D, s).
    """
    def __init__(self, p, n, k, D, s, exact=False):
        params = TanakaParams(p, n, k, D, s);

        # The monoid Z/p^nZ x Z/p^(n-k)Z
        G = TanakaMonoid(params);

        # The subgroup C < G
        C = G.cyclic_abelian_subgroup(
            members = [
                g for g in G if G.norm(g) == 1
            ]
        );

        # The subgroup CL < C
        CL = C.cyclic_abelian_subgroup(
            generator = C.generator(), 
            members   = [
                # Note the content of CL is dependent on k
                c for c in C if 
                    (mod(c[0], p**(n-1)), mod(c[1], p**(n-1-k))) == (1,0)
            ],
        );

        # Attach algebraic structures to TanakaRepSpace 
        # after instantiation (allows reader to copy the
        # above initialization without editing out self.*).
        self.params = params;
        self.G      = G;
        self.C      = C;
        self.CL     = CL;
        self.exact  = exact;

        self.e_sigma = CyclicCharacter(p**n, power=s, exact=exact);

        # Compute W stuff
        self.W_constant = self._compute_W_constant(p, n, k, D, s);
        self.W_cached   = None;

    def _legendre(self, a):
        """
        Computes the value of the Legendre symbol (-/p).
        @a    : Integer numerator
        Return: Integer in {0, -1, 1} 
        """
        r = mod(a, self.params.p)**((self.params.p-1)/2);
        return Integer(r) if r != -1 else -1;

    def _compute_W_constant(self, p, n, k, D, s):
        """ 
        Computes the complicated constant used to make W 
        """
        # A fancy way to write i
        i = UCF.zeta(4); 

        # Compute epsilon
        if   k%2 == 1 and self._legendre(-1) == 1:
            epsilon = 1;
        elif k%2 == 1 and self._legendre(-1) == -1:
            epsilon = -i;
        else:
            epsilon = -1**n;

        # Compute square root of p 
        zeta_p = UCF.zeta(p);
        sqrt_p = sum(self._legendre(a)*(zeta_p**a) for a in range(1, p));
        if p % 4 == 3:
            sqrt_p = i * sqrt_p; 

        # The final constant 
        return sqrt_p**(k -2*n)              \
             * self._legendre(D)**(n-k)      \
             * self._legendre(s)**k*epsilon; 

    def _get_orbits(self):
        """
        Computes list of representatives for orbit of action of C on G.
        """
        orbit_reps = [];
        G_members  = list(self.G.members); # Copying this sucks 
        X          = dict(self.X);

        while len(G_members) != 0:
            g          = G_members[0];
            orbit_of_g = {};
            for c in self.C:
                gc = self.G.mult(g, c);
                if gc not in orbit_of_g:
                    orbit_of_g[gc] = c; 
                else: 
                    d = orbit_of_g[gc]; # gc == gd
                    if X[c] != X[d]:
                        break;
            else:
                orbit_reps.append(self.G.mult(g, self.C[0]));

            for u in orbit_of_g:
                G_members.remove(u);

        return orbit_reps;

    def get_primitive_characters(self):
        """
        Returns list of character indices, NOT TanakaCharacter-s
        """
        primitives = [];
        for x in range(0, self.C.order()):
            if not all(cx == 0 for (c, cx) in TanakaCharacter(self.CL, x)):
                primitives.append(x);
        return primitives;

    def set_primitive_character(self, chi):
        """
        @chi:
        Return:
        """
        self.X          = TanakaCharacter(self.C, chi, exact=self.exact);
        self.orbit_reps = self._get_orbits();
        self.W_cached   = None; # Reset the cache

    def A(self, a):
        """ 
        Compute the action A
        @a    : Integer in
        Return: Matrix
        """
        M = [];
        for o1 in self.orbit_reps:
            V   = [];
            ao1 = (a*o1[0], a*o1[1]);
            for o2 in self.orbit_reps:
                for (c, xc) in self.X:
                    if self.G.mult(c, o2) == ao1:
                        V.append(
                            (self._legendre(a)**self.params.k)*self.X.eval(xc)
                        );
                        break;
                else:
                    V.append(0);
            M.append(V);
        return matrix(M).transpose()

    def B(self, b):
        """ 
        Compute the action B 
        @b    : Integer in 
        Return: Matrix
        """
        M = [];
        for o1 in self.orbit_reps:
            V = []
            for o2 in self.orbit_reps:            
                if o1 == o2: # Along diagonal
                    V.append(self.e_sigma.eval(ZZ(b * self.G.norm(o2))));
                else:
                    V.append(0);
            M.append(V);
        return matrix(M)

    def W(self):
        """ 
        Compute the action W 
        Return: A matrix 
        NOTE
        The result is cached.
        """
        if self.W_cached == None:
            M = [];
            for o1 in self.orbit_reps:        
                V = [];
                H = [(xc, self.G.mult(c, o1)) for (c,xc) in self.X];
                for o2 in self.orbit_reps:
                    Sum = 0
                    for (xc, co1) in H:
                        Sum += self.e_sigma.eval(-2*self.G.traceconj(o2,co1)) \
                             * self.X.eval(xc);
                    V.append(Sum);
                M.append(V)
            self.W_cached = self.W_constant*matrix(M);
        return self.W_cached;

    def rep(self, a, b, c, d):
        """ 
        Implements Bruhat decomposition 
        @a    :
        @b    :
        @c    :
        @d    :
        Return:
        """
        ring = Integers(self.params.p**self.params.n);

        a = ring(a);
        b = ring(b);
        c = ring(c);
        d = ring(d);

        if c % self.params.p == 0:
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
    def __init__(self, p, n, exact=False):
        self.exact = exact;
        self.p     = p;
        self.n     = n;

    def representation_space(self, k, D, s):
        return TanakaRepSpace(self.p, self.n, k, D, s, exact=self.exact);


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
#T = TanakaSystem(p=17, n=2);
#R = T.representation_space(k=1, D=1, s=1);

T = TanakaSystem(p=5, n=2, exact=True);
R = T.representation_space(k=1, D=9, s=7);

X = R.get_primitive_characters();

R.set_primitive_character(4);


from jason_sage import test_output;
from jason_sage import test_primitives;

out = test_output();

p1 = test_primitives();
p2 = R.get_primitive_characters();

print(p1 == p2);




print(R.B(3) == out[0]);
print((R.B(3) - out[0]).norm());
print(R.A(4) ==  out[1]);
print((R.A(4) - out[1]).norm());
print(R.W() ==  out[2]);
print((R.W() - out[2]).norm());

exit();


#print(R.rep(1, 3, 0, -1) == R.B(3));

#for x in X:
    #R.set_primitive_character(x);
    #print(R.B(3));

#T = TanakaSystem(p=5, n=2, exact=True);
#R = T.representation_space(k=1, D=9, s=7);
#X = R.get_primitive_characters();


#R.set_primitive_character(4);

#print(R.W());
#exit();

#from jason_base_no_basis import getit;

#B3  = R.A(3);
#B3p = getit();

#print(B3)
#print(B3p)

#from jason_sage import test_output;
#out = test_output();

#print(R.B(3) == out[0]);
#print((R.B(3) - out[0]).norm());
#print(R.A(4) ==  out[1]);
#print((R.A(4) - out[1]).norm());
#print(R.W() ==  out[2]);
#print((R.W() - out[2]).norm());
