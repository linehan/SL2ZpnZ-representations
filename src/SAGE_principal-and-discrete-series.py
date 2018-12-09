#!/usr/bin/sage -python
###############################################################################
# ------------------- 
# PROGRAM DESCRIPTION 
# ------------------- 
# This program implements the construction of the complex
# representations of the groups 
#
#   SL_2(Z/(p^n)Z)      (p: odd prime, n: positive integer)
#
# using the method pioneered by Kloosterman in 1946 and 
# extended by Tanaka in 1967. 
# 
# The notation introduced in [KLOO46] and later in [TANA67] 
# is used where appropriate to aid in cross-reference of
# this program with the theory.
#
# A paper [TODO] detailing the aspects of the computation 
# and reviewing the techniques in [KLOO46] and [TANA67]
# is slated to accompany the release of this program.
#
# ---------- 
# REFERENCES 
# ---------- 
# 1. [KLOO46] Kloosterman, H. (1946). The Behavior of General 
# Theta Functions Under the Modular Group and the Characters of 
# Binary Modular Congruence Groups. II. Annals of Mathematics, 
# 47(3), second series, 317-375. doi:10.2307/1969082
#
# 2. [TANA67] Tanaka, S. (1967). Irreducible representations of 
# the binary modular congruence groups mod p^l. J. Math. Kyoto 
# Univ. 7, no. 2, 123-132. doi:10.1215/kjm/1250524272. 
#
###############################################################################
import sys
from sage.all import *

# GLOBAL STATE 
#------------------------------------------------------------------------------
UCF = UniversalCyclotomicField();

# CLASSES
#------------------------------------------------------------------------------
class TanakaParams():
    """
    Validate parameters used to instantiate a TanakaMonoid.
    
    NOTE
    Validation is implemented as a class (rather than a function) 
    to make it simpler for receivers to explicitly verify inputs 
    by testing:

        isinstance(P, TanakaParams).

    rather than verifying each of the individual parameters.

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
    Implements the monoid G = Z/(p^n)Z x Z/(p^(n-k))Z from [TANA67].

    NOTES
    1. This is the main algebraic object used in the construction.
    It is a monoid under its multiplication operation, and its 
    cyclic subgroups (cf. the derived TanakaCyclicSubgroup class) 
    have this multiplication as their group operation. 

    2. For efficiency and compatibility, members are stored as
    pairs of plain Integers, not as pairs of e.g. elements from 
    a SAGE IntegerModRing. This prevents expensive casting during 
    internal arithmetic, and gives the caller a clean expectation
    of what to give and receive.

    So when we say that 
        "u = (u0,u1) is in Z/(p^n)Z x Z/(p^(n-k))Z",
    or
        "u = (u0,u1) is a member of the monoid",
    we mean that u0, u1 are integers such that
        0 <= u0 < p^n 
        0 <= u1 < p^(n-k).
    """
    def __init__(self, params, members=None):
        """
        Instantiate a TanakaMonoid object
        @params : TanakaParams object
        @members: List of Integer tuples in Z/(p^n)Z x Z/(p^(n-k))Z
        Return  : NONE 
        """
        if not isinstance(params, TanakaParams):
            raise TypeError("Expected TanakaParams");

        self.params   = params;
        self.delta    = params.D*(params.p**params.k);
        self.pn       = params.p**params.n;
        self.pnk      = params.p**(params.n-params.k);
        self.identity = (1,0);  # Multiplicative identity.

        if members == None:
            self.members = [
                (u0,u1) for u0 in range(0,self.pn) for u1 in range(0,self.pnk)
            ];
        else:
            # Caller's tuples are coerced to Integer (see NOTES, 2.)
            self.members = [
                (Integer(u0), Integer(u1)) for (u0,u1) in members
            ];

    def __iter__(self):
        return iter(self.members);
    def __str__(self):
        return str(self.members);
    def __getitem__(self, i):
        return self.members[i];

    def order(self): 
        """
        Determine the order of the monoid.
        @NONE
        Return: Integer length of the members array. 
        """
        return len(self.members);

    def norm(self, u):
        """
        Implements the norm defined on the monoid (cf. [TODO]). 
        @u    : Integer tuple (u0,u1) in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer in Z/(p^n)Z
        """
        result = u[0]**2 + self.delta*(u[1]**2);
        return Integer(mod(result, self.pn));
    
    def traceconj(self, u, v):
        """
        Implements the trace conjugate defined on the monoid (cf. [TODO])
        @u    : Integer tuple (u0,u1) in Z/(p^n)Z x Z/(p^(n-k))Z
        @v    : Integer tuple (v0,v1) in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer in Z/(p^n)Z
        """
        result = u[0]*v[0] + self.delta*u[1]*v[1];
        return Integer(mod(result, self.pn));

    def mult(self, u, v):
        """
        Implements the multiplication defined on the monoid (cf. [TODO])
        @u    : Integer tuple (u0,u1) in Z/(p^n)Z x Z/(p^(n-k))Z
        @v    : Integer tuple (v0,v1) in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer tuple (w0,w1) in Z/(p^n)Z x Z/(p^(n-k))Z
        """
        w0 = u[0]*v[0] - self.delta*u[1]*v[1];
        w1 = u[0]*v[1] + u[1]*v[0];

        return (Integer(mod(w0, self.pn)), Integer(mod(w1, self.pnk)));

    def scalar_mult(self, a, u):
        """
        Implements scalar multiplication defined on the monoid.
        @a    : Integer scalar
        @u    : Integer tuple (u0, u1)  in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer tuple (au0,au1) in Z/(p^n)Z x Z/(p^(n-k))Z
        """
        return (Integer(mod(a*u[0], self.pn)), Integer(mod(a*u[1], self.pnk)));

    def cyclic_subgroup(self, members, generator=None):
        """
        Instantiates a TanakaCyclicSubgroup derived from the monoid.
        @members  : List of Integer tuples in Z/(p^n)Z x Z/(p^(n-k))Z
        @generator: Integer tuple (g0, g1) in Z/(p^n)Z x Z/(p^(n-k))Z 
        Return    : TanakaCyclicSubgroup object
        """
        return TanakaCyclicSubgroup(
                params    = self.params,
                members   = members, 
                generator = generator,
        );


class TanakaCyclicSubgroup(TanakaMonoid):
    """
    Implements the subgroups C < G and CL < C from [TANA67]. 

    CAUTION 
    No explicit test for cyclic properties is performed,
    however if the group is not cyclic, you will know soon 
    enough because the method generator() will complain.
    """
    def __init__(self, params, members, generator=None):
        """
        Instantiate a TanakaCyclicSubgroup object.

        @params   : TanakaParams object 
        @members  : List of Integer tuples in Z/(p^n)Z x Z/(p^(n-k))Z
        @generator: Integer tuple (g0, g1) in Z/(p^n)Z x Z/(p^(n-k))Z 
        Return    : NONE
        """
        TanakaMonoid.__init__(self, params=params, members=members);
        self._generator = generator; 
        self._powers    = [];

    def power_of(self, c, d): 
        """
        Determine the power of c that equals d. 
        @c    : Integer tuple (c0,c1) in Z/(p^n)Z x Z/(p^(n-k))Z
        @d    : Integer tuple (d0,d1) in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer pwr smallest such that (@c)^pwr == @d.
        """
        tmp = c;
        pwr = 1;

        while tmp != d:
            tmp = self.mult(tmp, c);
            pwr = pwr + 1;

            if pwr > self.order():
                print("Power of %s equal to %s: %d, order:%d" % (c, d, pwr, self.order()));
                print("[EMATH]: TanakaCyclicSubgroup is not a group");
                
        return pwr; 

    def order_of(self, c):
        """
        Determine order of a group member.
        @c    : Integer tuple (c0,c1) in Z/(p^n)Z x Z/(p^(n-k))Z
        Return: Integer pwr smallest such that (@c)^pwr == self.identity
        """
        return self.power_of(c, self.identity);

    def generator(self):
        """
        Determine a generator for the group
        @NONE
        Return: Integer tuple (c0,c1) in Z/(p^n)Z x Z/(p^(n-k))Z.
        """
        # CACHES
        if self._generator == None:
            for c in self.members:
                if self.order_of(c) == self.order(): 
                    self._generator = c;
                    break;
            else:
                print("[EMATH] TanakaCyclicSubgroup is not cyclic");

        return self._generator;

    def powers(self):
        """
        Provide list associating each member with a power of the generator.
        @NONE
        Return: List of tuples (c,j) where c=self._generator^j. 
        """
        gen = self.generator();

        # CACHES
        if len(self._powers) == 0:
            for c in self:
                self._powers.append((c, self.power_of(gen, c)));
        return self._powers;

    def powers2(self):

        gen = self.generator();

        if len(self._powers) == 0:
            tmp = gen;
            pwr = 1;

            while tmp != c:
                tmp = self.mult(tmp, gen);
                pwr = pwr + 1;
            self._powers.append((c, pwr));
        return self._powers;


class CyclicCharacter():
    """
    Implements a character on a cyclic group.
    
    NOTES
    1. If C is a finite cyclic group of order m, then C is 
    isomorphic to the group of m-th roots of unity under 
    multiplication, and this isomorphism is the character.

    The primitive m-th roots of unity are the complex numbers
        e^(2*pi*i*j*(1/m))      for 1<=j<=m, gcd(j,m)=1.
    Let 
        zeta_m = e^(2*pi*i*(1/m)). 
        
    Then the m-th roots of unity are the distinct powers 
        zeta_m^j                for 1<=j<=m.

    Fix a generator <g>=C. The character must map g to 
    a primitive m-th root of unity. There are phi(m) many
    primitive m-th roots of unity, where phi() is Euler's
    totient function, giving us phi(m) distinct characters
    on C.

    Define the character 
        X_j:C-->Complex         for 1<=j<=m, gcd(j,m)=1
            g-->zeta_m^j

    Once we have defined where the generator goes, the rest
    of the elements follow by multiplication:

        X_j(u=g^k) = (zeta_m^j)^k.

    2. These facts mean that to specify a character we need
    only specify the integer order (m) of the group and the
    integer power (j) of zeta_m to which we send the fixed 
    generator of C. 
    """
    def __init__(self, order, power=1, exact=True):
        """
        Instantiate a CyclicCharacter object.
        @order: Integer specifying primitive @order-th root of unity.
        @power: Integer power of zeta_@order (selects generator/primitive RoU).
        @exact: True-Store as Cyclotomic polynomial; False-Store as Complex.
        Return: NONE

        CAUTION
        If gcd(@order, @power) != 1, then zeta_@order^@power 
        is not primitive and hence won't generate the group.
        """
        self.base = UCF.zeta(order);
        self.power = power;
        self.order = order;
        if exact == False:
            self.base = complex(self.base.real(), self.base.imag());

    def eval(self, j):
        """
        Implement the mapping u=g^j --> (self.base)^j.
        @j    : Integer power to take of the base.
        Return: Cyclotomic polynomial or Complex value (see @exact, __init__())

        """
        return self.base**mod(self.power*Integer(j), self.order);


#------------------------------------------------------------------------------
# Below this is work in progress -Jason
#------------------------------------------------------------------------------

def legendre(a, p):
    """
    Computes the value of the Legendre symbol (-/p).
    @a    : Integer numerator
    @p    : Denominator of the symbol
    Return: Integer in {0, -1, 1} 
    """
    val = mod(a, p)**((p-1)/2);
    return -1 if val == -1 else Integer(val); 


def compute_W_constant(p, n, k, D, s):
    """ 
    Computes the complicated constant used to make W 
    """
    # A fancy way to write i
    i = UCF.zeta(4); 

    # Compute epsilon
    if   k%2 == 1 and legendre(-1, p) == 1:
        epsilon = 1;
    elif k%2 == 1 and legendre(-1, p) == -1:
        epsilon = -i;
    else:
        epsilon = -1**n;

    # Compute square root of p 
    zeta_p = UCF.zeta(p);
    sqrt_p = sum(legendre(a, p)*(zeta_p**a) for a in range(1, p));
    if p % 4 == 3:
        sqrt_p = i * sqrt_p; 

    # The final constant 
    return sqrt_p**(k -2*n)          \
         * legendre(D,p)**(n-k)      \
         * legendre(s,p)**k*epsilon; 


def get_orbits(G, C):
    """
    Build list of class reps. of distinct orbits of action of @C on @G
    @G    : TanakaMonoid object
    @C    : TanakaCyclicSubgroup of @G
    Return: List of Integer tuples (u0,u1) in Z/(p^n)Z x Z/(p^(n-k))Z
    """
    ## TODO: Revise this to be
    ## 1. in-place rather than requiring copy
    ## 2. use the powers of generator rather than 'X'
    ## 3. Re-name variables to be clear
    orbit_reps = [];
    G_members  = list(G.members); 
    X          = dict(C.powers());

    while len(G_members) != 0:
        g          = G_members[0];
        orbit_of_g = {};
        for c in C:
            gc = C.mult(g, c);
            if gc not in orbit_of_g:
                orbit_of_g[gc] = c; 
            else: 
                d = orbit_of_g[gc]; # gc == gd
                if X[c] != X[d]:
                    # NOTE:
                    # This leaves orbit_of_g incomplete, so when 
                    # we remove orbit_of_g from G_members below,
                    # there may be *something* there that needs to
                    # get removed. 
                    #
                    # If there is nothing to be removed... then we
                    # never halt...?
                    break;
        else:
            orbit_reps.append(C.mult(g, C.members[0]));

        for u in orbit_of_g:
            if u in G_members:
                G_members.remove(u);

    return orbit_reps;



class TanakaRepSpace():
    """
    Implement the representation space R_k(D, s).
    """
    def __init__(self, p, n, k, D, s, exact=False):
        params = self.params = TanakaParams(p, n, k, D, s);

        G = self.G = TanakaMonoid(params);
        C = self.C = G.cyclic_subgroup(
            members = [
                g for g in G if G.norm(g) == 1
            ]
        );
        CL = self.CL = C.cyclic_subgroup(
            generator = C.generator(), 
            members   = [
                # TODO: Is |CL| always == p? I have a sneaking suspicion we 
                # don't actually *need* this guy, but need proof to discard.
                c for c in C if 
                    (mod(c[0], p**(n-1)), mod(c[1], p**(n-1-k))) == C.identity
            ],
        );

        self.e_sigma = CyclicCharacter(order=p**n, power=s, exact=exact);

        # W stuff
        self.W_constant = compute_W_constant(p, n, k, D, s);
        self.W_cached   = None;

        # Get orbit representatives
        self.orbit_reps = get_orbits(self.G, self.C);

        # For setting the primitive character.
        self.exact = exact;

    def get_primitive_characters(self):
        """
        Determine which characters are primitive 
        Returns list of character indices, NOT TanakaCharacter-s
        """
        primitives = [];
        for chi in range(0, self.C.order()):
            if not all(mod(xc*chi, self.C.order()) == 0 for (c, xc) in self.CL.powers()):
                primitives.append(chi);
        return primitives;

    def set_primitive_character(self, x):
        """
        @x    :
        Return:
        """
        self.X = CyclicCharacter(
            order=self.C.order(),
            power=x,
            exact=self.exact
        );
        self.W_cached = None; # Reset the cache

    def A(self, a):
        """ 
        Compute the action A
        @a    : Integer in
        Return: Matrix
        """
        if self.X == None:
            print("[WARN] You must first call set_primitive_character()");

        M = [];
        p = self.params.p;
        k = self.params.k;
        for o1 in self.orbit_reps:
            V   = [];
            ao1 = self.C.scalar_mult(a, o1);
            for o2 in self.orbit_reps:
                for (c, xc) in self.C.powers():
                    if self.G.mult(c, o2) == ao1:
                        V.append(
                            (legendre(a,self.params.p)**self.params.k)*self.X.eval(xc)
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
        if self.X == None:
            print("[WARN] You must first call set_primitive_character()");

        if self.W_cached == None:
            M = [];
            for o1 in self.orbit_reps:        
                V = [];
                H = [(xc, self.G.mult(c, o1)) for (c,xc) in self.C.powers()];
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
        ring = IntegerModRing(self.params.p**self.params.n);

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


# CALLER INTERFACE
#------------------------------------------------------------------------------

def get_representations_of(a, b, c, d, p, n):
    """
    This will ultimately be one of the main interfaces.
    Will do something like below.
    """
    pass;


# TESTING
#------------------------------------------------------------------------------


T = TanakaSystem(p=3, n=2, exact=True);
R = T.representation_space(k=0, D=1, s=1);

print("|C|:"+str(R.C.order()));

if R.C.powers() == R.C.powers2():
    print("Ok!");

X = R.get_primitive_characters();
R.set_primitive_character(X[0]);


print(R.rep(1, 0, 0, 1));

