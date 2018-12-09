
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
    G_members  = G.members; # Copying this isn't efficient 

    while len(G) != 0:
        g = G_members[0];
        g_orbit = {};
        g_orbit_rep = None;
        
        for c in C:
            gc = G.mult(g, c);

            if g_orbit_rep is None:
                g_orbit_rep = gc;

            if gc not in g_orbit:
                g_orbit[gc] = c; 
            else: 
                d = g_orbit[gc]; # gc == gd
                if C.get_character_value(x, c) != C.get_character_value(x, d):
                    g_orbit_rep = None;
                    break;

        if rep != None:
            orbit_reps.append(g_orbit_rep);

        for u in g_orbit:
            G_members.remove(u);
            
        #for c in C.get_character(x):
            ##gc = group.mult(g, c);
            ##if (gc, xc) not in xg_orbit:
                ##xg_orbit.append((gc, xc));

        #if len(g_orbit) == len(xg_orbit):
            #"""
            #This happens precisely when there is v1, v2 in subgroup
            #such that for this u, 
                #u*v1 == u*v2 
            #and
                #chi(v1) != chi(v2).

            #In such a case, there will be fewer elements in orb_u
            #than in orb_xu.

            #This implies that f(u) = 0 (see our paper draft), and
            #this leads to badness ... ? So we do not count it.

            #(THIS IS SORT OF CRUDE -- CONSIDER REFIT)
            #"""
            #orbit_reps.append([g_orbit[0], len(g_orbit)])
        #else:
            #print("Orbit length not equal!");

        #G_elements = [u for u in G_elements if u not in g_orbit];

    return orbit_reps;


def _get_orbits_chi(group, subgroup, chi):
    """
    _get_orbit_representatives_IndVx()
    ----------------------------------
    Compute the orbits of the action of C on G and return a
    list of representatives for the orbits, i.e. a list of
    elements b_i in G, for each element [b_i] in the orbit space G/C. 

    @G: The monoid G
    @C: The subgroup C<G
    @x: The index of a primitive/principal character 

    Return: List of orbit representatives.

    """
    orbits = [];
    G      = group.members;
    while len(G) != 0:
        u      = G[0];
        orb_u  = [];
        orb_xu = [];

        for v in subgroup:
            uv = group.mult(u, v);
            if uv not in orb_u:
                orb_u.append(uv)

        for (v, xv) in subgroup.get_character(chi):
            uv = group.mult(u, v);
            if (uv, xv) not in orb_xu:
                orb_xu.append((uv, xv));

        if len(orb_u) == len(orb_xu):
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
            orbits.append([orb_u[0], len(orb_u)])
        else:
            print("Orbit length not equal!");

        G = [g for g in G if g not in orb_u];

    return orbits;


def _get_orbit_representatives_IndVx(G, C, x):
    """
    _get_orbit_representatives_IndVx()
    ----------------------------------
    Compute the orbits of the action of C on G by left-multiplication 

    @G: The monoid G
    @C: The subgroup C<G
    @x: The index of a primitive/principal character 

    Return: List of b_i in G, for each [b_i] in the orbit space G/C. 

    FOR THEORY: Appendix B.
    """
    orbit_reps = [];
    G_members  = list(G.members); # Copying this isn't my fav. thing ever. 

    while len(G_members) != 0:
        g = G_members[0];
        g_orbit = {};
        g_orbit_rep = None;
        
        for c in C:
            gc = G.mult(g, c);

            if g_orbit_rep == None:
                g_orbit_rep = gc;

            if gc not in g_orbit:
                g_orbit[gc] = c; 
            else: 
                d = g_orbit[gc]; # gc == gd
                if C.get_character_value(x, c) != C.get_character_value(x, d):
                    g_orbit_rep = None;
                    break;

        if g_orbit_rep != None:
            orbit_reps.append(g_orbit_rep);

        for u in g_orbit:
            G_members.remove(u);

    return orbit_reps;


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


def ActionA(a, C, chi):
    M = []

    for o1 in orbits_chi(G):
        V = []
        for o2 in orbits_chi(G):
            char = C.get_character(chi); # Where chi is chosen beforehand.

            H = [C.mult(x[0], o2) for x in char];

            # TODO: This scalar multiplication should be defined on the class.
            if (a*o1[0], a*o1[1]) in H:
                g = H.index((a*o1[0], a*o1[1]))                

                V.append(legendre(a)**C.monoid.param.k*echar(Z(char[g][1])));
            else:
                V.append(0)

        M.append(V)    

    return matrix(M).transpose()



class TanakaCharacter(CyclicCharacter):
    """
    Implements a character on a TanakaCyclicAbelianSubgroup.

    NOTES
    A TanakaCyclicAbelianSubgroup C is always finite, so 
    the character is stored as a list of tuples

        (c, xc) in C x Integers,

    where c is the element of C, and if C=<g> and c=g^j, 
    then xc=j.
    
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




#print(R.G.norm((0,15)));
#print(R.G.norm((0,12)));
#print(R.G.norm((0,26)));
#print(R.G.norm((0,1)));

#for g in R.G:
    #if R.G.norm(g) == 9:
        #print(g);



#exit();

#print(R.G.norm((0,6)));
#print(R.G.norm((0,8)));
#print(R.G.norm((0,3)));
#print(R.G.norm((0,1)));
#print(R.G.mult((0,3),(0,8)));
#print(R.C.mult((0,3),(0,8)));

#exit();

#def _mult_(u, v, p, n, k, D):
    #"""
    #Implements the special multiplication defined on the monoid.
    #@u    : Integer tuple (u0,u1) in Z/(p^n)Z x Z/(p^(n-k))Z
    #@v    : Integer tuple (v0,v1) in Z/(p^n)Z x Z/(p^(n-k))Z
    #Return: Integer tuple (w0,w1) in Z/(p^n)Z x Z/(p^(n-k))Z
    #"""
    #w0 = u[0]*v[0] - (D*(p**k))*u[1]*v[1];
    #w1 = u[0]*v[1] + u[1]*v[0];

    #print(w0);
    #print(w1);

    #print(mod(w0, p**n));
    #print(mod(w1, p**(n-k)));

    #res = (Integer(mod(w0, p**n)), Integer(mod(w1, p**(n-k))));

    #print(res);

#_mult_((0,6),(0,1), 3, 3, 0, 2);
#print('---');
#_mult_((0,3),(0,8), 3, 3, 0, 2);




exit();
#R.set_primitive_character(X[0]);


#for chi in R.get_primitive_characters():
    #R.set_primitive_character(chi);

    ##print(R.B(3));
    ##print(R.A(4));
    #print(R.W());
    #print("----");



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
    #orbit_reps = [];
    #G_members  = list(G.members); 
    #X          = dict(C.powers());

    #print(C.powers());

    #while len(G_members) != 0:
        #g          = G_members[0];
        #orbit_of_g = {};
        #for c in C:
            #gc = C.mult(g, c);
            #if gc not in orbit_of_g:
                #orbit_of_g[gc] = c; 
            #else: 
                #d = orbit_of_g[gc]; # gc == gd
                #if X[c] != X[d]:
                    #break;
        #else:
            #orbit_reps.append(C.mult(g, C[0]));

        #for u in orbit_of_g:
            #G_members.remove(u);

    #return orbit_reps;
