#!/usr/bin/sage -python

import sys
from sage.all import *

p = 3;
n = 3;
k = 1;

Zpn  = IntegerModRing(p**n);
Zpnk = IntegerModRing(p**(n-k));

Zpn1  = IntegerModRing(p**(n-1));
Zpn1k = IntegerModRing(p**(n-1-k));

p1 = p**(n-1);
p2 = p**(n-1-k);

M = [(a,b) for a in Zpn for b in Zpnk];
            
#A = [g for g in M if Zpn1(g[0])==Zpn1(1) and Zpn1k(g[1])==Zpn1k(0)];
#B = [g for g in M if mod(g[0], p1) == mod(1, p1) and mod(g[1], p2) == mod(0, p2)];
#C = [g for g in M if mod(g[0], p1) == 1 and mod(g[1], p2) == 0];
D = [g for g in M if mod(g[0], p**(n-1)) == 1 and mod(g[1], p**(n-1-k)) == 0];


E = [
    # Note the content of CL is dependent (via R1, R2) on k
    (g0,g1) for (g0,g1) in M if 
        (mod(g0, p**(n-1)), mod(g1, p**(n-1-k))) == (1,0)
],


print(E);
#print(A == B);
#print(B == C);
#print(C == D);
