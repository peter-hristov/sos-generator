import itertools

import sympy
from collections import defaultdict
import time
import sys
from sympy import symbols, diff, IndexedBase, latex, sign

import methods
import schemes

variables = symbols("i j k l u v")

pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

p = [pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2]
e = IndexedBase('e')

schemes.getEvaluationTableSosNew(p, e, variables)
