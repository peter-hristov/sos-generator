# Build the matrix
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign
import time
import sys

import methods

# Symbolic indexed bases
p = IndexedBase('p')
e = IndexedBase('e')

# Symbolic variables
variables = symbols('i j k l u v')
i, j, k, l, u, v = variables

# The input to our method is the indices and coordinates of all points

# Indices of all points
indexSubstitution = {i: 0, j: 1, k: 2, l: 3, u: 4, v: 5}

# Coordinates of all points

# Blue
pl = (2, 9)
pi = (9, 2)

# Orange
pv = (2, 4)
pj = (10, 6)

# Green
pu = (4, 2)
pk = (9, 5)

print("Computing initial expression...")
# Affine

initialExpressionStart = time.time()

# det = methods.orientationTest(p,
        # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2]), 
        # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2]), 
        # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2])
        # )

# det = orientationTestHomogenious(
        # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2], p[i, 3] + e[i, 3]), 
        # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2], p[j, 3] + e[j, 3]), 
        # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2], p[k, 3] + e[k, 3])
        # )

# det = orientationTestHomogenious4D(
        # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2], p[i, 3] + e[i, 3], p[i, 4] + e[i, 4]), 
        # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2], p[j, 3] + e[j, 3], p[j, 4] + e[j, 4]), 
        # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2], p[k, 3] + e[k, 3], p[k, 4] + e[k, 4]),
        # (p[l, 1] + e[l, 1], p[l, 2] + e[l, 2], p[l, 3] + e[l, 3], p[l, 4] + e[l, 4])
        # )

det = methods.dualizeAndOrient(p, e, variables)
# det = methods.parametrizeAndOrder(p, e, variables)


# Expand, simplify and subsisute
expression = simplify(expand(det)).subs(indexSubstitution)
expressionTermsOrdered = expression.as_ordered_terms()
print(f"Total terms {len(expressionTermsOrdered)}")

initialExpressionEnd = time.time()


print("Arranging terms into different levels...")
levelsStart = time.time()
allTerms = []
# the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
for index in range(0, 100000):
    terms = sum([
        t for t in expressionTermsOrdered
        if methods.count_indexed_with_base(t, e) == index
    ])

    print(f"Index: {index}, terms: {terms}")

    if (terms == 0):
        break

    allTerms.append(terms)

levelsEnd = time.time()

print("Here are all the terms in the expression:")
for t in expressionTermsOrdered:
    print(t)

print("\n\nHere are all the terms grouped by the number of epsilons:")
for t in range(0, len(allTerms)):
    print(t)
    print(allTerms[t])

print(f"We have this many types of mixed expressions {len(allTerms)}")

sosStart = time.time()
pExpressions, eExpressions = methods.printSosTable(allTerms, p, e, indexSubstitution.values())
sosEnd = time.time()

start = time.time()
expressionSign = methods.evaluateExpresisonSign(pExpressions, eExpressions, allTerms, p, e, variables, indexSubstitution, pl, pi, pv, pj, pu, pk)
end = time.time()

print(f"Time for initial expression         : {initialExpressionEnd - initialExpressionStart:.6f} seconds")
print(f"Time for levels                     : {levelsEnd - levelsStart:.6f} seconds")
print(f"Time for sos table generation       : {sosEnd - sosStart:.6f} seconds")
print(f"Time for expression sign evaluation : {end - start:.6f} seconds")
print(f"\nThe final sign is {expressionSign}")


