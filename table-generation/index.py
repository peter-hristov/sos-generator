# Build the matrix
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex
from collections import defaultdict
import time
import sys
import matplotlib.pyplot as plt

import methods
import schemes

if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)

# Symbolic indexed bases
p = IndexedBase('p')
e = IndexedBase('e')

# Symbolic variables for the indices
variables = symbols("i j k l u v")
i, j, k, l, u, v = variables

# I'm not using this much any more, it used to be used to subs
# indexSubstitution = {i: 0, j: 1, k: 2, l: 3, u: 4, v: 5}

# Compute the evaluation table
pExpressions, eExpressions = schemes.getEvaluationTableSos(p, e, variables)

# Substitute with some symbols so that we can evalute this later efficiently
pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

# Replace the indexed base with symbols which can be evaluated
replacements = {
    p[l, 1]: pl1,
    p[l, 2]: pl2,
    p[i, 1]: pi1,
    p[i, 2]: pi2,
    p[v, 1]: pv1,
    p[v, 2]: pv2,
    p[j, 1]: pj1,
    p[j, 2]: pj2,
    p[u, 1]: pu1,
    p[u, 2]: pu2,
    p[k, 1]: pk1,
    p[k, 2]: pk2,
}

pExpressions_substituted = [expr.subs(replacements) for expr in pExpressions]


signs = []
depths = []
depthsHistogram = defaultdict(int)

for iteration in range(1, n):

    start = time.time()

    print(f"---------------------------------------------------------- At iteration {iteration}")
    # Blue
    # pl = generateRandomRationalCirclePoint()
    # pi = (-pl[0], -pl[1])

    # # Orange
    # pv = generateRandomRationalCirclePoint()
    # pj = (-pv[0], -pv[1])

    # # Green
    # pu = generateRandomRationalCirclePoint()
    # pk = (-pu[0], -pu[1])

    pl, pi, pv, pj, pu, pk = methods.generateSegments((1, 1000))

    subs_dict = {
        pl1: pl[0], pl2: pl[1],
        pi1: pi[0], pi2: pi[1],
        pv1: pv[0], pv2: pv[1],
        pj1: pj[0], pj2: pj[1],
        pu1: pu[0], pu2: pu[1],
        pk1: pk[0], pk2: pk[1],
        }


    expressionSign = 0

    # Evaluate on sign on our input
    depth = -1

    for index in range(len(pExpressions_substituted)):
        pExpression = pExpressions_substituted[index]
        eExpression = eExpressions[index]


        expressionSign = sign(pExpression.subs(subs_dict).evalf())
        depth += 1

        # print(f"-------------------------------------------- Derivative order {eExpression}:")
        # # print(f"{deriv}")
        # print(f"Depth {depth}")
        # print(f"The numeric value is {expressionSign}")
        # print("\n\n")

        if expressionSign != 0:
            break


    # print("The input points are:")
    # print(f"pl = {pl}")
    # print(f"pi = {pi}")
    # print(f"pv = {pv}")
    # print(f"pj = {pj}")
    # print(f"pu = {pu}")
    # print(f"pk = {pk}")

    # Evaluate on sign on our input
    # start = time.time()
    # expressionSign, depth = methods.evaluateExpresisonSign(pExpressions, eExpressions, p, variables, indexSubstitution, pl, pi, pv, pj, pu, pk)
    # end = time.time()

    # assert depth > 0

    signs.append(expressionSign)
    depths.append(depth)
    depthsHistogram[depth]+=1

    end = time.time()
    print(f"The final sign is {expressionSign} at depth {depth}")
    print(f"Time for expression sign evaluation : {end - start:.6f} seconds")

# After the loop
n = len(depths)
maxD = max(depths)
mean = sum(depths) / n
variance = sum((d - mean) ** 2 for d in depths) / (n - 1)  # sample stddev
stddev = variance ** 0.5

# plt.hist(depths, bins=100, edgecolor='black')
# plt.xlabel('Count')
# plt.ylabel('Frequency')
# plt.title('Histogram of Counts')
# plt.grid(True)
# plt.show()

print(f"Average depth: {mean:.3f}")
print(f"Max depth: {maxD:.3f}")
print(f"Standard deviation: {stddev:.3f}")

print(f"Here's the histogram:")
for key, value in depthsHistogram.items():
    print(f"Depth: {key}, count: {value}")
