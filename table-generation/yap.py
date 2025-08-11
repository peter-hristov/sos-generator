import sympy
from collections import defaultdict
import time
import sys
from sympy import symbols, diff, IndexedBase, latex, sign

import methods
import schemes


if len(sys.argv) != 3:
    print("Usage: python your_script.py <num_iterations> <ordering type \in {'lex', 'total'}")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)

orderingType = sys.argv[2]

pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

pExpressions, eExpressions = schemes.getEvaluationTableYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, orderingType)

print(pExpressions[0])
exit()

signs = []
depths = []
depthsHistogram = defaultdict(int)

for iteration in range(1, n):

    print(f"---------------------------------------------------------- At iteration {iteration}")

    pl, pi, pv, pj, pu, pk = methods.generateSegments((1, 1000))

    subs_dict = {
        pl1: pl[0], pl2: pl[1],
        pi1: pi[0], pi2: pi[1],
        pv1: pv[0], pv2: pv[1],
        pj1: pj[0], pj2: pj[1],
        pu1: pu[0], pu2: pu[1],
        pk1: pk[0], pk2: pk[1],
        }

    # print("The input points are:")
    # print(f"pl = {pl}")
    # print(f"pi = {pi}")
    # print(f"pv = {pv}")
    # print(f"pj = {pj}")
    # print(f"pu = {pu}")
    # print(f"pk = {pk}")

    expressionSign = 0

    # Evaluate on sign on our input
    start = time.time()
    depth = -1

    for index in range(len(pExpressions)):
        pExpression = pExpressions[index]
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

    end = time.time()

    # assert depth > 0

    signs.append(expressionSign)
    depths.append(depth)
    depthsHistogram[depth]+=1

    print(f"Time for expression sign evaluation : {end - start:.6f} seconds")
    print(f"The final sign is {expressionSign} at depth {depth}")

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
