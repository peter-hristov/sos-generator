import sys
import time
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex

from collections import defaultdict
import schemes
import methods

if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)

def printStats(depths):
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

    print(f"Average: {float(mean):.3f}")
    print(f"Max: {float(maxD):.3f}")
    print(f"Standard deviation: {float(stddev):.3f}")

def evaluateTable(pExpressions, pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk):

    subs_dict = {
            pl1: pl[0], pl2: pl[1],
            pi1: pi[0], pi2: pi[1],
            pv1: pv[0], pv2: pv[1],
            pj1: pj[0], pj2: pj[1],
            pu1: pu[0], pu2: pu[1],
            pk1: pk[0], pk2: pk[1],
            }

    for i, pExpression in enumerate(pExpressions):
        expressionSign = sign(pExpression.subs(subs_dict).evalf())

        if (expressionSign != 0):
            return expressionSign, i

    raise ValueError("Table could not be evaluated") 


#
# Get Yap
#
pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

print("Computing tables for Yap Lex")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTableYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, "lex")

print("Computing tables for Yap Total")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTableYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, "total")



#
# Get SoS
#
# Symbolic indexed bases
p = IndexedBase('p')
e = IndexedBase('e')

# Symbolic variables for the indices
variables = symbols("i j k l u v")
i, j, k, l, u, v = variables

print("Computing tables for SoS")
pExpressionsSos, eExpressionsSos = schemes.getEvaluationTableSos(p, e, variables)

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

pExpressionsSoS_substituted = [expr.subs(replacements) for expr in pExpressionsSos]

# Count the number of operations in the tables for stats
operationCountYapLex = [methods.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [methods.count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [methods.count_ops(p) for p in pExpressionsSoS_substituted]

signsYapL = []
depthsYapL = []
operationsYapL = []
depthsHistogramYapL = defaultdict(int)

signsYapT = []
depthsYapT = []
operationsYapT = []
depthsHistogramYapT = defaultdict(int)

signsSoS = []
depthsSoS = []
operationsSoS = []
depthsHistogramSoS = defaultdict(int)

for iteration in range(0, n):

    print(f"---------------------------------------------------------- At iteration {iteration}")
    start = time.time()

    # Generate concurrent points
    pl, pi, pv, pj, pu, pk = methods.generateSegments((1, 1000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)
    signSoS, depthSoS = evaluateTable(pExpressionsSoS_substituted, pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)

    signsYapL.append(signYapL)
    depthsYapL.append(depthYapL)
    operationsYapL.append(sum(operationCountYapLex[:depthYapL]))
    depthsHistogramYapL[depthYapL]+=1

    signsYapT.append(signYapT)
    depthsYapT.append(depthYapT)
    operationsYapT.append(sum(operationCountYapTotal[:depthYapT]))
    depthsHistogramYapT[depthYapT]+=1

    signsSoS.append(signSoS)
    depthsSoS.append(depthSoS)
    operationsSoS.append(sum(operationCountSoS[:depthSoS]))
    depthsHistogramSoS[depthSoS]+=1

    end = time.time()

    print(f"Time for expression sign evaluation : {end - start:.6f} seconds")
    # print(f"The sign is {expressionSign} at depth {depth}")
    # print(f"The sign is {expressionSign} at depth {depth}")
    # print(f"The sign is {expressionSign} at depth {depth}")



print("\n\n------------------------------------------------------------------- Yap Lex")
print("Here are the depth stats for Yap Lex")
printStats(depthsYapL)
print("\nHere are the operations stats for Yap Lex")
printStats(operationsYapL)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramYapL.items():
    print(f"Depth: {key}, count: {value}")

print("\n\n------------------------------------------------------------------- Yap Total")
print("Here are the depth stats for Yap Total")
printStats(depthsYapT)
print("\nHere are the operations stats for Yap Total")
printStats(operationsYapT)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramYapT.items():
    print(f"Depth: {key}, count: {value}")

print("\n\n------------------------------------------------------------------- Sos")
print("Here are the depth stats for SoS")
printStats(depthsSoS)
print("\nHere are the operations stats for SoS")
printStats(operationsSoS)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramSoS.items():
    print(f"Depth: {key}, count: {value}")
