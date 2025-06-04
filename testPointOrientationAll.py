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

def evaluateTable(pExpressions, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk):

    subs_dict = {
            pi1: pi[0], pi2: pi[1],
            pj1: pj[0], pj2: pj[1],
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
pi1, pi2 = symbols("pi1, pi2")
pj1, pj2 = symbols("pj1, pj2")
pk1, pk2 = symbols("pk1, pk2")

print("Computing tables for Yap Lex")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTableYapOrient(pi1, pi2, pj1, pj2, pk1, pk2, "lex")

for index in range(len(pExpressionsYapLex)):
    print(f"Index: {index}")
    print(f"expression: {pExpressionsYapLex[index]}")


print("Computing tables for Yap Total")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTableYapOrient(pi1, pi2, pj1, pj2, pk1, pk2, "total")

for index in range(len(pExpressionsYapTotal)):
    print(f"Index: {index}")
    print(f"expression: {pExpressionsYapTotal[index]}")



#
# Get SoS
#
# Symbolic indexed bases
p = IndexedBase('p')
e = IndexedBase('e')

# Symbolic variables for the indices
variables = symbols("i j k l u v")
i, j, k, l, u, v = variables

# The expression
det = methods.orientationTest(p,
        (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2]), 
        (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2]), 
        (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2])
        )

print("Computing tables for SoS")
pExpressionsSos, eExpressionsSos = schemes.getEvaluationTableSos(p, e, variables, det)

# Replace the indexed base with symbols which can be evaluated
replacements = {
    p[i, 1]: pi1,
    p[i, 2]: pi2,
    p[j, 1]: pj1,
    p[j, 2]: pj2,
    p[k, 1]: pk1,
    p[k, 2]: pk2,
}

pExpressionsSos_substituted = [expr.subs(replacements) for expr in pExpressionsSos]



# Count the number of operations in the tables for stats
operationCountYapLex = [methods.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [methods.count_ops(p) for p in pExpressionsYapTotal]
operationCountSos = [methods.count_ops(p) for p in pExpressionsSos_substituted]

# print("\n\n--------------------------------------------------------------------------- Yap Lex operations")
# for index in range(len(operationCountYapLex)):
    # print("--------------------------")
    # print(f"Index: {index}")
    # print(f"Expression: {pExpressionsYapLex[index]}")
    # print(f"#Terms: {len(pExpressionsYapLex[index].as_ordered_terms())}")
    # print(f"#Operations: {operationCountYapLex[index]}")

# print("\n\n--------------------------------------------------------------------------- Yap Total operations")
# for index in range(len(operationCountYapTotal)):
    # print("--------------------------")
    # print(f"Index: {index}")
    # print(f"Expression: {pExpressionsYapTotal[index]}")
    # print(f"#Terms: {len(pExpressionsYapTotal[index].as_ordered_terms())}")
    # print(f"#Operations: {operationCountYapTotal[index]}")

# print("\n\n--------------------------------------------------------------------------- Sos operations")
# for index in range(len(operationCountSos)):
    # print("--------------------------")
    # print(f"Index: {index}")
    # print(f"Expression: {pExpressionsSos[index]}")
    # print(f"#Terms: {len(pExpressionsSos[index].as_ordered_terms())}")
    # print(f"#Operations: {operationCountSos[index]}")

# exit()

signsYapL = []
depthsYapL = []
operationsYapL = []
depthsHistogramYapL = defaultdict(int)

signsYapT = []
depthsYapT = []
operationsYapT = []
depthsHistogramYapT = defaultdict(int)

signsSos = []
depthsSos = []
operationsSos = []
depthsHistogramSos = defaultdict(int)

for iteration in range(0, n):

    print(f"---------------------------------------------------------- At iteration {iteration}")
    start = time.time()

    # Generate concurrent points
    pi, pj, pk = methods.generateColinearPoints((1, 1000000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)
    signSos, depthSos = evaluateTable(pExpressionsSos_substituted, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)

    signsYapL.append(signYapL)
    depthsYapL.append(depthYapL)
    operationsYapL.append(sum(operationCountYapLex[:depthYapL]))
    depthsHistogramYapL[depthYapL]+=1

    signsYapT.append(signYapT)
    depthsYapT.append(depthYapT)
    operationsYapT.append(sum(operationCountYapTotal[:depthYapT]))
    depthsHistogramYapT[depthYapT]+=1

    signsSos.append(signSos)
    depthsSos.append(depthSos)
    operationsSos.append(sum(operationCountSos[:depthSos]))
    depthsHistogramSos[depthSos]+=1

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
printStats(depthsSos)
print("\nHere are the operations stats for SoS")
printStats(operationsSos)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramSos.items():
    print(f"Depth: {key}, count: {value}")
