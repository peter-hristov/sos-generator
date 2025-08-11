import sys
import time
from sympy import symbols, sign 
from collections import defaultdict

# Local imports
from . import geometry, stats
from table_generation import schemes, utility 

# Handle input arguments
if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)

# Evaluate the table for a given scheme with concrete numbers to output +-
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



# Set up symbols for the evaluation tables
pi1, pi2 = symbols("pi1, pi2")
pj1, pj2 = symbols("pj1, pj2")
pk1, pk2 = symbols("pk1, pk2")

# Compute the evaluation tables for each scheme
print("Computing tables for Yap Lex...")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTablePointOrientationYap(pi1, pi2, pj1, pj2, pk1, pk2, "lex")

print("Computing tables for Yap Total...")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTablePointOrientationYap(pi1, pi2, pj1, pj2, pk1, pk2, "total")

print("Computing tables for SoS...")
pExpressionsSos, eExpressionsSos = schemes.getEvaluationTablePointOrientationSoS(pi1, pi2, pj1, pj2, pk1, pk2)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [utility.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [utility.count_ops(p) for p in pExpressionsYapTotal]
operationCountSos = [utility.count_ops(p) for p in pExpressionsSos]


# Print for debugging purposes 
# 
for index in range(len(pExpressionsSos)):
    print(f"Index: {index}")
    print(f"expression: {pExpressionsSos[index]}")
    print(f"e-Term: {eExpressionsSos[index]}")
    print(f"operations: {operationCountSos[index]}")

# for index in range(len(pExpressionsYapLex)):
    # print(f"Index: {index}")
    # print(f"expression: {pExpressionsYapLex[index]}")
    # print(f"operations: {operationCountYapLex[index]}")

# print("\n\n")
# for index in range(len(pExpressionsYapTotal)):
    # print(f"Index: {index}")
    # print(f"expression: {pExpressionsYapTotal[index]}")
    # print(f"operations: {operationCountYapTotal[index]}")


# Set up arrays to hold the results for each test
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
    pi, pj, pk = geometry.generateColinearPoints((1, 1000000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)
    signSos, depthSos = evaluateTable(pExpressionsSos, pi1, pi2, pj1, pj2, pk1, pk2, pi, pj, pk)

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


# Output stats over all tests
print("\n\n------------------------------------------------------------------- Yap Lex")
print("Here are the depth stats for Yap Lex")
stats.printStats(depthsYapL)
print("\nHere are the operations stats for Yap Lex")
stats.printStats(operationsYapL)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramYapL.items():
    print(f"Depth: {key}, count: {value}")

print("\n\n------------------------------------------------------------------- Yap Total")
print("Here are the depth stats for Yap Total")
stats.printStats(depthsYapT)
print("\nHere are the operations stats for Yap Total")
stats.printStats(operationsYapT)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramYapT.items():
    print(f"Depth: {key}, count: {value}")

print("\n\n------------------------------------------------------------------- Sos")
print("Here are the depth stats for SoS")
stats.printStats(depthsSos)
print("\nHere are the operations stats for SoS")
stats.printStats(operationsSos)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramSos.items():
    print(f"Depth: {key}, count: {value}")
