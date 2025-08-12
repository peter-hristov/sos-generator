import sys
import time
from collections import defaultdict

from sympy import symbols, sign, count_ops

# Local imports
from . import geometry, stats
from table_generation import schemes 

# Evaluate the table for a given scheme with concrete numbers to output +-
def evaluateTable(pExpressions, pi, pj, pk, pi_c, pj_c, pk_c):

    subs_dict = {
            pi[0]: pi_c[0], pi[1]: pi_c[1],
            pj[0]: pj_c[0], pj[1]: pj_c[1],
            pk[0]: pk_c[0], pk[1]: pk_c[1],
            }

    for i, pExpression in enumerate(pExpressions):
        expressionSign = sign(pExpression.subs(subs_dict).evalf())

        if (expressionSign != 0):
            return expressionSign, i

    raise ValueError("Table could not be evaluated") 


# Handle input arguments
if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)

# Set up symbolic coordinates of the input points
pi = symbols("pi1, pi2")
pj = symbols("pj1, pj2")
pk = symbols("pk1, pk2")

# Compute the evaluation tables for each scheme
print("Computing tables for Yap Lex...")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTablePointOrientationYap(pi, pj, pk, "lex")

print("Computing tables for Yap Total...")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTablePointOrientationYap(pi, pj, pk, "total")

print("Computing tables for SoS...")
pExpressionsSoS, eExpressionsSoS = schemes.getEvaluationTablePointOrientationSoS(pi, pj, pk)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [count_ops(p) for p in pExpressionsSoS]


# Print for debugging purposes 
# 
for index in range(len(pExpressionsSoS)):
    print(f"Index: {index}")
    print(f"expression: {pExpressionsSoS[index]}")
    print(f"e-Term: {eExpressionsSoS[index]}")
    print(f"operations: {operationCountSoS[index]}")

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

signsSoS = []
depthsSoS = []
operationsSoS = []
depthsHistogramSoS = defaultdict(int)

for iteration in range(0, n):

    print(f"---------------------------------------------------------- At iteration {iteration}")
    start = time.time()

    # Generate concurrent points
    pi_c, pj_c, pk_c = geometry.generateColinearPoints((1, 1000000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pi, pj, pk, pi_c, pj_c, pk_c)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pi, pj, pk, pi_c, pj_c, pk_c)
    signSoS, depthSoS = evaluateTable(pExpressionsSoS, pi, pj, pk, pi_c, pj_c, pk_c)

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

print("\n\n------------------------------------------------------------------- SoS")
print("Here are the depth stats for SoS")
stats.printStats(depthsSoS)
print("\nHere are the operations stats for SoS")
stats.printStats(operationsSoS)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramSoS.items():
    print(f"Depth: {key}, count: {value}")
