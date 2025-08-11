import sys
import time
from sympy import symbols, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

# Local imports
from . import geometry, stats
from table_generation import schemes, methods

if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)


# Evaluate the table for a given scheme with concrete numbers to output +-
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


def evaluate_iteration(iteration):
    # print(f"---------------------------------------------------------- At iteration {iteration}")
    start = time.time()

    pl, pi, pv, pj, pu, pk = geometry.generateSegments((1, 1000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pl1, pl2, pi1, pi2, pv1, pv2,
                                        pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pl1, pl2, pi1, pi2, pv1, pv2,
                                        pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)
    signSoS, depthSoS = evaluateTable(pExpressionsSoS, pl1, pl2, pi1, pi2, pv1, pv2,
                                      pj1, pj2, pu1, pu2, pk1, pk2, pl, pi, pv, pj, pu, pk)

    end = time.time()
    # print(f"Time for expression sign evaluation : {end - start:.6f} seconds")

    return {
        "signYapL": signYapL,
        "depthYapL": depthYapL,
        "signYapT": signYapT,
        "depthYapT": depthYapT,
        "signSoS": signSoS,
        "depthSoS": depthSoS
    }


# Set up symbols for the evaluation tables
pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

# Compute the evaluation tables for each scheme
print("Computing tables for Yap Lex...")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTableSegmentOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, "lex")

print("Computing tables for Yap Total...")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTableSegmentOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, "total")

print("Computing tables for SoS...")
pExpressionsSoS, eExpressionsSos = schemes.getEvaluationTableSegmentOrderSoS(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [methods.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [methods.count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [methods.count_ops(p) for p in pExpressionsSoS]


# for index in range(len(pExpressionsYapTotal)):
    # print(f"Index: {index}")
    # print(f"expression: {pExpressionsYapTotal[index]}")
    # print(f"operations: {operationCountYapTotal[index]}")

# print("\n\n")

# for index in range(len(pExpressionsSoS)):
    # print(f"Index: {index}")
    # print(f"expression: {pExpressionsSoS[index]}")
    # print(f"operations: {operationCountSoS[index]}")

# exit()


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



# Run all tests in paralle
with ProcessPoolExecutor() as executor:
    results = list(executor.map(evaluate_iteration, range(n)))

for r in results:
    signsYapL.append(r["signYapL"])
    depthsYapL.append(r["depthYapL"])
    operationsYapL.append(sum(operationCountYapLex[:r["depthYapL"]]))
    depthsHistogramYapL[r["depthYapL"]] += 1

    signsYapT.append(r["signYapT"])
    depthsYapT.append(r["depthYapT"])
    operationsYapT.append(sum(operationCountYapTotal[:r["depthYapT"]]))
    depthsHistogramYapT[r["depthYapT"]] += 1

    signsSoS.append(r["signSoS"])
    depthsSoS.append(r["depthSoS"])
    operationsSoS.append(sum(operationCountSoS[:r["depthSoS"]]))
    depthsHistogramSoS[r["depthSoS"]] += 1


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
stats.printStats(depthsSoS)
print("\nHere are the operations stats for SoS")
stats.printStats(operationsSoS)
print(f"\nHere's the histogram:")
for key, value in depthsHistogramSoS.items():
    print(f"Depth: {key}, count: {value}")
