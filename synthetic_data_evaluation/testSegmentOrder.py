import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from sympy import symbols, sign, count_ops

# Local imports
from . import geometry, stats
from table_generation import schemes 

# Evaluate the table for a given scheme with concrete numbers to output +-
def evaluateTable(pExpressions, pl, pi, pv, pj, pu, pk, pl_c, pi_c, pv_c, pj_c, pu_c, pk_c):

    subs_dict = {
            pl[0]: pl_c[0], pl[1]: pl_c[1],
            pi[0]: pi_c[0], pi[1]: pi_c[1],
            pv[0]: pv_c[0], pv[1]: pv_c[1],
            pj[0]: pj_c[0], pj[1]: pj_c[1],
            pu[0]: pu_c[0], pu[1]: pu_c[1],
            pk[0]: pk_c[0], pk[1]: pk_c[1],
            }

    for i, pExpression in enumerate(pExpressions):
        expressionSign = sign(pExpression.subs(subs_dict).evalf())

        if (expressionSign != 0):
            return expressionSign, i

    raise ValueError("Table could not be evaluated") 


def evaluate_iteration(iteration):
    # print(f"---------------------------------------------------------- At iteration {iteration}")
    start = time.time()

    pl_c, pi_c, pv_c, pj_c, pu_c, pk_c = geometry.generateSegments((1, 1000))

    signYapL, depthYapL = evaluateTable(pExpressionsYapLex, pl, pi, pv, pj, pu, pk, pl_c, pi_c, pv_c, pj_c, pu_c, pk_c)
    signYapT, depthYapT = evaluateTable(pExpressionsYapTotal, pl, pi, pv, pj, pu, pk, pl_c, pi_c, pv_c, pj_c, pu_c, pk_c)
    signSoS, depthSoS = evaluateTable(pExpressionsSoS, pl, pi, pv, pj, pu, pk, pl_c, pi_c, pv_c, pj_c, pu_c, pk_c)

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




if len(sys.argv) != 2:
    print("Usage: python your_script.py <num_iterations>")
    sys.exit(1)

try:
    n = int(sys.argv[1])
except ValueError:
    print("Invalid input: must be an integer")
    sys.exit(1)


# Set up symbols for the evaluation tables
pl = symbols("pl1, pl2")
pi = symbols("pi1, pi2")
pv = symbols("pv1, pv2")
pj = symbols("pj1, pj2")
pu = symbols("pu1, pu2")
pk = symbols("pk1, pk2")

# Compute the evaluation tables for each scheme
print("Computing tables for Yap Lex...")
pExpressionsYapLex, eExpressionsYapLex = schemes.getEvaluationTableSegmentOrderYap(pl, pi, pv, pj, pu, pk, "lex")

print("Computing tables for Yap Total...")
pExpressionsYapTotal, eExpressionsYapTotal = schemes.getEvaluationTableSegmentOrderYap(pl, pi, pv, pj, pu, pk, "total")

print("Computing tables for SoS...")
pExpressionsSoS, eExpressionsSoS = schemes.getEvaluationTableSegmentOrderSoS(pl, pi, pv, pj, pu, pk)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [count_ops(p) for p in pExpressionsSoS]

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
