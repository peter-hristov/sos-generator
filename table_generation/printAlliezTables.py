from sympy import (
    IndexedBase,
    symbols,
    simplify,
    factor,
    collect,
    latex
)

from . import predicates, allienz, utility, schemes





# Set up symbols for the evaluation tables
pl1, pl2 = symbols("pl1, pl2")
pi1, pi2 = symbols("pi1, pi2")
pv1, pv2 = symbols("pv1, pv2")
pj1, pj2 = symbols("pj1, pj2")
pu1, pu2 = symbols("pu1, pu2")
pk1, pk2 = symbols("pk1, pk2")

print("Computing Allienz evaluation table...")
# pExpressions, eExpressions = schemes.getEvaluationTablePointOrientationAllienz(pi1, pi2, pj1, pj2, pk1, pk2)
pExpressions, eExpressions = schemes.getEvaluationTableSegmentOrderAllienz(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)
print("-------------------------------------------------------------------")

print(pExpressions)
print(eExpressions)

for index in range(len(pExpressions)):
    pExpression = pExpressions[index] 
    eExpression = eExpressions[index] 
    # print(f"Index: {index}.")
    # # print(f"Expression: {pExpression}.")
    # print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
    # print(f"Number of operations: {methods.count_ops(pExpression)}")

    # print(f"{index}, {len(pExpression.as_ordered_terms())}, {utility.count_ops(pExpression)},")
    # print(f"Index: {index}.")
    # # print(f"Expression: {pExpression}.")
    # print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
    # print(f"Number of operations: {methods.count_ops(pExpression)}")
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressions[index])}")
    print(f"e-Term      : {latex(eExpressions[index])}")
    print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
    print(f"operations  : {utility.count_ops(pExpression)}")
    print("")














