from sympy import symbols, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex

from . import schemes, utility

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

print("Computing tables for Alliez...")
pExpressionsAlliez, eExpressionsAlliez = schemes.getEvaluationTableSegmentOrderAlliez(pl, pi, pv, pj, pu, pk)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [utility.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [utility.count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [utility.count_ops(p) for p in pExpressionsSoS]
operationCountAlliez = [utility.count_ops(p) for p in pExpressionsAlliez]

print("-------------------------------------------------------")
print("Printing the table for Yap Lex...")
print("-------------------------------------------------------")
for index in range(len(pExpressionsYapLex)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsYapLex[index])}")
    print(f"operations  : {operationCountYapLex[index]}")
    print("")

print("\n\n")
print("-------------------------------------------------------")
print("Printing the table for Yap Total...")
print("-------------------------------------------------------")
for index in range(len(pExpressionsYapTotal)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsYapTotal[index])}")
    print(f"operations  : {operationCountYapTotal[index]}")
    print("")

print("\n\n")
print("-------------------------------------------------------")
print("Printing the table for SoS")
print("-------------------------------------------------------")
for index in range(len(pExpressionsSoS)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsSoS[index])}")
    print(f"e-Term      : {latex(eExpressionsSoS[index])}")
    print(f"operations  : {operationCountSoS[index]}")
    print("")

print("\n\n")
print("-------------------------------------------------------")
print("Printing the table for Alliez")
print("-------------------------------------------------------")
for index in range(len(pExpressionsAlliez)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsAlliez[index])}")
    print(f"e-Term      : {latex(eExpressionsAlliez[index])}")
    print(f"operations  : {operationCountAlliez[index]}")
    print("")
