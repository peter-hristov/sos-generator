from sympy import symbols, count_ops, latex

# Local imports
from . import schemes 

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

print("Computing tables for Alliez...")
pExpressionsAlliez, eExpressionsAlliez = schemes.getEvaluationTablePointOrientationAlliez(pi, pj, pk)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [count_ops(p) for p in pExpressionsSoS]
operationCountAlliez = [count_ops(p) for p in pExpressionsAlliez]


print("-------------------------------------------------------")
print("Printing the table for Yap Lex...")
print("-------------------------------------------------------")
for index in range(len(pExpressionsYapLex)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsYapLex[index])}")
    print(f"terms       : {len(pExpressionsYapLex[index].as_ordered_terms())}")
    print(f"operations  : {operationCountYapLex[index]}")
    print("")

print("\n\n")
print("-------------------------------------------------------")
print("Printing the table for Yap Total...")
print("-------------------------------------------------------")
for index in range(len(pExpressionsYapTotal)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsYapTotal[index])}")
    print(f"terms       : {len(pExpressionsYapTotal[index].as_ordered_terms())}")
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
    print(f"terms       : {len(pExpressionsSoS[index].as_ordered_terms())}")
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
    print(f"terms       : {len(pExpressionsAlliez[index].as_ordered_terms())}")
    print(f"operations  : {operationCountAlliez[index]}")
    print("")
