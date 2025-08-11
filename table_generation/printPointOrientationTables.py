from sympy import symbols, latex  

# Local imports
from . import schemes, utility 


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
pExpressionsSoS, eExpressionsSoS = schemes.getEvaluationTablePointOrientationSoS(pi1, pi2, pj1, pj2, pk1, pk2)

print("Computing tables for Allienz...")
pExpressionsAllienz, eExpressionsAllienz = schemes.getEvaluationTablePointOrientationAllienz(pi1, pi2, pj1, pj2, pk1, pk2)

# Compute number of arithemtic operations for each row of the evaluation table for each scheme
operationCountYapLex = [utility.count_ops(p) for p in pExpressionsYapLex]
operationCountYapTotal = [utility.count_ops(p) for p in pExpressionsYapTotal]
operationCountSoS = [utility.count_ops(p) for p in pExpressionsSoS]
operationCountAllienz = [utility.count_ops(p) for p in pExpressionsAllienz]


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
print("Printing the table for Allienz")
print("-------------------------------------------------------")
for index in range(len(pExpressionsAllienz)):
    print(f"Index       : {index}")
    print(f"expression  : {latex(pExpressionsAllienz[index])}")
    print(f"e-Term      : {latex(eExpressionsAllienz[index])}")
    print(f"operations  : {operationCountAllienz[index]}")
    print("")
