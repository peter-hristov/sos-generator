# This code implements the symbolic perturbation scheme from the following paper Section 3.4
# 
# Alliez, P., Devillers, O. and Snoeyink, J., 2000. Removing degeneracies by perturbing the problem or perturbing the world. Reliable Computing, 6(1), pp.61-79.

from sympy import simplify, factor, collect 

# Local imports
from . import predicates 


def perturbPointAlliez(p, e):
    return [p[0] + e * p[1], p[1] + (e**2) * p[0] + (e**3) * (p[0]**2 + p[1]**2)]


def computeEvaluationTable(expression, e):
    
    pExpressions = []
    eExpressions = []

    # the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
    for index in range(0, 20):

        coefficient = e**index


        terms = sum([
            term for term in expression.as_ordered_terms()
            if term.as_coeff_exponent(e)[1] == index
        ])


        # There is no e, so there's nothing to factor
        if index == 0:
            pExpression = factor(simplify(terms))
        else:
            pExpression = factor(simplify(collect(terms, e**index).coeff(e**index)))

        # if (methods.count_ops(pExpression) == 0):
            # continue

        # print("-------------------------------------------------------------------")
        # print(f"Index: {index}.")
        # print(f"Expression: {pExpression}.")
        # print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
        # print(f"Number of operations: {utility.count_ops(pExpression)}")

        # print(f"Adding Expression: {pExpression}.")


        if not pExpression.is_zero:
            pExpressions.append(pExpression)
            eExpressions.append(e**index)

            if pExpression.is_constant():
                break


        # pExpressions.append(pExpression)
        # eExpressions.append(index)

    return pExpressions, eExpressions

