import time

from sympy import symbols, IndexedBase, latex
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex

# Local import
from . import predicates, sos, yap, allienz



def getEvaluationTableSegmentOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, orderingType):

    expression = predicates.dualizeAndOrientYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)
    # expression = predicates.parametrizeAndOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)

    # The 3 here is hardcoded because it is enough for the segment order predicate

    if (orderingType == 'total'):
        pExpressions, eExpressions = yap.all_partials_orderedTotal(expression, [pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2], 4)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = yap.all_partials_orderedLex(expression, [pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2], 4)

    # Filter out the zero expressions
    pExpressionsNonZero = []
    eExpressionsNonZero = []

    for pExpression, eExpression in zip(pExpressions, eExpressions):
        if not pExpression.is_zero:
            pExpressionsNonZero.append(pExpression)
            eExpressionsNonZero.append(eExpression)

            if pExpression.is_constant():
                break

    return pExpressionsNonZero, eExpressionsNonZero 


def getEvaluationTableSegmentOrderSoS(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

    # Symbolic indexed bases
    p = IndexedBase('p')
    e = IndexedBase('e')

    # Symbolic variables for the indices
    variables = symbols("i j k l u v")
    i, j, k, l, u, v = variables

    # print("Computing tables for SoS")

    det = predicates.dualizeAndOrient(p, e, variables)
    pExpressionsSos, eExpressionsSos = sos.getEvaluationTable(p, e, variables, det)

    # Replace the indexed base with symbols which can be evaluated
    replacements = {
        p[l, 1]: pl1,
        p[l, 2]: pl2,
        p[i, 1]: pi1,
        p[i, 2]: pi2,
        p[v, 1]: pv1,
        p[v, 2]: pv2,
        p[j, 1]: pj1,
        p[j, 2]: pj2,
        p[u, 1]: pu1,
        p[u, 2]: pu2,
        p[k, 1]: pk1,
        p[k, 2]: pk2,
    }

    pExpressionsSoS_substituted = [expr.subs(replacements) for expr in pExpressionsSos]

    return pExpressionsSoS_substituted, eExpressionsSos



def getEvaluationTableSegmentOrderAllienz(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

    # Symbolic indexed bases
    p = IndexedBase('p')
    e = IndexedBase('e')

    # Symbolic variables for the indices
    variables = symbols("i j k l u v")
    i, j, k, l, u, v = variables

    # print("Computing tables for SoS")

    det = predicates.dualizeAndOrient(p, e, variables)
    pExpressionsSos, eExpressionsSos = sos.getEvaluationTable(p, e, variables, det)

    # Replace the indexed base with symbols which can be evaluated
    replacements = {
        p[l, 1]: pl1,
        p[l, 2]: pl2,
        p[i, 1]: pi1,
        p[i, 2]: pi2,
        p[v, 1]: pv1,
        p[v, 2]: pv2,
        p[j, 1]: pj1,
        p[j, 2]: pj2,
        p[u, 1]: pu1,
        p[u, 2]: pu2,
        p[k, 1]: pk1,
        p[k, 2]: pk2,
    }

    pExpressionsSoS_substituted = [expr.subs(replacements) for expr in pExpressionsSos]

    return pExpressionsSoS_substituted, eExpressionsSos




















def getEvaluationTablePointOrientationYap(pi1, pi2, pj1, pj2, pk1, pk2, orderingType):

    expression = predicates.orientationTestYap([pi1, pi2], [pj1, pj2], [pk1, pk2])

    # print(f"The expressions is {expression}")

    # The 3 here is hardcoded because it is enough for the segment order predicate

    if (orderingType == 'total'):
        pExpressions, eExpressions = yap.all_partials_orderedTotal(expression, [pi1, pi2, pj1, pj2, pk1, pk2], 3)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = yap.all_partials_orderedLex(expression, [pi1, pi2, pj1, pj2, pk1, pk2], 3)

    # Filter out the zero expressions
    pExpressionsNonZero = []
    eExpressionsNonZero = []

    for pExpression, eExpression in zip(pExpressions, eExpressions):

        if not pExpression.is_zero:
            pExpressionsNonZero.append(pExpression)
            eExpressionsNonZero.append(eExpression)

            if pExpression.is_constant():
                break

    return pExpressionsNonZero, eExpressionsNonZero 


def getEvaluationTablePointOrientationSoS(pi1, pi2, pj1, pj2, pk1, pk2):

    # Symbolic indexed bases
    p = IndexedBase('p')
    e = IndexedBase('e')

    # Symbolic variables for the indices
    variables = symbols("i j k l u v")
    i, j, k, l, u, v = variables

    # The expression
    det = predicates.orientationTest(p,
            (pi1 + e[i, 1], pi2 + e[i, 2]), 
            (pj1 + e[j, 1], pj2 + e[j, 2]), 
            (pk1 + e[k, 1], pk2 + e[k, 2])
            )

    return sos.getEvaluationTable(p, e, variables, det)




def getEvaluationTableSegmentOrderAllienz(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):
    e = symbols('e')
    pi = allienz.perturbPointAllienz([pi1, pi2], e)
    pj = allienz.perturbPointAllienz([pj1, pj2], e)
    pk = allienz.perturbPointAllienz([pk1, pk2], e)
    pl = allienz.perturbPointAllienz([pl1, pl2], e)
    pu = allienz.perturbPointAllienz([pu1, pu2], e)
    pv = allienz.perturbPointAllienz([pv1, pv2], e)

    p = IndexedBase('p')
    expression = expand(predicates.dualizeAndOrientAllienz(p, pl, pi, pv, pj, pu, pk))

    return allienz.computeEvaluationTable(expression, e)

def getEvaluationTablePointOrientationAllienz(pi1, pi2, pj1, pj2, pk1, pk2):

    e = symbols('e')
    pi = allienz.perturbPointAllienz([pi1, pi2], e)
    pj = allienz.perturbPointAllienz([pj1, pj2], e)
    pk = allienz.perturbPointAllienz([pk1, pk2], e)

    p = IndexedBase('p')
    expression = predicates.orientationTest(p, pi, pj, pk)
    return allienz.computeEvaluationTable(expression, e)
