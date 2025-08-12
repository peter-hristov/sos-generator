import time

from sympy import symbols, IndexedBase, latex
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex

# Local import
from . import predicates, sos, yap, alliez

def getEvaluationTablePointOrientationYap(pi1, pi2, pj1, pj2, pk1, pk2, orderingType):

    expression = predicates.orientationTest([pi1, pi2], [pj1, pj2], [pk1, pk2])

    # The depth 3 here is hardcoded because it is enough for the segment order predicate
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

def getEvaluationTableSegmentOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, orderingType):

    expression = predicates.dualizeAndOrient(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)
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




def getEvaluationTablePointOrientationSoS(pi1, pi2, pj1, pj2, pk1, pk2):

    # Set up the perturbed expression
    e = IndexedBase('e')
    i, j, k = symbols("i j k")

    # The expression
    perturbedExpression = predicates.orientationTest(
            (pi1 + e[i, 1], pi2 + e[i, 2]), 
            (pj1 + e[j, 1], pj2 + e[j, 2]), 
            (pk1 + e[k, 1], pk2 + e[k, 2])
            )

    return sos.getEvaluationTable(perturbedExpression, e, [i, j, k])

def getEvaluationTableSegmentOrderSoS(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

    # Set up the perturbed expression
    e = IndexedBase('e')
    i, j, k, l, u, v = symbols("i j k l u v")

    perturbedExpression = predicates.dualizeAndOrient(
            pl1 + e[l, 1], pl2 + e[l, 2], 
            pi1 + e[i, 1], pi2 + e[i, 2], 
            pv1 + e[v, 1], pv2 + e[v, 2], 
            pj1 + e[j, 1], pj2 + e[j, 2], 
            pu1 + e[u, 1], pu2 + e[u, 2], 
            pk1 + e[k, 1], pk2 + e[k, 2]
            )

    return sos.getEvaluationTable(perturbedExpression, e, [i, j, k, l, u, v])




def getEvaluationTablePointOrientationAlliez(pi1, pi2, pj1, pj2, pk1, pk2):

    # Set up the perturbed expression
    e = symbols('e')
    pi = alliez.perturbPointAlliez([pi1, pi2], e)
    pj = alliez.perturbPointAlliez([pj1, pj2], e)
    pk = alliez.perturbPointAlliez([pk1, pk2], e)

    perturbedExpression = predicates.orientationTest(pi, pj, pk)
    return alliez.computeEvaluationTable(perturbedExpression, e)


def getEvaluationTableSegmentOrderAlliez(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

    # Set up the perturbed expression
    e = symbols('e')
    pi = alliez.perturbPointAlliez([pi1, pi2], e)
    pj = alliez.perturbPointAlliez([pj1, pj2], e)
    pk = alliez.perturbPointAlliez([pk1, pk2], e)
    pl = alliez.perturbPointAlliez([pl1, pl2], e)
    pu = alliez.perturbPointAlliez([pu1, pu2], e)
    pv = alliez.perturbPointAlliez([pv1, pv2], e)

    perturbedExpression = expand(predicates.dualizeAndOrient(pl[0], pl[1], pi[0], pi[1], pv[0], pv[1], pj[0], pj[1], pu[0], pu[1], pk[0], pk[1]))
    return alliez.computeEvaluationTable(perturbedExpression, e)
