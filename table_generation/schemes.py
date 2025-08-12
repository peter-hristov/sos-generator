from sympy import symbols, IndexedBase, expand

# Local imports
from . import predicates, sos, yap, alliez

def getEvaluationTablePointOrientationYap(pi, pj, pk, orderingType):

    expression = predicates.orientationTest(pi, pj, pk)

    # The depth 3 here is hardcoded because it is enough for the point orientation predicate
    if (orderingType == 'total'):
        pExpressions, eExpressions = yap.all_partials_orderedTotal(expression, [pi[0], pi[1], pj[0], pj[1], pk[0], pk[1]], 3)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = yap.all_partials_orderedLex(expression, [pi[0], pi[1], pj[0], pj[1], pk[0], pk[1]], 3)

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

def getEvaluationTableSegmentOrderYap(pl, pi, pv, pj, pu, pk, orderingType):

    expression = predicates.dualizeAndOrient(pl, pi, pv, pj, pu, pk)

    # The 4 here is hardcoded because it is enough for the segment order predicate
    if (orderingType == 'total'):
        pExpressions, eExpressions = yap.all_partials_orderedTotal(expression, [pl[0], pl[1], pi[0], pi[1], pv[0], pv[1], pj[0], pj[1], pu[0], pu[1], pk[0], pk[1]], 4)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = yap.all_partials_orderedLex(expression, [pl[0], pl[1], pi[0], pi[1], pv[0], pv[1], pj[0], pj[1], pu[0], pu[1], pk[0], pk[1]], 4)

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




def getEvaluationTablePointOrientationSoS(pi, pj, pk):

    # Set up the perturbed expression
    e = IndexedBase('e')
    i, j, k = symbols("i j k")

    pi = (pi[0] + e[i, 1], pi[1] + e[i, 2]) 
    pj = (pj[0] + e[j, 1], pj[1] + e[j, 2]) 
    pk = (pk[0] + e[k, 1], pk[1] + e[k, 2])

    perturbedExpression = predicates.orientationTest(pi, pj, pk)
    return sos.getEvaluationTable(perturbedExpression, e, [i, j, k])

def getEvaluationTableSegmentOrderSoS(pl, pi, pv, pj, pu, pk):

    # Set up the perturbed expression
    e = IndexedBase('e')
    i, j, k, l, u, v = symbols("i j k l u v")

    pl = (pl[0] + e[l, 1], pl[1] + e[l, 2]) 
    pi = (pi[0] + e[i, 1], pi[1] + e[i, 2]) 
    pv = (pv[0] + e[v, 1], pv[1] + e[v, 2]) 
    pj = (pj[0] + e[j, 1], pj[1] + e[j, 2]) 
    pu = (pu[0] + e[u, 1], pu[1] + e[u, 2]) 
    pk = (pk[0] + e[k, 1], pk[1] + e[k, 2])

    perturbedExpression = predicates.dualizeAndOrient(pl, pi, pv, pj, pu, pk)
    return sos.getEvaluationTable(perturbedExpression, e, [i, j, k, l, u, v])


def getEvaluationTablePointOrientationAlliez(pi, pj, pk):

    # Set up the perturbed expression
    e = symbols('e')
    pi = alliez.perturbPointAlliez(pi, e)
    pj = alliez.perturbPointAlliez(pj, e)
    pk = alliez.perturbPointAlliez(pk, e)

    perturbedExpression = predicates.orientationTest(pi, pj, pk)
    return alliez.computeEvaluationTable(perturbedExpression, e)


def getEvaluationTableSegmentOrderAlliez(pl, pi, pv, pj, pu, pk):

    # Set up the perturbed expression
    e = symbols('e')
    pi = alliez.perturbPointAlliez(pi, e)
    pj = alliez.perturbPointAlliez(pj, e)
    pk = alliez.perturbPointAlliez(pk, e)
    pl = alliez.perturbPointAlliez(pl, e)
    pu = alliez.perturbPointAlliez(pu, e)
    pv = alliez.perturbPointAlliez(pv, e)

    perturbedExpression = expand(predicates.dualizeAndOrient(pl, pi, pv, pj, pu, pk))
    return alliez.computeEvaluationTable(perturbedExpression, e)
