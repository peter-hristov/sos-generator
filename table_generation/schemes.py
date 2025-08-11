import time

from sympy import symbols, diff, IndexedBase, latex
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, latex

from . import methods

import itertools



def getEvaluationTableSosNew(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, e, variables):

    # Unpack variables
    i, j, k, l, u, v = variables

    # Perform an epislon expansion
    pl1e = pl1 + e[l, 1]
    pl2e = pl2 + e[l, 2]

    pi1e = pi1 + e[i, 1]
    pi2e = pi2 + e[i, 2]

    pv1e = pv1 + e[v, 1]
    pv2e = pv2 + e[v, 2]

    pj1e = pj1 + e[j, 1]
    pj2e = pj2 + e[j, 2]

    pu1e = pu1 + e[u, 1]
    pu2e = pu2 + e[u, 2]

    pk1e = pk1 + e[k, 1]
    pk2e = pk2 + e[k, 2]

    initialExpressionStart = time.time()

    expression = methods.dualizeAndOrientYap(pl1e, pl2e, pi1e, pi2e, pv1e, pv2e, pj1e, pj2e, pu1e, pu2e, pk1e, pk2e)
    expression = simplify(expand(expression))
    expressionTermsOrdered = expression.as_ordered_terms()

    initialExpressionEnd = time.time()
    print(f"Total terms {len(expressionTermsOrdered)}")


    # Break into terms and count e-products
    print("Arranging terms into different levels...")
    levelsStart = time.time()
    allTerms = []
    # the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
    for index in range(0, 100000):
        terms = sum([
            t for t in expressionTermsOrdered
            if methods.count_indexed_with_base(t, e) == index
        ])

        # print(f"Index: {index}, terms: {terms}")

        if (terms == 0):
            break

        allTerms.append(terms)

    levelsEnd = time.time()

    # print("Here are all the terms in the expression:")
    # for t in expressionTermsOrdered:
        # print(t)

    # print("\n\nHere are all the terms grouped by the number of epsilons:")
    # for t in range(0, len(allTerms)):
        # print(t)
        # print(allTerms[t])

    # print(f"We have this many types of mixed expressions {len(allTerms)}")

    sosStart = time.time()
    pExpressions, eExpressions = methods.printSosTable(allTerms, p, e, variables)
    sosEnd = time.time()

    print(f"Time for initial expression         : {initialExpressionEnd - initialExpressionStart:.6f} seconds")
    print(f"Time for levels                     : {levelsEnd - levelsStart:.6f} seconds")
    print(f"Time for sos table generation       : {sosEnd - sosStart:.6f} seconds")

    return pExpressions, eExpressions




def getEvaluationTableSos(p, e, variables, det):

    # Unpack variables
    i, j, k, l, u, v = variables

    print("Computing initial expression...")
    # Affine

    initialExpressionStart = time.time()

    # det = methods.orientationTest(p,
            # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2]), 
            # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2]), 
            # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2])
            # )

    # det = orientationTestHomogenious(
            # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2], p[i, 3] + e[i, 3]), 
            # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2], p[j, 3] + e[j, 3]), 
            # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2], p[k, 3] + e[k, 3])
            # )

    # det = orientationTestHomogenious4D(
            # (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2], p[i, 3] + e[i, 3], p[i, 4] + e[i, 4]), 
            # (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2], p[j, 3] + e[j, 3], p[j, 4] + e[j, 4]), 
            # (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2], p[k, 3] + e[k, 3], p[k, 4] + e[k, 4]),
            # (p[l, 1] + e[l, 1], p[l, 2] + e[l, 2], p[l, 3] + e[l, 3], p[l, 4] + e[l, 4])
            # )

    # det = methods.dualizeAndOrient(p, e, variables)
    # det = methods.parametrizeAndOrder(p, e, variables)


    # Expand, simplify and subsisute
    expression = simplify(expand(det))
    expressionTermsOrdered = expression.as_ordered_terms()
    print(f"Total terms {len(expressionTermsOrdered)}")

    initialExpressionEnd = time.time()


    print("Arranging terms into different levels...")
    levelsStart = time.time()
    allTerms = []
    # the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
    for index in range(0, 100000):
        terms = sum([
            t for t in expressionTermsOrdered
            if methods.count_indexed_with_base(t, e) == index
        ])

        # print(f"Index: {index}, terms: {terms}")

        if (terms == 0):
            break

        allTerms.append(terms)

    levelsEnd = time.time()

    # print("Here are all the terms in the expression:")
    # for t in expressionTermsOrdered:
        # print(t)

    # print("\n\nHere are all the terms grouped by the number of epsilons:")
    # for t in range(0, len(allTerms)):
        # print(t)
        # print(allTerms[t])

    # print(f"We have this many types of mixed expressions {len(allTerms)}")

    sosStart = time.time()
    pExpressions, eExpressions = methods.printSosTable(allTerms, p, e, variables)
    sosEnd = time.time()

    print(f"Time for initial expression         : {initialExpressionEnd - initialExpressionStart:.6f} seconds")
    print(f"Time for levels                     : {levelsEnd - levelsStart:.6f} seconds")
    print(f"Time for sos table generation       : {sosEnd - sosStart:.6f} seconds")

    return pExpressions, eExpressions



def all_partials_orderedTotal(f, vars, m):
    """
    Generate all partial derivatives of f with respect to vars
    up to total order m, ordered by total order and lex order:
    dx -> dy -> dz -> dx^2 -> dx dy -> dxdz ...

    Args:
      f: sympy expression
      vars: list of sympy symbols
      m: max total degree of derivatives

    Returns:
      List of tuples (orders, derivative_expr)
    """
    n = len(vars)
    pProducts = []
    eProducts = []

    def gen_multiindices_exact_sum(n, total_order, prefix=()):
        if len(prefix) == n:
            if sum(prefix) == total_order:
                yield prefix
            return
        for i in range(total_order + 1):
            if sum(prefix) + i <= total_order:
                yield from gen_multiindices_exact_sum(n, total_order, prefix + (i,))
            else:
                break

    for total_order in range(m + 1):
        for orders in gen_multiindices_exact_sum(n, total_order):
            deriv = f
            for var, order in zip(vars, orders):
                if order > 0:
                    deriv = diff(deriv, var, order)
            pProducts.append(deriv)
            eProducts.append(orders)

    return pProducts, eProducts



def all_partials_orderedLex(f, vars, m):
    """
    Generate all partial derivatives of f with respect to vars
    up to total order m, ordered by total order and lex order:
    dx -> dy -> dz -> dx^2 -> dx dy -> dxdz ...

    Args:
      f: sympy expression
      vars: list of sympy symbols
      m: max total degree of derivatives

    Returns:
      List of tuples (orders, derivative_expr)
    """
    n = len(vars)
    pProducts = []
    eProducts = []

    def gen_multiindices_lex(n, max_total_order):
        """
        Generate multi-indices (tuples of length n) with sum <= max_total_order,
        in lexicographic order: (0,0,1), (0,0,2), ..., (0,1,0), ...
        """
        for orders in itertools.product(range(max_total_order), repeat=n):
            if sum(orders) <= max_total_order:
                yield orders

    for orders in gen_multiindices_lex(n, m):
        # print(orders)
        deriv = f
        for var, order in zip(vars, orders):
            if order > 0:
                deriv = diff(deriv, var, order)
        pProducts.append(deriv)
        eProducts.append(orders)

    return pProducts, eProducts


def getEvaluationTableYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2, orderingType):

    expression = methods.dualizeAndOrientYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)
    # expression = methods.parametrizeAndOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2)

    # The 3 here is hardcoded because it is enough for the segment order predicate

    if (orderingType == 'total'):
        pExpressions, eExpressions = all_partials_orderedTotal(expression, [pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2], 4)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = all_partials_orderedLex(expression, [pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2], 4)

    # Filter out the zero expressions
    pExpressionsNonZero = []
    eExpressionsNonZero = []

    for pExpression, eExpression in zip(pExpressions, eExpressions):
        if not pExpression.is_zero:
            pExpressionsNonZero.append(pExpression)
            eExpressionsNonZero.append(eExpression)

    return pExpressionsNonZero, eExpressionsNonZero 





def getEvaluationTableYapOrient(pi1, pi2, pj1, pj2, pk1, pk2, orderingType):

    expression = methods.orientationTestYap([pi1, pi2], [pj1, pj2], [pk1, pk2])

    print(f"The expressions is {expression}")

    # The 3 here is hardcoded because it is enough for the segment order predicate

    if (orderingType == 'total'):
        pExpressions, eExpressions = all_partials_orderedTotal(expression, [pi1, pi2, pj1, pj2, pk1, pk2], 3)
    elif (orderingType == 'lex'):
        pExpressions, eExpressions = all_partials_orderedLex(expression, [pi1, pi2, pj1, pj2, pk1, pk2], 3)

    # Filter out the zero expressions
    pExpressionsNonZero = []
    eExpressionsNonZero = []

    for pExpression, eExpression in zip(pExpressions, eExpressions):
        if not pExpression.is_zero:
            pExpressionsNonZero.append(pExpression)
            eExpressionsNonZero.append(eExpression)

    return pExpressionsNonZero, eExpressionsNonZero 


def getEvaluationTableSoSOrient(pi1, pi2, pj1, pj2, pk1, pk2):

    # Symbolic indexed bases
    p = IndexedBase('p')
    e = IndexedBase('e')

    # Symbolic variables for the indices
    variables = symbols("i j k l u v")
    i, j, k, l, u, v = variables

    # The expression
    det = methods.orientationTest(p,
            (pi1 + e[i, 1], pi2 + e[i, 2]), 
            (pj1 + e[j, 1], pj2 + e[j, 2]), 
            (pk1 + e[k, 1], pk2 + e[k, 2])
            )

    return getEvaluationTableSos(p, e, variables, det)
