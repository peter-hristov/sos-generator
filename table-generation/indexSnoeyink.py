from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign, Rational, factor

import methods



def perturbPoint(p, e):
    return [p[0] + e * p[1], p[1] + (e**2) * p[0] + (e**3) * (p[0]**2 + p[1]**2)]
    # return [p[0] + e*p[1], p[1] + p[0]*e**2]
    # return [p[0], p[1]]

def dualizeAndOrientAllienz(p_l, p_i, p_v, p_j, p_u, p_k):

    # p_li homogenious
    p_li = (
        p_l[1] - p_i[1],
        p_i[0] - p_l[0],
        p_l[0] * p_i[1] - p_l[1] * p_i[0]
    )

    # p_vj homogenious
    p_vj = (
        p_v[1] - p_j[1],
        p_j[0] - p_v[0],
        p_v[0] * p_j[1] - p_v[1] * p_j[0]
    )

    # p_uk
    p_uk = (
        p_u[1] - p_k[1],
        p_k[0] - p_u[0],
        p_u[0] * p_k[1] - p_u[1] * p_k[0]
    )

    return methods.orientationTestHomogenious(p, p_li, p_vj, p_uk)




def getEvaluationTableAllienz(p, e, variables):

    pi = perturbPoint([p[i, 1], p[i, 2]], e)
    pj = perturbPoint([p[j, 1], p[j, 2]], e)

    pk = perturbPoint([p[k, 1], p[k, 2]], e)
    pl = perturbPoint([p[l, 1], p[l, 2]], e)

    pu = perturbPoint([p[u, 1], p[u, 2]], e)
    pv = perturbPoint([p[v, 1], p[v, 2]], e)

    det = methods.orientationTest(p, pi, pj, pk)

    print("The expression is:")
    print(det)

    # det = expand(dualizeAndOrientAllienz(pl, pi, pv, pj, pu, pk))

    # print(det)
    # num_terms = len(det.as_ordered_terms())
    # print("Number of terms in expanded det:", num_terms)

    pExpressions = []
    eExpressions = []

    # the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
    for index in range(0, 20):

        coefficient = e**index


        terms = sum([
            term for term in det.as_ordered_terms()
            if term.as_coeff_exponent(e)[1] == index
        ])


        # There is no e, so there's nothing to factor
        if index == 0:
            pExpression = factor(simplify(terms))
        else:
            pExpression = factor(simplify(collect(terms, e**index).coeff(e**index)))
            

        # if (methods.count_ops(pExpression) == 0):
            # continue

        print("-------------------------------------------------------------------")
        print(f"Index: {index}.")
        print(f"Expression: {pExpression}.")
        print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
        print(f"Number of operations: {methods.count_ops(pExpression)}")

        # print(f"Adding Expression: {pExpression}.")
        pExpressions.append(pExpression)
        eExpressions.append(index)

    return pExpressions, eExpressions



# Symbolic indexed bases
p = IndexedBase('p')
e = symbols('e')

variables = symbols("i j k l u v")
i, j, k, l, u, v = variables

print("Computing Allienz evaluation table...")
pExpressions, eExpressions = getEvaluationTableAllienz(p, e, variables)
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

    print(f"{index}, {len(pExpression.as_ordered_terms())}, {methods.count_ops(pExpression)},")
    # print(f"Index: {index}.")
    # # print(f"Expression: {pExpression}.")
    # print(f"Number of terms: {len(pExpression.as_ordered_terms())}")
    # print(f"Number of operations: {methods.count_ops(pExpression)}")














