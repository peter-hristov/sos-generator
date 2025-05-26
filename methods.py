from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode, sign

def evaluateExpression(expression, p, variables, indexSubstitution, p1, p2, p3, p4, p5, p6):
    i, j, k, l, u, v = variables

    # expression_index_subs = expression.subs(indexSubstitution)

    valueSubs = {
            p[indexSubstitution[i], 1]: p1[0], p[indexSubstitution[i], 2]: p1[1],
            p[indexSubstitution[j], 1]: p2[0], p[indexSubstitution[j], 2]: p2[1],
            p[indexSubstitution[k], 1]: p3[0], p[indexSubstitution[k], 2]: p3[1],
            p[indexSubstitution[l], 1]: p4[0], p[indexSubstitution[l], 2]: p4[1],
            p[indexSubstitution[u], 1]: p5[0], p[indexSubstitution[u], 2]: p5[1],
            p[indexSubstitution[v], 1]: p6[0], p[indexSubstitution[v], 2]: p6[1],
            }

    expression_value_subs = expression.subs(valueSubs)

    return expression_value_subs


def isInArray(expression, expressionArray):
    for e in expressionArray:
        # print(f"comparing {expression} and {e}")
        if expression == e or -expression == e:
            return True
    return False

def orientationTest(p, p1, p2, p3):
    M = Matrix([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1],
    ])

    return M.det()

def orientationTestHomogenious(p, p1, p2, p3):
    M = Matrix([
        [p1[0], p1[1], p1[2]],
        [p2[0], p2[1], p2[2]],
        [p3[0], p3[1], p3[2]],
    ])

    return M.det()

def orientationTestHomogenious4D(p, p1, p2, p3, p4):
    M = Matrix([
        [p1[0], p1[1], p1[2], p1[3]],
        [p2[0], p2[1], p2[2], p2[3]],
        [p3[0], p3[1], p3[2], p3[3]],
        [p4[0], p4[1], p4[2], p4[3]],
    ])

    return M.det()

def dualizeAndOrient(p, e, variables):
    # Homogenious

    i, j, k, l, u, v = variables

    # p_li
    p_l = (p[l, 1] + e[l, 1], p[l, 2] + e[l, 2])
    p_i = (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2])

    # p_vj
    p_v = (p[v, 1] + e[v, 1], p[v, 2] + e[v, 2])
    p_j = (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2])

    # p_uk
    p_u = (p[u, 1] + e[u, 1], p[u, 2] + e[u, 2])
    p_k = (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2])


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

    return orientationTestHomogenious(p, p_li, p_vj, p_uk)


def parametrizeAndOrder(p, e, variables):

    i, j, k, l, u, v = variables

    # p_l -> p_i
    r = Matrix([p[i,1] - p[l,1] + e[i,1] - e[l,1], p[i,2] - p[l,2] + e[i,2] - e[l,2]])

    # p_v -> p_j
    s1 = Matrix([p[j,1] - p[v,1] + e[j,1] - e[v,1], p[j,2] - p[v,2] + e[j,2] - e[v,2]])
    # p_u -> p_k
    s2 = Matrix([p[k,1] - p[u,1] + e[k,1] - e[u,1], p[k,2] - p[u,2] + e[k,2] - e[u,2]])

    # p_v -> p_l
    q1 = Matrix([p[l,1] - p[v,1] + e[l,1] - e[v,1], p[l,2] - p[v,2] + e[l,2] - e[v,2]])
    # p_u -> p_l
    q2 = Matrix([p[l,1] - p[u,1] + e[l,1] - e[u,1], p[l,2] - p[u,2] + e[l,2] - e[u,2]])


    q1xs1 = Matrix([
        q1.T.tolist()[0],  # convert row matrix to flat list
        s1.T.tolist()[0]
    ]).det()

    q2xs2 = Matrix([
        q2.T.tolist()[0],  # convert row matrix to flat list
        s2.T.tolist()[0]
    ]).det()

    q1xr = Matrix([
        q1.T.tolist()[0],  # convert row matrix to flat list
        r.T.tolist()[0]
    ]).det()

    q2xr = Matrix([
        q2.T.tolist()[0],  # convert row matrix to flat list
        r.T.tolist()[0]
    ]).det()

    rxs1 = Matrix([
        r.T.tolist()[0],  # convert row matrix to flat list
        s1.T.tolist()[0]
    ]).det()

    rxs2 = Matrix([
        r.T.tolist()[0],  # convert row matrix to flat list
        s2.T.tolist()[0]
    ]).det()


    return expand(q1xs1 * rxs2 - q2xs2 * rxs1)








def count_indexed_with_base(term, base):
    """Return the number of Indexed objects in a term with the given base."""
    return sum(1 for atom in term.atoms(Indexed) if atom.base == base)

def generate_sequence(e, variableIndices, dimensions=(1, 2)):

    # Sort the indies of all variables
    variableIndices = sorted(variableIndices)

    sequence = []
    # sequence.append([f"e[{variableIndices[0]}, {dimensions[1]}]", f"e[{variableIndices[0]}, {dimensions[0]}]"])

    firstRow = []
    for j in reversed(range(len(dimensions))):
        firstRow.append(e[variableIndices[0], dimensions[j]])

    sequence.append(firstRow)
    # sequence.append([e[variableIndices[0], dimensions[1]], e[variableIndices[0], dimensions[0]]])

    for i in range(1, len(variableIndices)):
        currentSequence = []
        for j in reversed(range(len(dimensions))):
            # currentSequence.append(f"e[{variableIndices[i]}, {dimensions[j]}]")
            currentSequence.append(e[variableIndices[i], dimensions[j]])
            for s in range(len(sequence)):
                previousSequence = sequence[s]
                for k in range(len(previousSequence)):
                    # appendString = f"e[{variableIndices[i]}, {dimensions[j]}]"
                    appendProduct = e[variableIndices[i], dimensions[j]] * previousSequence[k]
                    # if previousSequence[k] != "":
                    # appendProduct *= previousSequence[k]
                    currentSequence.append(appendProduct)
        sequence.append(currentSequence)

    return sequence


def printSosTable(allTerms, p, e, variableIndices):
    depth = 0
    pExpressions = [allTerms[0]]
    eExpressions = [[]]
    seq = generate_sequence(e, variableIndices)

    for line in seq:
        for eProduct in line:
            num_factors = len(eProduct.args) if eProduct.is_Mul else 1

            if (num_factors > len(allTerms) - 1):
                continue

            # # print(num_factors)

            collectionTerm = eProduct
            collectionExpression = allTerms[num_factors]
            depth+=1
            collected_expr = collect(collectionExpression, collectionTerm).coeff(collectionTerm)


            # Skip if it's not in the expression
            if not collectionExpression.has(collectionTerm):
                print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                # print(f"The number of factors is {num_factors}")
                # print(f"The collection expressions is {collectionExpression}")
                # print(f"The p-term of that product has {len(collected_expr.as_ordered_terms())} terms:")
                # print((collected_expr))
                # print(f"Here it is simplified...")
                # print(simplify(collected_expr))
                print("------------------------------------- Skipping")
                continue

            # if collectionExpression in s or (-1 * collectionExpression) in s:
            # if any(collectionExpression.equals(e) or (-collectionExpression).equals(e) for e in s):
                # print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                # print("------------------------------------- Skipping p-term is already there")


            if isInArray(collected_expr, pExpressions):
                print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                print(f"------------------------------------- Skipping {collected_expr} p-term is already there")
                continue

            pExpressions.append(collected_expr)
            eExpressions.append(eProduct)

            # print(f"\n\nThe current depth is {depth}.")
            # print("The current e-product is:")
            # print(eProduct)
            # print("The collected expression is:")
            # print(collected_expr)
            # print(f"Current num_factors = {num_factors}, compared to {len(allTerms) - 1}")

            print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
            print(f"The p-term of that product has {len(collected_expr.as_ordered_terms())} terms:")
            print((collected_expr))
            print(f"Here it is simplified...")
            print(simplify(collected_expr))

            if (num_factors == len(allTerms) - 1):
                print("FOUND THE CONSTANT FACTOR!!!")
                return (pExpressions, eExpressions)

    return (pExpressions, eExpressions)


def evaluateExpresisonSign(pExpressions, eExpressions, allTerms, p, e, variables, indexSubstitution, pl, pi, pv, pj, pu, pk):
    print(f"The final depth is {len(pExpressions)}")

    for t in range(0, len(pExpressions)):
        print(f"\n\n-------------------------------------------- At depth {t}")
        print(f"The p-expression is:")
        print(f"{(pExpressions[t])}")
        print(f"The p-expression evaluation is:")
        pExpressionValue = evaluateExpression(pExpressions[t], p, variables, indexSubstitution, pi, pj, pk, pl, pu, pv)
        print(f"{pExpressionValue}")
        print(f"The e-expression is:")
        print(f"{eExpressions[t]}")

        if (pExpressionValue != 0):
            print("Reached an answer.")
            return sign(pExpressionValue)

    # expr3 = expr2.subs({
        # p[0,1]: 4,
        # p[1,1]: 5,
        # p[2,1]: 6,
        # p[3,1]: 7
    # })  # Step 2


