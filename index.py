# Build the matrix
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify, ccode

p = IndexedBase('p', shape=(6, 2))
e = IndexedBase('e', shape=(6, 2))

# Symbolic indices
variables = symbols('i j k l u v')
i, j, k, l, u, v = variables

def evaluateExpression(expression, p1, p2, p3, p4, p5, p6):
    expression_index_subs = expression.subs({i: 0, j: 1, k: 2, l: 3, u: 4, v: 5})

    valueSubs = {
            p[0, 1]: p1[0], p[0, 2]: p1[1],
            p[1, 1]: p2[0], p[1, 2]: p2[1],
            p[2, 1]: p3[0], p[2, 2]: p3[1],
            p[3, 1]: p4[0], p[3, 2]: p4[1],
            p[4, 1]: p5[0], p[4, 2]: p5[1],
            p[5, 1]: p6[0], p[5, 2]: p6[1],
            }

    expression_value_subs = expression_index_subs.subs(valueSubs)

    return expression_value_subs


def isInArray(expression, expressionArray):
    for e in expressionArray:
        # print(f"comparing {expression} and {e}")
        if expression == e or -expression == e:
            return True
    return False

def orientationTest(p1, p2, p3):
    M = Matrix([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1],
    ])

    return M.det()

def orientationTestHomogenious(p1, p2, p3):
    M = Matrix([
        [p1[0], p1[1], p1[2]],
        [p2[0], p2[1], p2[2]],
        [p3[0], p3[1], p3[2]],
    ])

    return M.det()

def orientationTestHomogenious4D(p1, p2, p3, p4):
    M = Matrix([
        [p1[0], p1[1], p1[2], p1[3]],
        [p2[0], p2[1], p2[2], p2[3]],
        [p3[0], p3[1], p3[2], p3[3]],
        [p4[0], p4[1], p4[2], p4[3]],
    ])

    return M.det()

def count_indexed_with_base(term, base):
    """Return the number of Indexed objects in a term with the given base."""
    return sum(1 for atom in term.atoms(Indexed) if atom.base == base)

def generate_sequence(varNames, dimensions=(1, 2)):
    sequence = []
    # sequence.append([f"e[{varNames[0]}, {dimensions[1]}]", f"e[{varNames[0]}, {dimensions[0]}]"])

    firstRow = []
    for j in reversed(range(len(dimensions))):
        firstRow.append(e[varNames[0], dimensions[j]])

    sequence.append(firstRow)
    # sequence.append([e[varNames[0], dimensions[1]], e[varNames[0], dimensions[0]]])

    for i in range(1, len(varNames)):
        currentSequence = []
        for j in reversed(range(len(dimensions))):
            # currentSequence.append(f"e[{varNames[i]}, {dimensions[j]}]")
            currentSequence.append(e[varNames[i], dimensions[j]])
            for s in range(len(sequence)):
                previousSequence = sequence[s]
                for k in range(len(previousSequence)):
                    # appendString = f"e[{varNames[i]}, {dimensions[j]}]"
                    appendProduct = e[varNames[i], dimensions[j]] * previousSequence[k]
                    # if previousSequence[k] != "":
                    # appendProduct *= previousSequence[k]
                    currentSequence.append(appendProduct)
        sequence.append(currentSequence)
                    
                    

    return sequence

def dualizeAndOrient():
    # Homogenious

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

    return orientationTestHomogenious(p_li, p_vj, p_uk)


def parametrizeAndOrder():

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








print("Computing initial expression...")
# Affine

# det = orientationTest(
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

det = dualizeAndOrient()
# det = parametrizeAndOrder()


expression = simplify(expand(det))
expressionTermsOrdered = expression.as_ordered_terms()

print(f"Total terms {len(expressionTermsOrdered)}")


print("Arranging terms into different levels...")
allTerms = []
# the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
for index in range(0, 100000):
    terms = sum([
        t for t in expressionTermsOrdered
        if count_indexed_with_base(t, e) == index
    ])

    if (terms == 0):
        break

    allTerms.append(terms)

print("Here are all the terms in the expression:")
for t in expressionTermsOrdered:
    print(t)

print("\n\nHere are all the terms grouped by the number of epsilons:")
for t in range(0, len(allTerms)):
    print(t)
    print(allTerms[t])

print(f"We have this many types of mixed expressions {len(allTerms)}")

def printSosTable():
    print("\n\nHere are all the sequences:")
    depth = 0
    pExpressions = [allTerms[0]]
    eExpressions = [[]]

    seq = generate_sequence(variables)
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

pExpressions, eExpressions = printSosTable()
print(f"The final depth is {len(pExpressions)}")

# Blue
pl = (2, 9)
pi = (9, 2)

# Orange
pv = (2, 4)
pj = (10, 6)

# Green
pu = (4, 2)
pk = (9, 5)


for t in range(0, len(pExpressions)):
    print(f"\n\n-------------------------------------------- At depth {t}")
    print(f"The p-expression is:")
    print(f"{(pExpressions[t])}")
    print(f"The p-expression evaluation is:")
    pExpressionValue = evaluateExpression(pExpressions[t], pi, pj, pk, pl, pu, pv)
    print(f"{pExpressionValue}")
    print(f"The e-expression is:")
    print(f"{eExpressions[t]}")

    if (pExpressionValue != 0):
        print("Reaches an answer.")
        break



# expr3 = expr2.subs({
    # p[0,1]: 4,
    # p[1,1]: 5,
    # p[2,1]: 6,
    # p[3,1]: 7
# })  # Step 2


