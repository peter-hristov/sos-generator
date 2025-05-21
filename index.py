# Build the matrix
from sympy import symbols, Matrix, collect, expand, IndexedBase, Indexed, simplify

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

def count_indexed_with_base(term, base):
    """Return the number of Indexed objects in a term with the given base."""
    return sum(1 for atom in term.atoms(Indexed) if atom.base == base)

def generate_sequence(varNames, dimensions=(1, 2, 3)):
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

def getHomogeniousExpression():
    # Homogenious

    # p_li
    p_l = (p[l, 0] + e[l, 0], p[l, 1] + e[l, 1])
    p_i = (p[i, 0] + e[i, 0], p[i, 1] + e[i, 1])

    # p_vj
    p_v = (p[v, 0] + e[v, 0], p[v, 1] + e[v, 1])
    p_j = (p[j, 0] + e[j, 0], p[j, 1] + e[j, 1])

    # p_uk
    p_u = (p[u, 0] + e[u, 0], p[u, 1] + e[u, 1])
    p_k = (p[k, 0] + e[k, 0], p[k, 1] + e[k, 1])


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





p = IndexedBase('p')
e = IndexedBase('e')

# Symbolic indices
variables = symbols('i j k l u v')
i, j, k, l, u, v = variables

# Affine
det = orientationTestHomogenious(
        (p[i, 1] + e[i, 1], p[i, 2] + e[i, 2], p[i, 3] + e[i, 3]), 
        (p[j, 1] + e[j, 1], p[j, 2] + e[j, 2], p[j, 3] + e[j, 3]), 
        (p[k, 1] + e[k, 1], p[k, 2] + e[k, 2], p[k, 3] + e[k, 3])
        )

# det = getHomogeniousExpression()

# Optional: collect by e_i1 or any expression
# det_expanded = expand(det)
# det_collected = collect(det_expanded, e_i1)

expression = simplify(expand(det))
expressionTermsOrdered = expression.as_ordered_terms()

print(f"Total terms {len(expressionTermsOrdered)}")


allTerms = []
for i in range(0, 100):
    terms = sum([
        t for t in expressionTermsOrdered
        if count_indexed_with_base(t, e) == i
    ])

    if (terms == 0):
        break

    allTerms.append(terms)

print("Here are all the terms in the expression:")
for t in expressionTermsOrdered:
    print(t)

print("\n\nHere are all the terms grouped by the number of epsilons:")
for i in range(0, len(allTerms)):
    print(i)
    print(allTerms[i])

print(f"We have this many types of mixed expressions {len(allTerms)}")

def printSosTable():
    print("\n\nHere are all the sequences:")
    depth = 0
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
                print("FOUND THE CONSTANT FATOR!!!")
                return

printSosTable()
