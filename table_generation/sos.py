# This code implements the symbolic perturbation scheme from the following paper
# 
# Edelsbrunner, H. and Mücke, E.P., 1990. Simulation of simplicity: a technique to cope with degenerate cases in geometric algorithms. ACM Transactions on Graphics (tog), 9(1), pp.66-104.

from sympy import simplify, expand, collect, Indexed

def getEvaluationTable(expression, e, variables):

    expression = simplify(expand(expression))
    expressionTermsOrdered = expression.as_ordered_terms()

    # Break up the expression into terms and order them by the number of e-variables in them.
    # 
    # the upper limit on the range is an arbitrary large number, the loops really is only supposed to stop when we break
    allTerms = []
    for index in range(0, 1000000):
        terms = sum([
            t for t in expressionTermsOrdered
            if count_indexed_with_base(t, e) == index
        ])

        # print(f"Index: {index}, terms: {terms}")

        if (terms == 0):
            break

        allTerms.append(terms)

    return computeEvaluationTable(allTerms, e, variables)

def computeEvaluationTable(allTerms, e, variableIndices):
    depth = 0
    pExpressions = [allTerms[0]]
    eExpressions = [[]]
    seq = generate_sequence(e, variableIndices)

    for line in seq:
        for eProduct in line:

            num_factors = len(eProduct.args) if eProduct.is_Mul else 1

            if (num_factors > len(allTerms) - 1):
                continue

            collectionTerm = eProduct
            collectionExpression = allTerms[num_factors]
            depth+=1
            collected_expr = collect(collectionExpression, collectionTerm).coeff(collectionTerm)


            # Skip if it's not in the expression
            if not collectionExpression.has(collectionTerm):
                # print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                # print(f"The number of factors is {num_factors}")
                # print(f"The collection expressions is {collectionExpression}")
                # print(f"The p-term of that product has {len(collected_expr.as_ordered_terms())} terms:")
                # print((collected_expr))
                # print(f"Here it is simplified...")
                # print(simplify(collected_expr))
                # print("------------------------------------- Skipping")
                continue

            # if collectionExpression in s or (-1 * collectionExpression) in s:
            # if any(collectionExpression.equals(e) or (-collectionExpression).equals(e) for e in s):
                # print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                # print("------------------------------------- Skipping p-term is already there")


            if isInArray(collected_expr, pExpressions):
                # print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
                # print(f"------------------------------------- Skipping {collected_expr} p-term is already there")
                continue

            pExpressions.append(collected_expr)
            eExpressions.append(eProduct)

            # print(f"\n\nThe current depth is {depth}.")
            # print("The current e-product is:")
            # print(eProduct)
            # print("The collected expression is:")
            # print(collected_expr)
            # print(f"Current num_factors = {num_factors}, compared to {len(allTerms) - 1}")

            # print(f"\n\nThe epsilon product with order {depth}: {(eProduct)}")
            # print(f"The p-term of that product has {len(collected_expr.as_ordered_terms())} terms:")
            # print((collected_expr))
            # print(f"Here it is simplified...")
            # print(simplify(collected_expr))

            if (num_factors == len(allTerms) - 1):
                # print("FOUND THE CONSTANT FACTOR!!!")
                return (pExpressions, eExpressions)

    return (pExpressions, eExpressions)



def isInArray(expression, expressionArray):
    for e in expressionArray:
        # print(f"comparing {expression} and {e}")
        if expression == e or -expression == e:
            return True
    return False

def count_indexed_with_base(term, base):
    """Return the number of Indexed objects in a term with the given base."""
    return sum(1 for atom in term.atoms(Indexed) if atom.base == base)



def generate_sequence(e, variables, dimensions=(1, 2)):

    # Sort the indies of all variables
    # variables = sorted(variables)

    sequence = []
    # sequence.append([f"e[{variables[0]}, {dimensions[1]}]", f"e[{variables[0]}, {dimensions[0]}]"])

    firstRow = []
    for j in reversed(range(len(dimensions))):
        firstRow.append(e[variables[0], dimensions[j]])

    sequence.append(firstRow)
    # sequence.append([e[variables[0], dimensions[1]], e[variables[0], dimensions[0]]])

    for i in range(1, len(variables)):
        currentSequence = []
        for j in reversed(range(len(dimensions))):
            # currentSequence.append(f"e[{variables[i]}, {dimensions[j]}]")
            currentSequence.append(e[variables[i], dimensions[j]])
            for s in range(len(sequence)):
                previousSequence = sequence[s]
                for k in range(len(previousSequence)):
                    # appendString = f"e[{variables[i]}, {dimensions[j]}]"
                    appendProduct = e[variables[i], dimensions[j]] * previousSequence[k]
                    # if previousSequence[k] != "":
                    # appendProduct *= previousSequence[k]
                    currentSequence.append(appendProduct)
        sequence.append(currentSequence)

    return sequence


