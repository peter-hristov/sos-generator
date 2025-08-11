from sympy import sign, Add, Mul, Pow

def evaluateExpresisonSign(pExpressions, eExpressions, p, variables, indexSubstitution, pl, pi, pv, pj, pu, pk):
    # print(f"The maximum depth is {len(pExpressions)}")

    for t in range(0, len(pExpressions)):
        # print(f"\n\n-------------------------------------------- At depth {t}")
        # print(f"The p-expression is:")
        # print(f"{(pExpressions[t])}")
        # print(f"The p-expression evaluation is:")
        pExpressionValue = evaluateExpression(pExpressions[t], p, variables, indexSubstitution, pi, pj, pk, pl, pu, pv)
        # print(f"{pExpressionValue}")
        # print(f"The e-expression is:")
        # print(f"{eExpressions[t]}")

        if (pExpressionValue != 0):
            # print("Reached an answer.")
            return sign(pExpressionValue), t

    # We should never be here
    assert False

def count_ops(expr):
    """Count arithmetic operations in a SymPy expression."""
    if expr.is_Atom:
        return 0

    count = 0
    if isinstance(expr, Add):
        # n terms requires n - 1 additions
        count += len(expr.args) - 1
    elif isinstance(expr, Mul):
        # n terms requires n - 1 multiplications
        count += len(expr.args) - 1
    elif isinstance(expr, Pow):
        # one exponentiation
        count += 1

    # Recurse into arguments
    for arg in expr.args:
        count += count_ops(arg)

    return count

