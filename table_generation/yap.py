from sympy import diff
import itertools


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
