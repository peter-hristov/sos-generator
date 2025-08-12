import random
from sympy import symbols, Rational, simplify

def generateRandomRationalCirclePoint(randomRange):
    # Define symbols
    x_sym, y_sym, t = symbols("x y t")

    # Parametrize x and y
    x_expr = (1 - t**2) / (1 + t**2)
    y_expr = (2 * t) / (1 + t**2)

    a = random.randint(randomRange[0], randomRange[1])
    b = random.randint(randomRange[0], randomRange[1])

    # Define a rational value of t
    t_val = Rational(a, b)  # This is t = 3/4

    # Evaluate x and y at t = 3/4
    x_val = simplify(x_expr.subs(t, t_val))
    y_val = simplify(y_expr.subs(t, t_val))

    # Print results
    assert x_val**2 + y_val**2 == 1

    return [x_val, y_val]


# Make sure we generate a valid rational (not division by zero)
def safe_rand_scale(randomRange):
    while True:
        num = random.randint(*randomRange)
        denom = random.randint(*randomRange)
        if denom != 0:
            return Rational(num, denom)

# Scale and offset each point
def transform_point(p, offsetX, offsetY, randomRange):
    scale = safe_rand_scale(randomRange)
    return [p[0]*scale + offsetX, p[1]*scale + offsetY]
    # return [p[0]*scale, p[1]*scale]
    # return [p[0], p[1]]

def generateSegments(randomRange):
    # Generate original points as lists
    pl = generateRandomRationalCirclePoint(randomRange)
    pi = [-pl[0], -pl[1]]

    pv = generateRandomRationalCirclePoint(randomRange)
    pj = [-pv[0], -pv[1]]

    pu = generateRandomRationalCirclePoint(randomRange)
    pk = [-pu[0], -pu[1]]

    # Offset center
    centerOffsetX = safe_rand_scale(randomRange)
    centerOffsetY = safe_rand_scale(randomRange)

    # Transform all points
    pl = transform_point(pl, centerOffsetX, centerOffsetY, randomRange)
    pi = transform_point(pi, centerOffsetX, centerOffsetY, randomRange)
    pv = transform_point(pv, centerOffsetX, centerOffsetY, randomRange)
    pj = transform_point(pj, centerOffsetX, centerOffsetY, randomRange)
    pu = transform_point(pu, centerOffsetX, centerOffsetY, randomRange)
    pk = transform_point(pk, centerOffsetX, centerOffsetY, randomRange)

    return pl, pi, pv, pj, pu, pk

def generateColinearPoints(randomRange):
    # Generate original points as lists
    pi = generateRandomRationalCirclePoint(randomRange)
    pj = [0, 0]
    pk = [-pi[0], -pi[1]]

    # Offset center
    centerOffsetX = safe_rand_scale(randomRange)
    centerOffsetY = safe_rand_scale(randomRange)

    # Transform all points
    pi = transform_point(pi, centerOffsetX, centerOffsetY, randomRange)
    pj = transform_point(pj, centerOffsetX, centerOffsetY, randomRange)
    pk = transform_point(pk, centerOffsetX, centerOffsetY, randomRange)

    return pi, pj, pk

def evaluateExpression(expression, p, variables, indexSubstitution, p1, p2, p3, p4, p5, p6):
    i, j, k, l, u, v = variables

    expression_index_subs = expression.subs(indexSubstitution)

    valueSubs = {
            p[indexSubstitution[i], 1]: p1[0], p[indexSubstitution[i], 2]: p1[1],
            p[indexSubstitution[j], 1]: p2[0], p[indexSubstitution[j], 2]: p2[1],
            p[indexSubstitution[k], 1]: p3[0], p[indexSubstitution[k], 2]: p3[1],
            p[indexSubstitution[l], 1]: p4[0], p[indexSubstitution[l], 2]: p4[1],
            p[indexSubstitution[u], 1]: p5[0], p[indexSubstitution[u], 2]: p5[1],
            p[indexSubstitution[v], 1]: p6[0], p[indexSubstitution[v], 2]: p6[1],
            }

    expression_value_subs = expression_index_subs.subs(valueSubs)

    return expression_value_subs

