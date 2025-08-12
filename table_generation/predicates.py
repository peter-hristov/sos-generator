from sympy import Matrix, expand, IndexedBase, Symbol

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

def parametrizeAndOrder(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

    # p_l -> p_i
    r = Matrix([pi1 - pl1, pi2 - pl2])

    # p_v -> p_j
    s1 = Matrix([pj1 - pv1, pj2 - pv2])
    s2 = Matrix([pk1 - pu1, pk2 - pu2])

    # q1 = Matrix([pl1 - pv1, pl2 - pv2])
    # q2 = Matrix([pl1 - pu1, pl2 - pu2])

    q1 = Matrix([pv1 - pl1, pv2 - pl2])
    q2 = Matrix([pu1 - pl1, pu2 - pl2])

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



def dualizeAndOrient(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):
    # p_li
    p_l = (pl1, pl2)
    p_i = (pi1, pi2)

    # p_vj
    p_v = (pv1, pv2)
    p_j = (pj1, pj2)
    p_u = (pu1, pu2)
    p_k = (pk1, pk2)


    # p_li homogenious
    p_li = (
        p_i[1] - p_l[1],
        p_l[0] - p_i[0],
        p_i[0] * p_l[1] - p_l[0] * p_i[1]
    )

    # p_vj homogenious
    p_vj = (
        p_j[1] - p_v[1],
        p_v[0] - p_j[0],
        p_j[0] * p_v[1] - p_v[0] * p_j[1]
    )

    # p_uk
    p_uk = (
        p_k[1] - p_u[1],
        p_u[0] - p_k[0],
        p_k[0] * p_u[1] - p_u[0] * p_k[1]
    )

    return orientationTestHomogenious(p_li, p_vj, p_uk)
