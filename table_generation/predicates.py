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

def parametrizeAndOrder(pl, pi, pv, pj, pu, pk):

    # p_l -> p_i
    r = Matrix([pi[0] - pl[0], pi[1] - pl[1]])

    # p_v -> p_j
    s1 = Matrix([pj[0] - pv[0], pj[1] - pv[1]])
    s2 = Matrix([pk[0] - pu[0], pk[1] - pu[1]])

    # q1 = Matrix([pl1 - pv1, pl2 - pv2])
    # q2 = Matrix([pl1 - pu1, pl2 - pu2])

    q1 = Matrix([pv[0] - pl[0], pv[1] - pl[1]])
    q2 = Matrix([pu[0] - pl[0], pu[1] - pl[1]])

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



def dualizeAndOrient(pl, pi, pv, pj, pu, pk):

    # p_li homogenious
    p_li = (
        pi[1] - pl[1],
        pl[0] - pi[0],
        pi[0] * pl[1] - pl[0] * pi[1]
    )

    # p_vj homogenious
    p_vj = (
        pj[1] - pv[1],
        pv[0] - pj[0],
        pj[0] * pv[1] - pv[0] * pj[1]
    )

    # p_uk
    p_uk = (
        pk[1] - pu[1],
        pu[0] - pk[0],
        pk[0] * pu[1] - pu[0] * pk[1]
    )

    return orientationTestHomogenious(p_li, p_vj, p_uk)
