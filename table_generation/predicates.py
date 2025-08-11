from sympy import Matrix, expand, IndexedBase, Symbol

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

def orientationTestYap(p1, p2, p3):
    M = Matrix([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1],
    ])

    return M.det()


def parametrizeAndOrderYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):

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



def dualizeAndOrientYap(pl1, pl2, pi1, pi2, pv1, pv2, pj1, pj2, pu1, pu2, pk1, pk2):
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

    p = IndexedBase('p')
    return orientationTestHomogenious(p, p_li, p_vj, p_uk)

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


def dualizeAndOrientAllienz(p, p_l, p_i, p_v, p_j, p_u, p_k):

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



