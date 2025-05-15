from sympy import symbols

# Define symbols like p_{i,1}, e_{j,2}, etc.
p_i1, p_i2, p_i3 = symbols('p_i1 p_i2 p_i3')
e_i1, e_i2, e_i3 = symbols('e_i1 e_i2 e_i3')
p_j1, p_j2, p_j3 = symbols('p_j1 p_j2 p_j3')
e_j1, e_j2, e_j3 = symbols('e_j1 e_j2 e_j3')
p_k1, p_k2, p_k3 = symbols('p_k1 p_k2 p_k3')
e_k1, e_k2, e_k3 = symbols('e_k1 e_k2 e_k3')

# Build the matrix
from sympy import Matrix

M = Matrix([
    [p_i1 + e_i1, p_i2 + e_i2, p_i3 + e_i3],
    [p_j1 + e_j1, p_j2 + e_j2, p_j3 + e_j3],
    [p_k1 + e_k1, p_k2 + e_k2, p_k3 + e_k3],
])

# Compute determinant
det = M.det()

# Optional: collect by e_i1 or any expression
from sympy import collect, expand
det_expanded = expand(det)
det_collected = collect(det_expanded, e_i1)

print(det_collected)
