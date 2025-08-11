from sympy import symbols, IndexedBase

# Define indexed base
p = IndexedBase('p')

# Define symbolic indices
i, j = symbols('i j')

# Construct expression
expr = p[i, 1] + p[j, 2]
print("Original symbolic expression:")
print(expr)

# Substitute specific index values for i and j
expr_index_subs = expr.subs({i: 0, j: 1})
print("\nAfter substituting indices i=0, j=1:")
print(expr_index_subs)

# Now substitute actual numerical values for the indexed elements
# For example: p[0,1] = 5, p[1,2] = 3
expr_value_subs = expr_index_subs.subs({p[0, 1]: 5, p[1, 2]: 3})
print("\nAfter substituting actual values:")
print(expr_value_subs)
