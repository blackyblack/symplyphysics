#!/usr/bin/env python3
"""
Bricks are stacked one on top of another without any glue, so that each successive brick extends
slightly beyond the one below it. What is the maximum distance by which the right edge of the top
brick can extend beyond the right edge of the bottom brick that serves as the base of the entire
stack?
"""

import matplotlib.pyplot as plt
from sympy import Idx, solve, Eq, IndexedBase, Indexed, Sum, Wild, Expr, rot_axis3
from sympy.matrices.dense import DenseMatrix
from symplyphysics import symbols, print_expression, global_index, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    as_coordinate_vector)
from symplyphysics.core.experimental.vectors import VectorDot, VectorCross
from symplyphysics.core.experimental.solvers import vector_equals
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy_law
from symplyphysics.laws.kinematics.vector import center_of_mass_for_system_of_particles as center_of_mass_def

## Part 1. Creating the necessary symbols

# Let us number each block so that the topmost block is `1`, the one below it is `2`, and so on.

# Index of the system bound to the center of block `n`
system_index = Idx("n")

stop_block_index = Idx("k")

block_index = Idx("i")

# `C[n, i]` is the position of the center of block `i` in the system bound to the center of block
# `n`.
center = IndexedBase("C", (system_index, block_index), positive=True)

# `X[n, k]` is the position of the center of mass of blocks `1` through `k` in the system bound to
# the center of block `n`.
mass_center = IndexedBase("X", (system_index, stop_block_index), positive=True)

# `S[n, k] = Sum(X[n, i], (i, 1, k))` is the total sum of the centers of blocks from `1` to `k` in
# the system bound to the center of block `n`.
sum_of_centers = IndexedBase("S", (system_index, stop_block_index), positive=True)

## Part 2. Deriving the expression for the mass center in a one-dimensional case

block_mass = clone_as_symbol(symbols.mass, positive=True)

block_index_ranged = Idx("i", (1, stop_block_index))

e_x = CoordinateVector([1, 0, 0], CARTESIAN)
center_vector_expr = center[system_index, block_index_ranged] * e_x

mass_center_vector_expr = center_of_mass_def.law.rhs.subs({
    global_index: block_index_ranged,
}).subs({
    center_of_mass_def.mass[block_index_ranged]: block_mass,
    center_of_mass_def.position_vector[block_index_ranged]: center_vector_expr,
}).doit().simplify()
mass_center_vector_expr = as_coordinate_vector(mass_center_vector_expr)

assert vector_equals(VectorCross(mass_center_vector_expr, e_x), 0)
mass_center_expr = VectorDot(mass_center_vector_expr, e_x).doit()

sum_of_centers_eqn = Eq(
    sum_of_centers[system_index, stop_block_index],
    Sum(center[system_index, block_index_ranged], block_index_ranged),
)

# `X[n, k]`
mass_center_expr = solve(
    (Eq(mass_center[system_index, stop_block_index], mass_center_expr), sum_of_centers_eqn),
    (mass_center[system_index, stop_block_index], sum_of_centers_eqn.rhs),
    dict=True,
)[0][mass_center[system_index, stop_block_index]]

## Part 3. Deriving the relation for the mass centers.

# `X[n, k] = S[n, k] / k`. See derived result above.
mass_center_eqn = Eq(
    mass_center[system_index, stop_block_index],
    mass_center_expr,
)

# `S[n, k]`
sum_of_centers_expr = solve(
    mass_center_eqn,
    sum_of_centers[system_index, stop_block_index],
)[0]

# `S[n, n - 1]`
prev_sum_of_centers_in_curr_expr = sum_of_centers_expr.subs({
    stop_block_index: system_index - 1,
})

# `S[n, n]`. Note that `C[n, n] = 0`, which simply means that the center of block `n` is the
# origin of the system bound to it.
curr_sum_of_centers_in_curr_expr = (prev_sum_of_centers_in_curr_expr +
    center[system_index, system_index])
curr_sum_of_centers_in_curr_expr = curr_sum_of_centers_in_curr_expr.subs({
    center[system_index, system_index]: 0,
})

# XXX Apparently, `sympy` has an issue here. The `sum_of_centers` symbol gets replaced in the
# following substitution and we have to retrieve it ourselves.

# `X[n, n]`
curr_mass_center_in_curr_expr = mass_center_eqn.rhs.subs(stop_block_index, system_index)

(new_sum_of_centers,) = (atom for atom in curr_mass_center_in_curr_expr.atoms(IndexedBase)
    if str(atom) == str(sum_of_centers))

curr_mass_center_in_curr_expr = curr_mass_center_in_curr_expr.subs({
    new_sum_of_centers[system_index, system_index]: curr_sum_of_centers_in_curr_expr,
})

# `X[n + 1, n] = C[n + 1, n] + X[n, n]`. This represents the fact that the system bound to block
# `n` is shifted by distance `C[n + 1, n]` relative to the system bound to the center of block
# `n + 1`.
curr_mass_center_in_next_expr = (center[system_index + 1, system_index] +
    curr_mass_center_in_curr_expr)

# Finally, we obtain a recurrent relationship for the mass centers.
curr_mass_center_in_next_eqn = Eq(
    mass_center[system_index + 1, system_index],
    curr_mass_center_in_next_expr,
)

## Part 4. Deriving the stability condition for stacked blocks.

# Replace blocks from `1` to `n - 1` with a single block of the same mass, the same position of
# the mass center, and the same overhang relative to block number `n`. Now, we can rotate the top
# block around the pivot point (which is the top right corner of the bottom block) by a small
# angle and see if the potential energy of the top block decreases or increases.

bottom_block_length = clone_as_symbol(symbols.length, positive=True)
bottom_block_height = clone_as_symbol(symbols.height, positive=True)
top_mass = clone_as_symbol(symbols.mass, positive=True, display_symbol="M")
top_mass_center_x = mass_center[system_index, system_index - 1]
top_mass_center_y = clone_as_symbol(symbols.position, positive=True, display_symbol="Y[n, n - 1]")
angle = clone_as_symbol(symbols.angle, positive=True)

# In the system bound to the center of the bottom block:
pivot_vector_expr = CoordinateVector(
    [bottom_block_length / 2, bottom_block_height / 2, 0],
    CARTESIAN,
)
mass_center_vector_expr = CoordinateVector(
    [top_mass_center_x, top_mass_center_y, 0],
    CARTESIAN,
)


def matmul(matrix: DenseMatrix, vector: Expr) -> CoordinateVector:
    vector = as_coordinate_vector(vector)

    if vector == 0:
        return vector

    components = vector.components
    rotated_components = matrix * components
    rotated_vector = CoordinateVector(rotated_components, vector.system, vector.point)

    return rotated_vector


# For `rot_axis3`, positive argument means the rotation is clockwise as we need.
rotation_matrix = rot_axis3(angle)

# shift -> rotate -> unshift
rotated_center_vector_expr = (pivot_vector_expr +
    matmul(rotation_matrix, mass_center_vector_expr - pivot_vector_expr))

center_displacement_vector_expr = rotated_center_vector_expr - mass_center_vector_expr

e_z = CoordinateVector([0, 1, 0], CARTESIAN)
height_change_expr = VectorDot(center_displacement_vector_expr, e_z).doit()

potential_energy_change_expr = potential_energy_law.law.rhs.subs({
    potential_energy_law.mass: top_mass,
    potential_energy_law.height: height_change_expr,
})

# The angle of rotation is small, so we can apply the Taylor series expansion and cut off higher
# terms.
potential_energy_change_expr = potential_energy_change_expr.series(angle, 0, 2).removeO().simplify()

# The equilibrium is stable if the potential energy of the system increases when disturbed. That
# is, the equilibrium is the local minimum of the potential energy.
stability_condition = solve(potential_energy_change_expr >= 0, top_mass_center_x)

print("Let `", print_expression(bottom_block_length), "` be the length of each block.", sep="")
print(
    "Let `",
    print_expression(top_mass_center_x),
    "` be the horizontal position of the mass center of blocks from 1 to `n - 1` relative to the center of block `n`.",
    sep="",
    end="\n\n",
)

print(
    "Stability condition of the blocks `1` through `n`:",
    print_expression(stability_condition),
    sep="\n\n",
    end="\n\n",
)

## Part 5. Combine the results

block_length = bottom_block_length

# The stability condition is that `X[n, n - 1] <= l / 2` for all `n`. That is, the mass center of
# blocks from `1` to `n - 1` lies no further than the end of block `n`. The maximum overhang can
# be obtained when `X[n, n - 1] = l / 2`.

maximum_mass_center_expr = solve(
    Eq(stability_condition.lhs, stability_condition.rhs),
    top_mass_center_x,
)[0]

# XXX Same issue as above. For some reason the `mass_center` symbol gets replaced.

curr_mass_center_in_next_eqn_subs = curr_mass_center_in_next_eqn
wild_index = Wild("wild_index")

# Substituting `X[n, n - 1]` with `maximum_mass_center_expr`.
for atom in curr_mass_center_in_next_eqn_subs.atoms(Indexed):
    base = atom.args[0]
    if str(base) != str(mass_center):
        continue

    matched = atom.match(base[wild_index, wild_index - 1])
    if not matched:
        continue

    curr_mass_center_in_next_eqn_subs = curr_mass_center_in_next_eqn_subs.subs({
        atom: maximum_mass_center_expr,
    })

# `C[n + 1, n]`
curr_center_position_in_next_expr = solve(
    curr_mass_center_in_next_eqn_subs,
    center[system_index + 1, system_index],
)[0]

print(
    "The position of the center of block `n` relative to the center of block `n + 1`:",
    print_expression(curr_center_position_in_next_expr),
    sep="\n\n",
    end="\n\n",
)

block_count = symbols.positive_number

# `C[n, 1] = C[2, 1] + C[n, 2]
#          = C[2, 1] + C[3, 2] + C[n, 3]
#          = ...
#          = C[2, 1] + C[3, 2] + ... + C[n, n - 1]`
maximum_overhang_expr = Sum(
    curr_center_position_in_next_expr,
    (system_index, 1, block_count - 1),
).doit().simplify()

print(
    "The maximum overhang in a system of `N` blocks (counting the bottom block as well):",
    print_expression(maximum_overhang_expr),
    sep="\n\n",
    end="\n\n",
)

maximum_overhang_expr_subs = maximum_overhang_expr.subs(block_count, 4)
assert expr_equals(maximum_overhang_expr_subs, 11 * block_length / 12)

print(
    "For a system of `N = 4` blocks, the maximum overhang is",
    print_expression(maximum_overhang_expr.subs(block_count, 4)),
    sep="\n\n",
    end="\n\n",
)

## Part 6. Plot the results

block_counts = list(range(2, 100))

maximum_overhangs = [
    maximum_overhang_expr.subs({
    block_length: 1,
    block_count: block_count_,
    }).evalf() for block_count_ in block_counts
]

plt.semilogx(block_counts, maximum_overhangs, "b.-")
plt.title("Maximum overhang as a function of block count")
plt.xlabel("number of blocks")
plt.ylabel("maximum overhang per unit block length")
plt.grid(which="both")
plt.show()
