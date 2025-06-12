#!/usr/bin/env python3

from sympy import solve, symbols as sym_symbols
from symplyphysics import print_expression, convert_to, units, Quantity
from symplyphysics.laws.dynamics.fields import (
    conservative_force_is_gradient_of_potential_energy as potential_law,)

from symplyphysics.core.experimental.vectors import VectorNorm
from symplyphysics.core.experimental.coordinate_systems import (CoordinateScalar, CARTESIAN,
    CoordinateVector)

# Description
## The potential energy of attraction between atoms and molecules can be modeled in the form of the
## Lennard-Jones potential: U = A/r^12 - B/r^6, where A and B are positive constants. For argon,
## the corresponding coefficients are A = 1.58e-134 J*m^12, and B = 1.02e-77 J*m^6. Find the equilibrium
## separation, i.e. the distance between the atoms at which the force on each atom is zero.

# Reference frame:
## We assume two atoms lying on the x axis, and the first atom is located at the origin.
## Thus, the position of the second atom is described by its position on the x axis.

x, _, _ = CARTESIAN.base_scalars
A, B = sym_symbols("A B", positive=True)

values = {
    A: Quantity(1.58e-134 * units.joule * units.meter**12),
    B: Quantity(1.02e-77 * units.joule * units.meter**6),
}

potential_energy_field = CoordinateScalar(A / x**12 - B / x**6, CARTESIAN)

intermolecular_force_vector = potential_law.law.rhs.subs(
    potential_law.potential_energy,
    potential_energy_field,
).doit()
intermolecular_force_vector = CoordinateVector.from_expr(intermolecular_force_vector)

intermolecular_force_magnitude = VectorNorm(intermolecular_force_vector).simplify()

# Equation has two solutions, first one is negative
equilibrium_separation = solve(intermolecular_force_magnitude, x)[1]

equilibrium_separation_in_nm = convert_to(
    Quantity(equilibrium_separation.subs(values)),
    units.nanometer,
).evalf(3)

print(f"Formula of equilibrium:\n{print_expression(equilibrium_separation)}\n")
print(f"Equilibrium separation for argon atoms is {equilibrium_separation_in_nm} nm.")
