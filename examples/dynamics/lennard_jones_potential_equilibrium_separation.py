from sympy import solve, Eq, symbols
from symplyphysics import (
    print_expression,
    vector_magnitude,
    convert_to,
    units,
    Quantity,
)
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.points import cartesian_point
from symplyphysics.laws.dynamics.fields import (
    conservative_force_is_gradient_of_potential_energy as potential_law,
)

# Description
## The potential energy of attraction between atoms and molecules can be modeled in the form of the
## Lennard-Jones potential: U = A/r^12 - B/r^6, where A and B are positive constants. For argon,
## the corresponding coefficients are A = 1.58e-134 J*m^12, and B = 1.02e-77 J*m^6. Find the equilibrium
## separation, i.e. the distance between the atoms at which the force on each atom is zero.

# Reference frame:
## We assume two atoms lying on the x axis, and the first atom is located at the origin.
## Thus, the position of the second atom is described by its position on the x axis.

x = symbols("x", real=True)
A, B = symbols("A B", positive=True)

values = {
    A: Quantity(1.58e-134 * units.joule * units.meter**12),
    B: Quantity(1.02e-77 * units.joule * units.meter**6),
}


def potential_energy_function(point):
    return A / point.x**12 - B / point.x**6


potential_energy_field = ScalarField(potential_energy_function)
intermolecular_force_vector = potential_law.law(potential_energy_field)

intermolecular_force_magnitude = (
    vector_magnitude(intermolecular_force_vector)
    .subs(potential_energy_field.coordinate_system.coord_system.base_scalars()[0], x)
    .simplify()
)

# Equation has two solutions, first one is negative
equilibrium_separation = solve(intermolecular_force_magnitude, x)[1]

equilibrium_separation_in_nm = convert_to(
    Quantity(equilibrium_separation.subs(values)),
    units.nanometer,
).evalf(3)

print(f"Formula of equilibrium:\n{print_expression(equilibrium_separation)}\n")
print(f"Equilibrium separation for argon atoms is {equilibrium_separation_in_nm} nm.")
