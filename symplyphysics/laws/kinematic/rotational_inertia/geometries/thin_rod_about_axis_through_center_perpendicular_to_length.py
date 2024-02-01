from sympy import Eq, Symbol as SymSymbol, solve, integrate
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic.rotational_inertia import (
    rotational_inertia_of_particle as rotational_inertia_def,
)

# Description
## A thin rod is rotating about an axis that passes through its center perpendicular to its length.

# Law: I = 1/12 * M * L**2
## I - rotational inertia of thin rod
## M - mass of rod
## L - length of rod

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
mass = Symbol("mass", units.mass)
length = Symbol("length", units.length)

law = Eq(rotational_inertia, mass * length**2 / 12)


# Derive law from definition of rotational inertia

# We will consider the thin rod as continuously distributed elements. We can calculate
# the total rotational inertia of the rod as the sum of rotational inertias of each of
# its parts, namely, by integrating over the rod's length.

infinitesimal_length = SymSymbol("infinitesimal_length", positive=True)
infinitesimal_mass = SymSymbol("infinitesimal_mass", positive=True)

# Because the rod is uniform, the ratio of mass to length is equal for all elements of
# the rod and for the rod itself. This equation basically states that the local linear
# density of the rod is constant and equal to the ratio of rod's mass to its length.

uniformity_eqn = Eq(infinitesimal_mass / infinitesimal_length, mass / length)

# Express the mass of the element through its length.
infinitesimal_mass = solve(uniformity_eqn, infinitesimal_mass)[0]

# The rod is oriented along the x-axis, and the origin is located in its center due to
# the fact that the rotational axis passes through the center of the rod.

x = SymSymbol("x", real=True)

infinitesimal_rotational_inertia = rotational_inertia_def.law.rhs.subs({
    rotational_inertia_def.radius: x,
    rotational_inertia_def.mass: infinitesimal_mass
})

rotational_inertia_derived = integrate(
    infinitesimal_rotational_inertia.subs(infinitesimal_length, 1),
    (x, -1 * length / 2, length / 2),
).simplify()

assert expr_equals(law.rhs, rotational_inertia_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, length_=length)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        length: length_,
    })
    return Quantity(result)
