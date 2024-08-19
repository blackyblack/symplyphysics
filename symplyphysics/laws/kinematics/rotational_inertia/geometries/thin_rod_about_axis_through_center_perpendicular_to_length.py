from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics.rotational_inertia.geometries import (
    slab_about_perpendicular_axis_through_center as slab_formula)

# Description
## A thin rod is rotating about an axis that passes through its center perpendicular to its length.

# Law: I = 1/12 * M * L**2
## I - rotational inertia of thin rod
## M - mass of rod
## L - length of rod

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
length = Symbol("length", units.length)
rod_mass = clone_symbol(symbols.basic.mass, "rod_mass")

law = Eq(rotational_inertia, rod_mass * length**2 / 12)

# Derive law from formula for a slab rotating about the axis perpendicular to its lenght and width
# passing through its center. The thin rod is a particular case of it, when the width of the slab
# approaches zero.

rotational_inertia_derived = slab_formula.law.rhs.subs({
    slab_formula.slab_mass: rod_mass,
    slab_formula.length: length,
    slab_formula.width: 0,
})

assert expr_equals(law.rhs, rotational_inertia_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=rod_mass, length_=length)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rod_mass: mass_,
        length: length_,
    })
    return Quantity(result)
