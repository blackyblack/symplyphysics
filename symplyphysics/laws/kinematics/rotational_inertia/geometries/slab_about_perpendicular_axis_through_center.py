from sympy import Eq, sqrt
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
from symplyphysics.laws.kinematics.rotational_inertia import (
    rotational_inertia_cartesian_integral as integral_law,)
from symplyphysics.definitions import density_from_mass_volume as density_def

# Description
## A solid slab is rotating about an axis that is perpendicular to the slab and goes through its center.

# Law: I = 1/12 * M * (a**2 + b**2)
## I - rotational inertia of slab
## M - mass of slab
## a, b - sizes of slab in the directions perpendicular to rotation axis

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
length = Symbol("length", units.length)
width = Symbol("width", units.length)
slab_mass = clone_symbol(symbols.basic.mass, "slab_mass")

law = Eq(rotational_inertia, slab_mass * (length**2 + width**2) / 12)

# Derive this law from the integral definition of rotational inertia in cartesian coordiantes.
# Condition: density of slab is constant.

# Reference frame:
## z-axis is parallel to the rotational axis in question (height of slab)
## x-axis and y-axis are perpendicular to the rotational axis (length and width of slab)

height = Symbol("height", units.length)
volume = length * width * height

density = density_def.definition.rhs.subs({
    density_def.mass: slab_mass,
    density_def.volume: volume,
})

distance_to_axis = sqrt(integral_law.x**2 + integral_law.y**2)

rotational_inertia_derived = integral_law.law.rhs.subs({
    integral_law.density(integral_law.x, integral_law.y, integral_law.z): density,
    integral_law.distance_to_axis(integral_law.x, integral_law.y, integral_law.z): distance_to_axis,
    integral_law.x_start: -1 * length / 2,
    integral_law.y_start: -1 * width / 2,
    integral_law.z_start: 0,
    integral_law.x_end: length / 2,
    integral_law.y_end: width / 2,
    integral_law.z_end: height,
}).doit().simplify()

assert expr_equals(law.rhs, rotational_inertia_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=slab_mass, length_=length, width_=width)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity, width_: Quantity) -> Quantity:
    result = law.rhs.subs({
        slab_mass: mass_,
        length: length_,
        width: width_,
    })
    return Quantity(result)
