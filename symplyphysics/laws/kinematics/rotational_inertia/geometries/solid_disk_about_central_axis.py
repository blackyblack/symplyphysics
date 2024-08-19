from sympy import Eq, pi
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
from symplyphysics.laws.kinematics.rotational_inertia import rotational_inertia_cylindrical_integral as integral_law
from symplyphysics.definitions import density_from_mass_volume as density_def

# Description
## A solid disk (cylinder) rotates about its central axis (axis of cylindrical symmetry).

# Law: I = 1/2 * M * R**2
## I - rotational inertia
## M - mass of disk
## R - radius of disk

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
radius = Symbol("radius", units.length)
disk_mass = clone_symbol(symbols.basic.mass, "disk_mass")

law = Eq(rotational_inertia, disk_mass * radius**2 / 2)

# Derive law from general integral in cylindrical coordinates

length = Symbol("length", units.length)
volume = pi * radius**2 * length

density = density_def.definition.rhs.subs({
    density_def.mass: disk_mass,
    density_def.volume: volume,
})

rotational_inertia_derived = integral_law.law.rhs.subs({
    integral_law.density(integral_law.radius, integral_law.polar_angle, integral_law.height):
        density,
    integral_law.radius_start:
        0,
    integral_law.radius_end:
        radius,
    integral_law.polar_angle_start:
        0,
    integral_law.polar_angle_end:
    2 * pi,
    integral_law.height_start:
        0,
    integral_law.height_end:
        length,
}).doit()

assert expr_equals(rotational_inertia_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=disk_mass, radius_=radius)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        disk_mass: mass_,
        radius: radius_,
    })
    return Quantity(result)
