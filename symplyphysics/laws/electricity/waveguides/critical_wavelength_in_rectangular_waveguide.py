from sympy import Eq, solve, sqrt
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## There is a critical wavelength. Signals with a wavelength greater than the critical one are attenuated and
## do not propagate in the waveguide.
## The critical frequency depends on the indices of the specific propagation mode and the dimensions of the
## waveguide cross-section.
## The first index shows how many half-wave lengths fit horizontally across the cross section. The second index
## shows how many half-wave lengths fit vertically across the cross section.

## Law is: L = 2 / sqrt((m / a)^2 + (n / b)^2), where
## L - critical wavelength in rectangular waveguide,
## m - first index,
## n - second index,
## a - width of the waveguide cross section,
## b - height of the waveguide cross section.

critical_wavelength = Symbol("critical_wavelength", units.length)

first_index = Symbol("first_index", dimensionless)
second_index = Symbol("second_index", dimensionless)
width = Symbol("width", units.length)
height = Symbol("height", units.length)

law = Eq(critical_wavelength, 2 / sqrt((first_index / width)**2 + (second_index / height)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(first_index_=first_index, second_index_=second_index, width_=width, height_=height)
@validate_output(critical_wavelength)
def calculate_critical_wavelength(first_index_: float, second_index_: float, width_: Quantity,
    height_: Quantity) -> Quantity:
    if first_index_ < 0 or second_index_ < 0:
        raise ValueError("Indexes should not be negative")
    result_velocity_expr = solve(law, critical_wavelength, dict=True)[0][critical_wavelength]
    result_expr = result_velocity_expr.subs({
        first_index: first_index_,
        second_index: second_index_,
        width: width_,
        height: height_
    })
    return Quantity(result_expr)
