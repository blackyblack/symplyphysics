#!/usr/bin/env python3
r"""
Two slits are separated by :math:`0.2 \, \text{mm}` and are located :math:`1.5 \, \text{m}` from a
screen. A beam of monochromatic light (:math:`\lambda = 500 \, \text{nm}`) from a distant source
falls on the slits. Find the distance between adjacent interference fringes.

**Source:** `(in Russian) <https://tasksall.ru/ovr0/watchbookfominyear1997num6.707.htm>`__.
"""

from sympy import solve, Symbol
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to, assert_equal
from symplyphysics.laws.optics import (
    interference_from_two_slits as two_slits_law,
    interference_maximum as interference_maximum_law,
)

distance_to_picture = two_slits_law.distance_to_picture
distance_between_slits = two_slits_law.distance_between_slits
wavelength = interference_maximum_law.wavelength
maximum_number = interference_maximum_law.integer_factor
optical_distance_difference = two_slits_law.travel_difference
distance_between_maxima = Symbol("D", real=True)

optical_distance_difference_expr = solve(two_slits_law.law, optical_distance_difference)[0]

maximum_condition_eqn = interference_maximum_law.law.subs(
    interference_maximum_law.travel_difference,
    optical_distance_difference_expr,
)
maximum_coordinate_expr = solve(maximum_condition_eqn, two_slits_law.coordinate)[0]

first_maximum_coordinate_expr = maximum_coordinate_expr.subs(maximum_number, 1)
second_maximum_coordinate_expr = maximum_coordinate_expr.subs(maximum_number, 2)

# The distance between the maxima is the difference in their coordinates in the interference pattern:
# `dx = x_2 - x_1`
distance_between_maxima_expr = second_maximum_coordinate_expr - first_maximum_coordinate_expr

print(
    "Equation for distance between first and second maxima:",
    print_expression(distance_between_maxima_expr),
    sep="\n",
)

distance_between_maxima_subs = distance_between_maxima_expr.subs({
    distance_to_picture: Quantity(1.5 * units.meters),
    distance_between_slits: Quantity(0.2 * prefixes.milli * units.meters),
    wavelength: Quantity(500 * prefixes.nano * units.meters),
})
assert_equal(distance_between_maxima_subs, 3.75 * units.milli * units.meter)

distance_between_maxima_float = convert_to(
    distance_between_maxima_subs,
    units.milli * units.meter,
).evalf(3)

print("distance between first and second maxima is:", distance_between_maxima_float, "mm")
