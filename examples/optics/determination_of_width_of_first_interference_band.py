#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.laws.optics import interference_from_two_slits as two_slits_law
from symplyphysics.laws.optics import interference_maximum as maximum_law

# Example from https://tasksall.ru/ovr0/watchbookfominyear1997num6.707.htm

distance_to_pattern = Symbol("distance_to_pattern")
distance_between_slits = Symbol("distance_between_slits")
wavelength = Symbol("wavelength")

number_of_maximum = Symbol("number_of_maximum")
optical_distance_difference = Symbol("optical_distance_difference")

distance_between_maxima = Symbol("distance_between_maxima")

optical_distance_difference_eq = two_slits_law.law.subs({
    two_slits_law.distance_to_picture: distance_to_pattern,
    two_slits_law.distance_between_slits: distance_between_slits,
    two_slits_law.travel_difference: optical_distance_difference
})
optical_distance_difference_value = solve(optical_distance_difference_eq, optical_distance_difference, dict=True)[0][optical_distance_difference]
maximum_condition_eq = maximum_law.law.subs({
    maximum_law.wave_length: wavelength,
    maximum_law.travel_difference: optical_distance_difference_value,
    maximum_law.number_maximum: number_of_maximum
})

coordinate_value = solve(maximum_condition_eq, two_slits_law.coordinate, dict=True)[0][two_slits_law.coordinate]

coordinate_for_first_maximum_value = coordinate_value.subs({
    number_of_maximum: 1
})
coordinate_for_second_maximum_value = coordinate_value.subs({
    number_of_maximum: 2
})

# The distance between the maxima is the difference in the coordinates
# of these maxima in the interference pattern
# dx = x_2 - x_1
distance_between_maxima_value = coordinate_for_second_maximum_value - coordinate_for_first_maximum_value
answer = Eq(distance_between_maxima, distance_between_maxima_value)
print(f"Equation for distance between first and second maxima:\n{print_expression(answer)}")

distance_between_maxima_ = answer.subs({
    distance_to_pattern: Quantity(1.5 * units.meters),
    distance_between_slits: Quantity(0.2 * prefixes.milli * units.meters),
    wavelength: Quantity(500 * prefixes.nano * units.meters),
}).rhs
print(f"distance between first and second maxima is: {print_expression(
    convert_to(Quantity(distance_between_maxima_), prefixes.milli * units.meters).evalf(3)
)} mm")
