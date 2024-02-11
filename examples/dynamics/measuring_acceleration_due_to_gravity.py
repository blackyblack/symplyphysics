#!/usr/bin/env python3

from sympy import Symbol as SymSymbol, solve
from symplyphysics import units, print_expression
from symplyphysics.laws.dynamics import (
    period_of_physical_pendulum as pendulum_law,
)
from symplyphysics.laws.kinematic.rotational_inertia.geometries import (
    thin_rod_about_axis_through_center_perpendicular_to_length as rod_inertia_law
)
from symplyphysics.laws.kinematic.rotational_inertia import (
    rotational_inertia_about_axis_and_through_center_of_mass as parallel_axis_theorem,
)

# Description
## Take a pendulum to be a uniform rod suspended from one end. For such a pendulum,
## the distance between the pivot point and the center of mass is half of its length.
## We can estimate the value of the acceleration due to gravity by measuring the period
## of oscillations of this pendulum.

rod_length = SymSymbol("rod_length")
distance_to_pivot = rod_length / 2

rod_rotational_inertia_through_com = rod_inertia_law.law.rhs.subs(
    rod_inertia_law.length, rod_length
)

rod_rotational_inertia_through_pivot = parallel_axis_theorem.law.rhs.subs({
    parallel_axis_theorem.rotational_inertia_through_com: rod_rotational_inertia_through_com,
    parallel_axis_theorem.mass: rod_inertia_law.mass,
    parallel_axis_theorem.distance_between_axes: distance_to_pivot,
})

acceleration_due_to_gravity = solve(
    pendulum_law.law, units.acceleration_due_to_gravity
)[0].subs({
    pendulum_law.rotational_inertia: rod_rotational_inertia_through_pivot,
    pendulum_law.pendulum_mass: rod_inertia_law.mass,
    pendulum_law.distance_to_pivot: distance_to_pivot,
})

print(
    "Acceleration due to gravity from period of pendulum oscillations:\n"
    f"{print_expression(acceleration_due_to_gravity)}\n"
)
