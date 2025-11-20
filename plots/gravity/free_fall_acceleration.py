#!/usr/bin/env python3
"""
Plot the free fall acceleration on Earth as a function of elevation (i.e. perpendicular distance
to Earth's surface).
"""

from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, quantities, units
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as acceleration_law

print(
    f"Free fall acceleration `a` as a function of elevation `h`:",
    print_expression(acceleration_law.law),
    sep="\n",
    end="\n\n",
)

elevation = acceleration_law.elevation

gravity_constant_value = quantities.gravitational_constant
EARTH_MASS = 5.9722e24 * units.kilogram
EARTH_RADIUS = 6.371e6 * units.meter

free_fall_acceleration_expr = acceleration_law.law.rhs.subs({
    acceleration_law.planet_mass: EARTH_MASS,
    acceleration_law.planet_radius: EARTH_RADIUS,
    acceleration_law.elevation: elevation,
})
free_fall_acceleration_expr = evaluate_expression(free_fall_acceleration_expr)

p1 = plot(
    free_fall_acceleration_expr,
    (elevation, 0, 10e3),  # meters
    ylim=(9.74, 9.85),
    axis_center=(0.0, 9.75),
    line_color="red",
    title="Free fall acceleration $g$ versus elevation $h$",
    xlabel=r"$h, \text{m}$",
    ylabel=r"$g, \text{m}/\text{s}^2$",
    backend=MatplotlibBackend,
    show=False,
)
p1.show()
