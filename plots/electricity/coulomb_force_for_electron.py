#!/usr/bin/env python3
"""
Plot the electrostatic force between two electrons as a function of distance
between them.
"""

from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, quantities, convert_to_si, units
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as coulomb_law

print(f"Formula is:\n{print_expression(coulomb_law.law)}")

distance = coulomb_law.distance

force_expression = coulomb_law.law.rhs.subs({
    coulomb_law.first_charge: -1 * quantities.elementary_charge,
    coulomb_law.second_charge: -1 * quantities.elementary_charge,
})
force_expression = evaluate_expression(force_expression)

min_distance = convert_to_si(1e-14 * units.meter)
max_distance = convert_to_si(5e-14 * units.meter)

p1 = plot(
    force_expression,
    (distance, min_distance, max_distance),
    line_color="blue",
    title="Coulomb Law",
    xlabel="r, m",
    ylabel="F, N",
    backend=MatplotlibBackend,
    label="Electrostatic force",
    legend=True,
    show=False,
)
p1.show()
