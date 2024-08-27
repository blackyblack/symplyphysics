#!/usr/bin/env python3
from sympy import Symbol
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as coulomb_law

ELECTROSTATIC_CONSTANT = 9 * 1E9
CHARGE_OF_ELECTRON = -1.6 * 1E-19

print(f"Formula is:\n{print_expression(coulomb_law.law)}")

distance = Symbol("distance")

force_equation = coulomb_law.law.subs({
    coulomb_law.first_charge: CHARGE_OF_ELECTRON,
    coulomb_law.second_charge: CHARGE_OF_ELECTRON,
    coulomb_law.distance: distance,
    units.coulomb_constant: ELECTROSTATIC_CONSTANT
}).rhs

p1 = plot(force_equation, (distance, 1, 5),
    line_color="blue",
    title="Coulomb Law",
    xlabel="r, m",
    ylabel="F, N",
    backend=MatplotlibBackend,
    label="Electrostatic force",
    legend=True,
    show=False)
p1.show()
