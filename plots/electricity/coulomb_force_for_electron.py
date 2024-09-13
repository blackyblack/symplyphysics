#!/usr/bin/env python3
from sympy import Symbol
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, quantities
from symplyphysics.core.convert import convert_to_si
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as coulomb_law

CHARGE_OF_ELECTRON = -1.6 * 1E-19
VACUUM_PERMITTIVITY = convert_to_si(quantities.vacuum_permittivity)

print(f"Formula is:\n{print_expression(coulomb_law.law)}")

distance = Symbol("distance")

force_equation = coulomb_law.law.subs({
    coulomb_law.first_charge: CHARGE_OF_ELECTRON,
    coulomb_law.second_charge: CHARGE_OF_ELECTRON,
    coulomb_law.distance: distance,
    quantities.vacuum_permittivity: VACUUM_PERMITTIVITY,
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
