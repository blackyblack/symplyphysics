#!/usr/bin/env python3
from sympy import solve, symbols, Symbol
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, Quantity, convert_to)

from symplyphysics.laws.electricity import force_from_charge_and_distance as coulomb_law

ELECTROSTATIC_CONSTANT = 9 * 1E9
CHARGE_OF_ELECTRON = -1.6 * 1E-19

print(f"Formula is:\n{coulomb_law.print_law()}")

distance = Symbol("distance")

force_equation = coulomb_law.law.subs({
    coulomb_law.first_charge: CHARGE_OF_ELECTRON,
    coulomb_law.second_charge: CHARGE_OF_ELECTRON,
    coulomb_law.distance: distance,
    coulomb_law.units.coulomb_constant: ELECTROSTATIC_CONSTANT
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
