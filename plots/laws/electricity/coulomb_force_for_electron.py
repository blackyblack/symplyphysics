#!/usr/bin/env python3
from sympy import solve, symbols, Symbol
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, Quantity, convert_to)

from symplyphysics.laws.electricity import force_from_charge_and_distance as coulomb_law

CHARGE_OF_ELECTRON = -1.6 * 1E-19

print(f"Formula is:\n{coulomb_law.print_law()}")

distance = Symbol("dictance")
electrostatic_constant = Symbol("electrostatic_constant")

force_equation = coulomb_law.law.subs({
    coulomb_law.first_charge: CHARGE_OF_ELECTRON,
    coulomb_law.second_charge: CHARGE_OF_ELECTRON,
    coulomb_law.distance: distance,
    coulomb_law.units.coulomb_constant:
    convert_to(Quantity(units.coulomb_constant),
               units.newton * (units.meters / units.coulomb) ** 2)
}).rhs

p1 = plot(force_equation, distance,
          line_color="blue",
          title="Coluomb Law",
          xlabel="r, m",
          ylabel="F, N",
          backend=MatplotlibBackend,
          legend=True,
          show=False)
p1.show()
