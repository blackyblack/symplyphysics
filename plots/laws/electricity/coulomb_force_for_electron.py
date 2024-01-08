#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import (print_expression, units, Quantity, convert_to)

from symplyphysics.laws.electricity import force_from_charge_and_distance as coulomb_law

CHARGE_OF_ELECTRON = -1.6 * 1E-19

print(f"Formula is:\n{coulomb_law.print_law()}")

distance = symbols("dictance")

force_equation = coulomb_law.law.subs({
    coulomb_law.first_charge: CHARGE_OF_ELECTRON,
    coulomb_law.second_charge: CHARGE_OF_ELECTRON,
    coulomb_law.distance: distance
})

p1 = plot(force_equation, (distance, 0, 50),
          line_color="blue",
          title="Coluomb Law",
          xlabel="r, m",
          ylabel="F, N",
          backend=MatplotlibBackend,
          legend=True,
          show=False)
p1.show()
