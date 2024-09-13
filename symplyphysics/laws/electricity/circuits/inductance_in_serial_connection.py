"""
Inductance in serial connection
===============================

The total inductance of the circuit whose components are connected in series is the sum
of the inductances of individual components.

**Conditions:**

#. Components are connected in series.
#. Inductors are not magnetically coupled.
"""

from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    inductance_is_magnetic_flux_over_current as _inductance_law
)

total_inductance = Symbol("total_inductance", units.inductance)
"""
Total inductance of the circuit.

Symbol:
    :code:`L`
"""

inductance = SymbolIndexed("inductance", units.inductance)
r"""
Inductance of the :math:`i`-th component.

Symbol:
    :code:`L_i`

Latex:
    :math:`L_i`
"""

law = Eq(total_inductance, SumIndexed(inductance[global_index], global_index))
r"""
:code:`L = Sum(L_i, i)`

Latex:
    .. math::
        L = \sum_i L_i
"""

# Derive law for two inductors. The current is the same in components connected
# in parallel.

_magnetic_flux_expr = solve(_inductance_law.law, _inductance_law.magnetic_flux)[0]

_magnetic_flux_1 = _magnetic_flux_expr.subs(_inductance_law.inductance, inductance[1])

_magnetic_flux_2 = _magnetic_flux_expr.subs(_inductance_law.inductance, inductance[2])

_total_magnetic_flux = _magnetic_flux_1 + _magnetic_flux_2

_total_inductance = _inductance_law.law.rhs.subs(
    _inductance_law.magnetic_flux, _total_magnetic_flux
)

_local_idx = Idx("local_idx", (1, 2))
_inductance_from_law = law.rhs.subs(global_index, _local_idx).doit()

assert expr_equals(_total_inductance, _inductance_from_law)


@validate_input(inductances_=inductance)
@validate_output(total_inductance)
def calculate_serial_inductance(inductances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(inductances_)))
    inductances_law = law.subs(global_index, local_index)
    inductances_law = inductances_law.doit()
    solved = solve(inductances_law, total_inductance, dict=True)[0][total_inductance]
    for i, v in enumerate(inductances_):
        solved = solved.subs(inductance[i + 1], v)
    return Quantity(solved)
