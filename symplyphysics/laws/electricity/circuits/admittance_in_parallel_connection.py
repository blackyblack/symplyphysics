"""
Admittance in parallel connection
=================================

The total admittance of the circuit whose components are connected in parallel is the sum
of the admittances of individual components.

**Conditions:**

#. Components are connected in parallel.

**Links:**

#. `Electronics Tutorials, "Admittance of a Parallel RLC Circuit" <https://www.electronics-tutorials.ws/accircuits/parallel-circuit.html>`__.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (units, Quantity, validate_input, validate_output, global_index,
    IndexedSum, symbols)
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import admittance_is_inverse_impedance as _admittance_def
from symplyphysics.laws.electricity import current_is_voltage_over_impedance as _ohms_law
from symplyphysics.laws.electricity.circuits import (
    sum_of_currents_through_junction_is_zero as _kirchhoffs_law,)

total_admittance = symbols.admittance
"""
Total :symbols:`admittance` of the circuit.
"""

admittance = clone_as_indexed(symbols.admittance)
"""
:symbols:`admittance` of :math:`i`-th circuit.
"""

law = Eq(total_admittance, IndexedSum(admittance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Ohm's law and Kirchhoff laws for a two-component system. The expression for more
# than two components can be found in a similar way.

_current = _ohms_law.current

_admittance = _admittance_def.admittance

# NOTE: voltage is the same across all components of a parallel connection.
_admittance_eqn = _admittance_def.definition.subs({
    _admittance_def.impedance: _ohms_law.impedance,
})

_current_expr = solve(
    (_ohms_law.law, _admittance_eqn),
    (_current, _ohms_law.impedance),
    dict=True,
)[0][_current]

_first_current_expr = _current_expr.subs(_admittance, admittance[1])
_second_current_expr = _current_expr.subs(_admittance, admittance[2])

# This is the current if we were to replace the two components with a single component that yields
# the same current and voltage.
_total_current_expr = _current_expr.subs(_admittance, total_admittance)

_total_current_eqn = _kirchhoffs_law.law.subs(
    _kirchhoffs_law.index,
    Idx("i", (1, 3)),
).doit().subs({
    _kirchhoffs_law.current[1]: _first_current_expr,  # flowing in, positive sign
    _kirchhoffs_law.current[2]: _second_current_expr,  # flowing in, positive sign
    _kirchhoffs_law.current[3]: -1 * _total_current_expr,  # flowing out, negative sign
})

_total_admittance_derived = solve(_total_current_eqn, total_admittance)[0]

_total_admittance_expected = law.rhs.subs(global_index, Idx("i", (1, 2))).doit()

assert expr_equals(_total_admittance_derived, _total_admittance_expected)


@validate_input(admittances_=admittance)
@validate_output(units.conductance)
def calculate_parallel_admittance(admittances_: list[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(admittances_)))
    admittances_law = law.subs(global_index, local_index)
    admittances_law = admittances_law.doit()
    solved = solve(admittances_law, total_admittance, dict=True)[0][total_admittance]
    for i, v in enumerate(admittances_):
        solved = solved.subs(admittance[i + 1], v)
    return Quantity(solved)
