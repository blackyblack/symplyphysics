r"""
Energy via constant power and time
==================================

When power is constant, energy can be expressed as a product of power and time.

**Conditions:**

#. Power is constant, i.e. :math:`\frac{d P}{d t} = 0`.
"""

from sympy import Eq, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import power_is_energy_derivative as power_def

energy = Symbol("energy", units.energy)
"""
Energy consumed or released during time :math:`t`.

Symbol:
    :code:`E`
"""

power = Symbol("power", units.power)
"""
Power.

Symbol:
    :code:`P`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(energy, power * time)
r"""
:code:`E = P * t`

Latex:
    .. math::
        E = P t
"""

# Derive from definition of power from energy

_power_eqn = (
    power_def.definition
    .replace(power_def.power, lambda _: power)
    .subs(power_def.time, time)
)

_energy_expr = dsolve(
    _power_eqn,
    power_def.energy(time),
    ics={power_def.energy(0): 0},
).rhs

assert expr_equals(_energy_expr, law.rhs)


@validate_input(
    power_=power,
    time_=time,
)
@validate_output(energy)
def calculate_energy(
    power_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        power: power_,
        time: time_,
    })
    return Quantity(result)
