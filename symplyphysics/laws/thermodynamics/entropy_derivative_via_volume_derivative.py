r"""
Entropy derivative via volume derivative
========================================

Maxwell relations are a set of equations that unite the most common thermodynamic quantities between
one another. They are derived from the fundamental themodynamic relations featuring differentials
of thermodynamic potentials, and this method of derivation is called the method of thermodynamic
potentials.

**Conditions:**

#. Particle count must be constant.
"""

from sympy import (
    Eq,
    Derivative,
    symbols as sym_symbols,
    Function as SymFunction,
    solve,
)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_differential

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

entropy = clone_as_function(symbols.entropy, [temperature, pressure])
"""
:symbols:`entropy` of the system as a function of temperature and pressure.
"""

volume = clone_as_function(symbols.volume, [temperature, pressure])
"""
:symbols:`volume` of the system as a function of temperature and pressure.
"""

law = Eq(Derivative(entropy(temperature, pressure), pressure),
    -1 * Derivative(volume(temperature, pressure), temperature))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from expression of Gibbs energy differential

_gibbs_energy_expr = sym_symbols("gibbs_energy", cls=SymFunction)(temperature, pressure)

# Since Gibbs energy is a differentiable function of all its parameters, one can apply
# the symmetry of second derivatives, which is that for any sufficiently differentiable
# function one can exchange the order of taking partial derivatives: d(df/dx)/dy = d(df/dy)/dx.
# The expression for the second derivative of Gibbs energy must be the same whether one first
# takes it w.r.t. temperature and then w.r.t. pressure or the other way round.

_gibbs_energy_diff_temperature = gibbs_energy_differential.law.rhs.subs({
    gibbs_energy_differential.temperature_change: 1,
    gibbs_energy_differential.pressure_change: 0,
    gibbs_energy_differential.particle_count_change: 0,
    gibbs_energy_differential.entropy: entropy(temperature, pressure)
})

_gibbs_energy_diff_temperature_diff_pressure = _gibbs_energy_diff_temperature.diff(pressure)

_eqn_temperature_pressure = Eq(_gibbs_energy_expr.diff(temperature, pressure),
    _gibbs_energy_diff_temperature_diff_pressure)

_gibbs_energy_diff_pressure = gibbs_energy_differential.law.rhs.subs({
    gibbs_energy_differential.temperature_change: 0,
    gibbs_energy_differential.pressure_change: 1,
    gibbs_energy_differential.particle_count_change: 0,
    gibbs_energy_differential.volume: volume(temperature, pressure)
})

_gibbs_energy_diff_pressure_diff_temperature = _gibbs_energy_diff_pressure.diff(temperature)

_eqn_pressure_temperature = Eq(
    _gibbs_energy_expr.diff(pressure, temperature),
    _gibbs_energy_diff_pressure_diff_temperature,
)

_entropy_diff_pressure = solve(
    (_eqn_temperature_pressure, _eqn_pressure_temperature),
    (entropy(temperature, pressure).diff(pressure), _gibbs_energy_expr.diff(temperature, pressure)),
    dict=True,
)[0][entropy(temperature, pressure).diff(pressure)]

assert expr_equals(_entropy_diff_pressure, law.rhs)


@validate_input(
    volume_change_=volume,
    temperature_change_=temperature,
)
@validate_output(units.energy / (units.temperature * units.pressure))
def calculate_entropy_differential(
    volume_change_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    volume_ = (volume_change_ / temperature_change_) * temperature
    result = law.rhs.subs(volume(temperature, pressure), volume_).doit()
    return Quantity(result)
