r"""
Average kinetic energy of ideal gas from temperature
====================================================

The kinetic theory of ideal gases allows us to determine the *average kinetic energy* for an ideal gas.
It states that the average kinetic energy for all ideal gases is directly proportional to the absolute
temperature of the gas and only depends on the temperature.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is ideal.
"""

from sympy import Eq, solve, Rational, stats, Interval, S
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.operations.average import Average
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import energy_distribution

average_kinetic_energy = Average(symbols.kinetic_energy, positive=True)
"""
Average :symbols:`kinetic_energy` of the ideal gas.
"""

equilibrium_temperature = clone_as_symbol(
    symbols.temperature,
    positive=True,
)
"""
Equilibrium :symbols:`temperature` of the gas.
"""

law = Eq(average_kinetic_energy,
    Rational(3, 2) * quantities.boltzmann_constant * equilibrium_temperature)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Maxwell-Boltzmann energy distribution

_distribution = energy_distribution.law.rhs.subs(
    energy_distribution.equilibrium_temperature,
    equilibrium_temperature,
)

_random_energy_variable = stats.ContinuousRV(
    energy_distribution.energy,
    _distribution,
    set=Interval(0, S.Infinity),
)

_average_energy_derived = stats.E(_random_energy_variable)

assert expr_equals(_average_energy_derived, law.rhs)


@validate_input(temperature_=symbols.temperature)
@validate_output(average_kinetic_energy)
def calculate_average_kinetic_energy(temperature_: Quantity) -> Quantity:
    result_expr = solve(law, average_kinetic_energy, dict=True)[0][average_kinetic_energy]
    result_average_kinetic_energy = result_expr.subs(equilibrium_temperature, temperature_)
    return Quantity(result_average_kinetic_energy)
