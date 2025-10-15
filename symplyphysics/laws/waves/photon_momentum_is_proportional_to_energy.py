"""
Photon momentum is proportional to energy
=========================================

The momentum of a photon is its energy divided by the speed of light.

See :doc:`laws.waves.photon_energy_is_proportional_to_angular_frequency` or
:doc:`laws.waves.photon_energy_is_proportional_to_linear_frequency` for energy
expressions.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Conditions:**

#. Works in a **vacuum** environment.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Photon#Energy_and_momentum>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import (
    angular_wavenumber_is_inverse_wavelength as _wavenumber_def,
    period_from_angular_frequency as _period_def,
)
from symplyphysics.laws.waves import (
    photon_momentum_is_proportional_to_angular_wavenumber as _momentum_law,
    photon_energy_is_proportional_to_angular_frequency as _energy_law,
    wavelength_from_phase_speed_and_period as _wavelength_law,
)

momentum = symbols.momentum
"""
:symbols:`momentum` of a photon.
"""

energy = symbols.energy
"""
:symbols:`energy` of a photon.
"""

law = Eq(momentum, energy / quantities.speed_of_light)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_momentum_eqn = _momentum_law.law.subs({
    _momentum_law.momentum: momentum,
})

_energy_eqn = _energy_law.law.subs({
    _energy_law.energy: energy,
})

_wavenumber_eqn = _wavenumber_def.definition.subs({
    _wavenumber_def.angular_wavenumber: _momentum_law.angular_wavenumber,
})

_wavelength_eqn = _wavelength_law.law.subs({
    _wavelength_law.wavelength: _wavenumber_def.wavelength,
    _wavelength_law.phase_velocity: quantities.speed_of_light,  # only holds in vacuum
})

_period_eqn = _period_def.law.subs({
    _period_def.period: _wavelength_law.period,
    _period_def.angular_frequency: _energy_law.angular_frequency,
})

_momentum_derived = solve(
    (_momentum_eqn, _energy_eqn, _wavenumber_eqn, _wavelength_eqn, _period_eqn),
    (momentum, _momentum_law.angular_wavenumber, _energy_law.angular_frequency,
    _wavenumber_def.wavelength, _wavelength_law.period),
    dict=True,
)[0][momentum]

assert expr_equals(_momentum_derived, law.rhs)


@validate_input(photon_energy_=energy)
@validate_output(momentum)
def calculate_momentum(photon_energy_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_momentum_expr.subs({energy: photon_energy_})
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 414
