"""
Instantaneous energy of resonator
=================================

A rectangular resonator consists of metal walls and a material filling it. The resonator
is capable of storing energy. The instantaneous value of the resonator energy depends on
its quality factor, initial energy value, time and angular frequency.

..
    TODO: find link
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

instantaneous_energy = symbols.energy
"""
Instantaneous :symbols:`energy` of the resonator.
"""

initial_energy = clone_as_symbol(symbols.energy, subscript="0")
"""
Initial :symbols:`energy` of the resonator.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~instantaneous_energy` is measured.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the current.
"""

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the resonator.
"""

law = Eq(instantaneous_energy, initial_energy * exp(-angular_frequency * time / quality_factor))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(initial_energy_=initial_energy,
    time_=time,
    angular_frequency_=angular_frequency,
    quality_factor_=quality_factor)
@validate_output(instantaneous_energy)
def calculate_instantaneous_energy(initial_energy_: Quantity, time_: Quantity,
    angular_frequency_: Quantity, quality_factor_: float) -> Quantity:
    result_velocity_expr = solve(law, instantaneous_energy, dict=True)[0][instantaneous_energy]
    result_expr = result_velocity_expr.subs({
        initial_energy: initial_energy_,
        time: time_,
        angular_frequency: angular_frequency_,
        quality_factor: quality_factor_
    })
    return Quantity(result_expr)
