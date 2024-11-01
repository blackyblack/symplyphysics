r"""
Phase of traveling wave
=======================

The phase of a wave or other periodic function of some real variable :math:`t` is an angle-like
quantity representing the fraction of the cycle covered up to :math:`t`. As the variable :math:`t`
completes a full period, the phase increases by :math:`360^\circ` or :math:`2 \pi`.

If a function :math:`h(x, t)` describes a traveling wave, then position :math:`x` and time :math:`t`
can only appear in the form of the wave phase described below.

**Notes:**

#. :math:`\omega = (\vec \omega \cdot {\vec e}_x)`, i.e. the angular frequency is a positive quantity
   if the wave travels in the positive direction of the :math:`x`-axis. Here :math:`{\vec e}_x` is the
   unit vector pointing in the positive direction of the :math:`x`-axis.

**Conditions:**

#. This law applies to a 1-dimensional traveling wave.
#. The constant phase shift is not taken into account.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

wave_phase = symbols.phase
"""
:symbols:`phase` of the wave.
"""

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of the wave.
"""

position = symbols.position
"""
:symbols:`position`, or spatial coordinate.
"""

angular_frequency = symbols.angular_frequency
r"""
:symbols:`angular_frequency` of the wave.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(wave_phase, angular_wavenumber * position - angular_frequency * time)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_wavenumber_=angular_wavenumber,
    position_=position,
    angular_frequency_=angular_frequency,
    time_=time,
)
@validate_output(wave_phase)
def calculate_wave_phase(
    angular_wavenumber_: Quantity,
    position_: Quantity,
    angular_frequency_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        angular_wavenumber: angular_wavenumber_,
        position: position_,
        angular_frequency: angular_frequency_,
        time: time_,
    })
    return Quantity(result)
