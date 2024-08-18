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
   if the wave travels in the positive direction of the :math:`x`-axis.

**Conditions:**

#. This law applies to a :math:`1`-dimensional traveling wave.
#. The constant phase shift is not taken into account.
"""

from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

wave_phase = Symbol("wave_phase", angle_type, real=True)
r"""
Phase of the wave.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
"""
Angular wavenumber of the wave.

Symbol:
    :code:`k`
"""

position = Symbol("position", units.length, real=True)
"""
Position, or spatial coordinate.

Symbol:
    :code:`x`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time, real=True)
r"""
Angular frequency of the wave.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

time = Symbol("time", units.time, positive=True)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(wave_phase, angular_wavenumber * position - angular_frequency * time)
r"""
:code:`phi = k * x - w * t`

Latex:
    .. math::
        \varphi = k x - \omega t
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
