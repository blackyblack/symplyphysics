"""
Transmission coefficient approximation of low-pass filter
=========================================================

The approximation of the power transmission coefficient of a normalized low-pass filter
is given by approximating functions of the order of :math:`n`. 

"""

from sympy import Eq, solve
from symplyphysics import (
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
)

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the signal.
"""

filter_function = Symbol("F", dimensionless)
"""
Function (of :attr:`~frequency`) of order :math:`n` which approximates the transfer
coefficient.
"""

bandwidth_distortion = Symbol("e", dimensionless)
"""
Bandwidth distortion determines the maximum distortion in the bandwidth.
"""

transfer_coefficient = Symbol("H", dimensionless)
"""
Transfer coefficient of the filter.
"""

law = Eq(transfer_coefficient, 1 / (1 + bandwidth_distortion**2 * filter_function**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(filter_function_at_frequency=filter_function,
    bandwidth_distortion_=bandwidth_distortion)
@validate_output(transfer_coefficient)
def calculate_coefficient(filter_function_at_frequency: float,
    bandwidth_distortion_: float) -> float:
    result_expr = law.subs({
        filter_function: filter_function_at_frequency,
        bandwidth_distortion: bandwidth_distortion_,
    })
    result = solve(result_expr, transfer_coefficient, dict=True)[0][transfer_coefficient]
    return convert_to_float(result)
