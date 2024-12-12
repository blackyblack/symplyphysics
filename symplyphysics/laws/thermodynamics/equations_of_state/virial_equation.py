r"""
Virial equation
===============

Also called the *virial expansion*, the *virial equation of state* expresses the compressibility factor
(and therefore the pressure) of a real gas in local equilibrium as a power series of molar density.

**Notes:**

#. The first virial coefficient :math:`A` is defined to be 1 in order to enforce that the equation
   reduces to the ideal gas equation as gas density approaches zero.
#. The :math:`n`-th virial coefficient represents non-additive :math:`n`-body interactions of
   particles and all mutual interactions of :math:`2` up to :math:`(n - 1)` particles.
#. In general, virial coefficients are functions of temperature.
#. :math:`O(\dots)` is the mathematical *Big O*.

**Conditions:**

#. Interactions between 4 and more bodies are quite rare to happen, so the expansion is truncated to contain only
   the second and third virial coefficients. Moreover, the latter have been extensively studied and tabulated
   for many fluids.
#. In this law the limit :math:`\rho \to 0` is assumed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Virial_expansion>`__.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

compressibility_factor = Symbol("compressibility_factor", dimensionless, positive=True)
"""
See :doc:`definitions.compressibility_factor_is_deviation_from_ideal_gas`.

Symbol:
    :code:`Z`
"""

second_virial_coefficient = Symbol("second_virial_coefficient",
    units.volume / units.amount_of_substance,
    real=True)
"""
Second virial coefficient correponding to pair interactions between particles.

Symbol:
    :code:`B`
"""

third_virial_coefficient = Symbol("third_virial_coefficient",
    (units.volume / units.amount_of_substance)**2,
    real=True)
"""
Third virial coefficient corresponding to 3-body interaction between particles.

Symbol:
    :code:`C`
"""

molar_density = Symbol("molar_density", units.amount_of_substance / units.volume, positive=True)
r"""
Molar density of the system, or as amount of substance per unit volume.
See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

law = Eq(
    compressibility_factor,
    1 + second_virial_coefficient * molar_density + third_virial_coefficient * molar_density**2)
r"""
:code:`Z = 1 + B * rho + C * rho^2 + O(rho^3)`

Latex:
    .. math::
        Z = 1 + B \rho + C \rho^2 + O \! \left( \rho^3 \right)
"""


@validate_input(
    second_virial_coefficient_=second_virial_coefficient,
    third_virial_coefficient_=third_virial_coefficient,
    molar_density_=molar_density,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    second_virial_coefficient_: Quantity,
    third_virial_coefficient_: Quantity,
    molar_density_: Quantity,
) -> float:
    result = law.rhs.subs({
        second_virial_coefficient: second_virial_coefficient_,
        third_virial_coefficient: third_virial_coefficient_,
        molar_density: molar_density_,
    })

    result_value = convert_to_float(result)
    if result_value < 0:
        raise ValueError(f"Compressibility factor cannot be negative, got {result_value} instead")
    return result_value
