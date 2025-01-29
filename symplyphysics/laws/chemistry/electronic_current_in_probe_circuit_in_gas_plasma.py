"""
Electron current in probe circuit in gas plasma
===============================================

A gas-discharge plasma is an ionized gas in which the concentrations of positively and
negatively charged particles are approximately equal to each other, and the Debye
shielding radius is significantly smaller than the characteristic size of the volume in
which the ionized gas is located.

A probe is an auxiliary metal electrode that is inserted into a plasma volume for its
examination. Flat, cylindrical and spherical probesare used. A potential is set relative
to one of the electrodes on the probe and the dependence of the incoming current on the
value of this potential is removed.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`electron_rest_mass`.
#. :quantity_notation:`boltzmann_constant`.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, exp, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import (
    elementary_charge,
    electron_rest_mass,
    boltzmann_constant,
)

probe_current = symbols.current
"""
Probe :symbols:`current`.
"""

probe_surface_area = symbols.area
"""
Probe surface :symbols:`area`.
"""

electron_concentration = symbols.number_density
"""
Electron concentration (:symbols:`number_density`) in plasma.
"""

plasma_temperature = symbols.temperature
"""
Plasma :symbols:`temperature`.
"""

floating_plasma_potential = clone_as_symbol(symbols.electric_potential,
    display_symbol="U_f",
    display_latex="U_\\mathbf{f}")
"""
Floating plasma :symbols:`electric_potential`.
"""

probe_potential = symbols.electric_potential
"""
:symbols:`electric_potential` at the location of the probe.
"""

law = Eq(probe_current, (0.25 * probe_surface_area * elementary_charge * electron_concentration) *
    sqrt(8 * boltzmann_constant * plasma_temperature /
    (pi * electron_rest_mass)) * exp(-elementary_charge *
    (floating_plasma_potential - probe_potential) / (boltzmann_constant * plasma_temperature)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(area_probe_surface_=probe_surface_area,
    electron_concentration_=electron_concentration,
    plasma_temperature_=plasma_temperature,
    floating_plasma_potential_=floating_plasma_potential,
    probe_potential_=probe_potential)
@validate_output(probe_current)
def calculate_current(area_probe_surface_: Quantity, electron_concentration_: Quantity,
    plasma_temperature_: Quantity, floating_plasma_potential_: Quantity,
    probe_potential_: Quantity) -> Quantity:
    result_expr = solve(law, probe_current, dict=True)[0][probe_current]
    result_expr = result_expr.subs({
        probe_surface_area: area_probe_surface_,
        electron_concentration: electron_concentration_,
        plasma_temperature: plasma_temperature_,
        floating_plasma_potential: floating_plasma_potential_,
        probe_potential: probe_potential_,
    })
    return Quantity(result_expr)
